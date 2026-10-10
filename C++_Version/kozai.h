#include <cstdio>
#include <cstdlib>      // stdlib.h with namespace std
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <cmath>        // math.h   with namespace std
#include <cstring>
#include <string>
//#include <fstream>

using namespace std;

#include "vector.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>

#define sqr(x) ((x)*(x))
#define DIMENSION 13

#define PRI(x) {for (int __pri__ = 0; __pri__ < x; __pri__++) cerr << " ";}
#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << endl

#define isnan(x) ((x) != (x))

//G's and c's (with some additional powers)
//All in CGS, just to annoy the relativists
//
//MODIFIED (Arakida+): these are no longer `const'.  With --units code the whole
//set is switched to the G = 1 system used in our Python single-/double-averaged
//runs (masses in Msun-like code units, lengths in a0-like code units, time in
//code units), and the speed of light becomes a free parameter (--clight).
//Call set_units() BEFORE any set_*() that multiplies by MSUN/AU/RSUN/YEAR.
double G=6.67407e-8; // cm^3 / g / s^2
double c=29979245800.; // cm / s
double Gc2=G/c/c; // cm / g
double G3=G*G*G;
double c5=c*c*c*c*c;
double G3c5 = G*G*G / (c*c*c*c*c);

//Conversion Factors
double YEAR=3.15569e7; // s
double AU=14959787070000.0; // cm
double RSUN=69570000000.0; // cm
double MSUN=1.9884754153381438e+33; // g
const double DEG=0.017453292519943295;

//Recompute everything that is built out of G and c
inline void recompute_derived(){
  Gc2  = G/c/c;
  G3   = G*G*G;
  c5   = c*c*c*c*c;
  G3c5 = G3/c5;
}

// unit_mode 0 = CGS/astro (original), 1 = code units (G = 1, Msun = AU = yr = 1)
inline void set_units(int unit_mode, double clight){
  if(unit_mode == 1){
    G = 1.0;
    YEAR = 1.0; AU = 1.0; MSUN = 1.0; RSUN = 1.0;
    // free parameter; a large value is the Newtonian limit
    c = (clight > 0.) ? clight : 1.0e6;
  }
  else if(clight > 0.) c = clight;
  recompute_derived();
}

//Some pN constants
const double c304o15=20.266666666666666;
const double c64o5=12.8;
const double c121o304=0.3980263157894737;
const double c73o24=3.0416666666666665;
const double c37o96=0.3854166666666667;

const double c29o19 = 29./19.;
const double c1181o2299 = 1181./2299.;
const double c12o19 = 12./19.;
const double c870o2299= 870./2299.;

const double SMALLNUMBER=1e-10;
const double TWOPI=6.28318530718;
const double PI=3.14159265359;

class kozai_struct;

int rhs(double t, const double y[], double f[], void *kozai_ptr);
void set_parameters(int argc, char **argv, kozai_struct *kozai, double& t_min, double& t_step, bool& IGNORE_GSL_ERROR,bool& PRINT_GSL);
void print_state(double t_next, kozai_struct *kozai, ostream& stream);
void print_header_and_initial_state(kozai_struct* kozai);
double peters_t(kozai_struct *kozai);
double peters_integral(double e, void *params);

class kozai_struct
{
	private:

		//triple configuration (initial)
		//
		//these are *not* returned by the getter
		//functions below; instead, they are
		//computed from the evolved y vector
		double a1;
		double a2;
		double e1;
		double e2;
		double inc; 
		double g1;
		double g2;
		double m1;
		double m2;
		double m3;
		double Omega1;
		double Omega2;
		double r1;
		double r2;

		//Things to compute ahead of time...
		double L1_no_a;
		double L2_no_a;
		double mu1;
		double mu2;
		double m1m2;
		double m2m1;
		double eta;

		//flags determining what terms to include
		bool quadrupole;
		bool octupole;
		bool pericenter;
    bool outerpericenter;
    bool crossterms;
		bool radiation;

    //---------------------------------------------------------------
    // MODIFIED (Arakida+): fine-grained control of the Lim & Rodriguez
    // (2020) cross terms, plus the two terms our analysis adds/corrects.
    //
    //  x_rn4   : grade (1,4)   Rn4   DN6+ID6+IDPNB + outer-period piece
    //            -> passes every canonicity test in our Paper I; keep on.
    //  x_rn2p5 : grade (1,5/2) Rn2p5 DN6 = de Sitter, LR Eq. (4.4); keep on.
    //  x_rn1p5 : grade (1,3/2) Rn1p5 IDNQ = LR Eq. (4.5) "libration" term.
    //            Our Paper I shows K^(1,3/2) == 0, i.e. this term is spurious.
    //  x_rn3p5 : grade (1,7/2) Rn3p5, absent from this C++ port.
    //            1 = our/KST value  K^(1,7/2) = -K_star
    //            2 = LR (A7)-(A11) reconstructed Hamiltonian (+K_star-prefactor/8,
    //                with LR's angular structure)
    //  lr_rn1p5_exact / lr_rn2p5_exact : use LR's PUBLISHED mass prefactors
    //            instead of the ones hard-wired in this port (see README.md).
    //  brown   : grade (0,9/2) omega_2-averaged Brown Hamiltonian C eps21 F21.
    //  restricted : freeze the outer orbit (j2, e2 = const) -> the restricted
    //            problem our Python runs integrate.
    //---------------------------------------------------------------
    bool x_rn4, x_rn2p5, x_rn1p5;
    int  x_rn3p5;
    bool lr_rn1p5_exact, lr_rn2p5_exact;
    bool brown;
    bool restricted;
    bool outer_frame;   // reference plane = outer orbit plane instead of invariant plane

		//The initial configuration of the system in j,e,s space
		//and the initial 18 dimensional vector (allocated outside class)
		vec j1_init;
		vec j2_init;
		vec e1_init;
		vec e2_init;

		double *y;


	public:

		//Sloppy af constructor
		kozai_struct(double ia1=365.11437, double ia2=2865.7776, double ie1=0.32229033, double ie2=0.35653592, double iinc=93.920405, 
			  double ig1=199.31831, double ig2=237.15944, double im1=24.22645, double im2=14.999986, double im3=21.509239,
			  double iOmega1=321.97666, double iOmega2=141.97666, double ir1=0, double ir2=0,
			  bool iquad=false, bool ioct=false, bool iperi=false, bool ioperi=false, bool icross=false,bool rad=false)
			{	
				a1=ia1*AU;
				a2=ia2*AU;
				e1=ie1;
				e2=ie2;
				inc=iinc*DEG;
				g1=ig1*DEG;
				g2=ig2*DEG;
				m1=im1*MSUN;
				m2=im2*MSUN;
				m3=im3*MSUN;
				Omega1=iOmega1*DEG;
				Omega2=iOmega2*DEG;
				r1=ir1*RSUN;
				r2=ir2*RSUN;
				quadrupole=iquad;
				octupole=ioct;
				pericenter=iperi;
        outerpericenter=ioperi;
        crossterms=icross;
				radiation=rad;
				//MODIFIED: defaults for the new switches.  With --cross_lim alone the
				//code reproduces the ORIGINAL behaviour exactly (all three LR blocks on,
				//port prefactors, no (1,7/2), no Brown).
				x_rn4=true; x_rn2p5=true; x_rn1p5=true; x_rn3p5=0;
				lr_rn1p5_exact=false; lr_rn2p5_exact=false;
				brown=false; restricted=false; outer_frame=false;
				y = new double[DIMENSION];
			}

		//destructor
		~kozai_struct() {delete y;}

		void initialize()
		{
			mu1 = m1*m2 / (m1+m2);
			mu2 = (m1+m2)*m3 / (m1+m2+m3);
			L1_no_a = mu1*sqrt(G*(m1+m2));
			L2_no_a = mu2*sqrt(G*(m1+m2+m3));
			m1m2 = m1/m2;
			m2m1 = m2/m1;
			eta = mu1 / (m1+m2);

			//If no radius specified, set to the Schwarzschild radius.
			//MODIFIED: the original line was `2*1476.*m1/MSUN', i.e. 2GM/c^2 with the
			//1476 written in METRES while every other length here is in cm (a factor
			//100 too small), and it is not unit-agnostic.  Written covariantly:
			if(r1 == 0.) r1 = 2.*G*m1/(c*c);
			if(r2 == 0.) r2 = 2.*G*m2/(c*c);

			//Using the angular momenta, figure out what the individual
			//inclinations are wrt the invariant plane
			//
			//NOTE: this means we're in a coordinate system where z-hat points
			//towards to total system angular momentum (before radiation
			//reaction)
			double ia = L1_no_a*sqrt(a1*(1-sqr(e1)));
			double oa = L2_no_a*sqrt(a2*(1-sqr(e2)));


			double total_ang = sqrt(sqr(ia)+sqr(oa)+2*ia*oa*cos(inc));
			double inc1 = acos((ia + oa*cos(inc))/total_ang);
			double inc2 = acos((oa + ia*cos(inc))/total_ang);

			//MODIFIED: --outer_frame puts z-hat along the OUTER orbit normal instead
			//of along the total angular momentum.  Then Omega_1 printed by the code is
			//the Omega_1 of our papers (node of the inner orbit on the outer plane),
			//and iota_1 == iota.  Use together with --restricted to reproduce the
			//restricted problem integrated by our Python runs.
			if(outer_frame){ inc1 = inc; inc2 = 0.; }

			//Define the actual vectors
			j1_init = sqrt(1 - sqr(e1)) * vec(sin(inc1) * sin(Omega1), -sin(inc1) * cos(Omega1), cos(inc1));
			j2_init = sqrt(1 - sqr(e2)) * vec(sin(inc2) * sin(Omega2), -sin(inc2) * cos(Omega2), cos(inc2));
			e1_init = e1 * vec(cos(g1)*cos(Omega1) - sin(g1)*cos(inc1)*sin(Omega1),
					cos(g1)*sin(Omega1) + sin(g1)*cos(inc1)*cos(Omega1),
					sin(g1)*sin(inc1));
			e2_init = e2 * vec(cos(g2)*cos(Omega2) - sin(g2)*cos(inc2)*sin(Omega2),
					cos(g2)*sin(Omega2) + sin(g2)*cos(inc2)*cos(Omega2),
					sin(g2)*sin(inc2));

			vec yhat = j1_init^e1_init;
			yhat /= abs(yhat);
			vec zhat = j1_init / abs(j1_init);

			//And set the components into the y vector (which will be passed to
			//GSL)
			for(int i=0; i<3; i++){
				y[i] = j1_init[i];
				y[i+3] = e1_init[i];
				y[i+6] = j2_init[i];
				y[i+9] = e2_init[i];
			}
			y[12] = a1;
		}

		//The getters and setters for the class;
		//better to use this than the constructor
		//
		//The inner binary
		void set_a1(double a1_i) {a1=a1_i;}
		double get_a1() {return y[12];}

		void set_ecc1(double e1_i) {e1=e1_i;}
		double get_ecc1() {return abs(this->get_e1());}

		void set_g1(double g1_i) {g1=g1_i;}
		double get_g1() {
			double g1 = asin(this->get_e1()[2] / (sin(this->get_inc1())) / this->get_ecc1());
			if((this->get_e1()^this->get_j1())[2] > 0) g1 = PI - g1;
			if(g1 < 0) g1 += TWOPI; // need to put the angle in the right quadrant
			return g1;
		}

		void set_Omega1(double Omega1_i) {Omega1=Omega1_i;}
		double get_Omega1() {
			double norm_xy = sqrt(sqr(this->get_j1()[0]) + sqr(this->get_j1()[1]));
			if(norm_xy == 0.)
				return 0.;
			else{
				if(this->get_j1()[0] >= 0)
					return acos(-this->get_j1()[1] / norm_xy);
				else
		        	return TWOPI - acos(-this->get_j1()[1] / norm_xy);
			}
		}

		void set_m1(double m1_i) {m1=m1_i;}
		double get_m1() {return m1;}

		void set_m2(double m2_i) {m2=m2_i;}
		double get_m2() {return m2;}

		//Outer Binary
		void set_a2(double a2_i) {a2=a2_i;}
		double get_a2() {return a2;}

		void set_ecc2(double e2_i) {e2=e2_i;}
		double get_ecc2() {return abs(this->get_e2());}

		void set_g2(double g2_i) {g2=g2_i;}
		double get_g2() {
			double g2 = asin(this->get_e2()[2] / (sin(this->get_inc2())) / this->get_ecc2());
			if((this->get_e2()^this->get_j2())[2] > 0) g2 = PI- g2;
			if(g2 < 0) g2 += TWOPI; // need to put the angle in the right quadrant
			return g2;
		}

		void set_Omega2(double Omega2_i) {Omega2=Omega2_i;}
		double get_Omega2() {
			double norm_xy = sqrt(sqr(this->get_j2()[0]) + sqr(this->get_j2()[1]));
			if(norm_xy == 0.)
				return 0.;
			else{
				if(this->get_j2()[0] >= 0)
					return acos(-this->get_j2()[1] / norm_xy);
				else
		        	return TWOPI - acos(-this->get_j2()[1] / norm_xy);
			}
		}

		void set_m3(double m3_i) {m3=m3_i;}
		double get_m3() {return m3;}

		void set_inc(double inc_i) {inc=inc_i;}
		double get_inc() {
			return acos(this->get_j1()*this->get_j2()/(abs(this->get_j1())*abs(this->get_j2())));
		}

		double get_inc1(){
			return acos(this->get_j1()[2] / abs(this->get_j1()));
		}

		double get_inc2(){
			return acos(this->get_j2()[2] / abs(this->get_j2()));
		}

		//Physical Properties of Inner Binary (including spins)
		void set_r1(double r1_i) {r1=r1_i*RSUN;}
		double get_r1() {return r1;}

		void set_r2(double r2_i) {r2=r2_i*RSUN;}
		double get_r2() {return r2;}

		double collision() {
			return (this->get_a1()*(1-this->get_ecc1()) < (this->get_r1() + this->get_r2()));
		}

		double get_forb() {
			return sqrt(G*(this->get_m1()+this->get_m2())/pow(this->get_a1(),3))/(2*M_PI);
		}

		// get the gravitational-wave frequency of the inner binary
		double gwave_freq() {
			double ecc = this->get_ecc1();
			double forb = this->get_forb();
			return 2*forb*pow(1+ecc,1.1954)/pow(1-ecc*ecc,1.5);
		}

		//Which terms to include in integration
		void set_quadrupole(bool quadrupole_i) {quadrupole=quadrupole_i;}
		bool get_quadrupole() {return quadrupole;}

		void set_octupole(bool octupole_i) {octupole=octupole_i;}
		bool get_octupole() {return octupole;}

		void set_pericenter(bool pericenter_i) {pericenter=pericenter_i;}
		bool get_pericenter() {return pericenter;}

        void set_outer_pericenter(bool opericenter_i) {outerpericenter=opericenter_i;}
        bool get_outer_pericenter() {return outerpericenter;}

        void set_cross_lim(bool cross_i) {crossterms=cross_i;}
        bool get_cross_lim() {return crossterms;}

		void set_radiation(bool radiation_i) {radiation=radiation_i;}
		bool get_radiation() {return radiation;}

		//----- MODIFIED (Arakida+): the new switches -----------------------
		void set_x_rn4(bool b)   {x_rn4=b;}    bool get_x_rn4()   {return x_rn4;}
		void set_x_rn2p5(bool b) {x_rn2p5=b;}  bool get_x_rn2p5() {return x_rn2p5;}
		void set_x_rn1p5(bool b) {x_rn1p5=b;}  bool get_x_rn1p5() {return x_rn1p5;}
		void set_x_rn3p5(int  i) {x_rn3p5=i;}  int  get_x_rn3p5() {return x_rn3p5;}
		void set_lr_rn1p5_exact(bool b){lr_rn1p5_exact=b;} bool get_lr_rn1p5_exact(){return lr_rn1p5_exact;}
		void set_lr_rn2p5_exact(bool b){lr_rn2p5_exact=b;} bool get_lr_rn2p5_exact(){return lr_rn2p5_exact;}
		void set_brown(bool b)   {brown=b;}    bool get_brown()   {return brown;}
		void set_restricted(bool b){restricted=b;} bool get_restricted(){return restricted;}
		void set_outer_frame(bool b){outer_frame=b;} bool get_outer_frame(){return outer_frame;}

		//Quantities used by the comparison output; all in the current unit system.
		double get_n0() {return sqrt(G*(m1+m2)/pow(this->get_a1(),3.));}
		double get_n3() {return sqrt(G*(m1+m2+m3)/pow(a2,3.));}
		// eps21 = (n3/n0)(m3/M3)(1 + 2 e3^2/3)/eta3^3   (Paper II, Eq. for eps21)
		double get_eps21() {
			double e2n = this->get_ecc2();
			return (this->get_n3()/this->get_n0())*(m3/(m1+m2+m3))
			       *(1. + 2.*sqr(e2n)/3.)/pow(sqrt(1.-sqr(e2n)),3.);
		}
		// jz = j1 . n2  (the z-component of the inner j in the OUTER orbit frame)
		double get_jz() {
			vec j2v = this->get_j2();
			return (this->get_j1()*j2v)/abs(j2v);
		}
		// Omega_1 via atan2.  Mathematically identical to get_Omega1() above, but
		// that one uses acos(), which loses half the digits whenever Omega_1 is
		// near 0 or pi.  Used by --raw and --audit; get_Omega1() is left untouched
		// so that the ordinary output is bit-for-bit the original code's.
		double get_Omega1_acc() {
			vec j1v = this->get_j1();
			double a = atan2(j1v[0], -j1v[1]);
			return (a < 0.) ? a + TWOPI : a;
		}

		// omega_1 measured from the mutual node, and Omega_1 about the outer normal
		double get_omega1_mutual() {
			vec n2v = this->get_j2()/abs(this->get_j2());
			vec n1v = this->get_j1()/abs(this->get_j1());
			vec node = n2v^n1v;            // points to the ascending node of 1 on plane 2
			double nn = abs(node);
			if(nn < 1e-14) return 0.;
			node /= nn;
			vec u1v = this->get_e1()/this->get_ecc1();
			double cw = u1v*node;
			double sw = ((node^u1v)*n1v);
			double w = atan2(sw,cw);
			return (w < 0.) ? w + TWOPI : w;
		}

		//precomputed constants
		//
		//except for L1, which changes as a1 shrinks
		//due to radiation reaction
		//MODIFIED: was `sqrt(y[18])' -- an out-of-bounds read left over from an
		//18-dimensional version of the state vector that included the spins.
		//With DIMENSION 13 the semi-major axis lives in y[12].
		double get_L1() {return L1_no_a*sqrt(y[12]);}
		double get_L1_no_a() {return L1_no_a;}
		double get_L2_no_a() {return L2_no_a;}
		double get_mu1() {return mu1;}
		double get_mu2() {return mu2;}
		double get_m1m2() {return m1m2;}
		double get_m2m1() {return m2m1;}
		double get_eta() {return eta;}

		//Pull the vectors back out of y
		vec get_j1() {return vec(y[0],y[1],y[2]);}
		vec get_e1() {return vec(y[3],y[4],y[5]);}
		vec get_j2() {return vec(y[6],y[7],y[8]);}
		vec get_e2() {return vec(y[9],y[10],y[11]);}

		//Return the pointer to y
		double *get_y() {return y;}

};

