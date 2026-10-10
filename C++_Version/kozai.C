#include "kozai.h"

const string help="\nCarl's Kozai Code\n\n\
For integrating the secular equations of motions\n\
for the three-body problem, including the quadrupole\n\
and octupole terms (see Tremaine+2009 and Liu+2015)\n\n\
Also includes relativistic terms for the inner binary, including:\n\
  Pericenter precession (ala Ford+2000)\n\
  Gravitational Wave emission (ala Blaes+2002)\n\
  Spin-Orbit coupling (ala me, following Barker and O'Connell 1975)\n\
  Spin-Spin coupling (same)\n\n\
Also includes cross terms from relativistic three-body effects (Lim and Rodriguez 2020)\n\
Options:\n\
  --no-int  don't integrate the system, just print info\n\
  --help    print this\n\n\
Parameters:\n\
  --m1      mass 1 of inner binary (MSUN)\n\
  --m2      mass 2 of inner binary (MSUN)\n\
  --m3      mass of tertiary (MSUN)\n\
  --a1      Semi-major axis of inner binary (AU)\n\
  --a2      Semi-major axis of outer binary (AU)\n\
  --e1      Eccentricity of inner binary\n\
  --e2      Eccentricity of outer binary\n\
  --g1      Argument of pericenter for inner binary (deg)\n\
  --g2      Argument of pericenter for outer binary (deg)\n\
  --omega1  Longitude of ascending node for inner binary (deg)\n\
  --omega2  Longitude of ascending node for outer binary (deg)\n\
  --inc     Mutual Inclination (deg)\n\
  --rad1    Radius of mass 1 (RSUN)\n\
  --rad2    Radius of mass 2 (RSUN)\n\
  --time    Time to integrate (years)\n\
  --dt      How often to print output\n\
                Positive dt outputs every dt years\n\
                Negative dt outputs every 10^{-dt} adaptive timesteps\n\
                    (e.g. dt=-2 outputs every 100 timesteps)\n\n\
Flags:\n\
  --quad       Include quadrupole terms \n\
  --oct        Include octupole Terms\n\
  --peri       Include inner pericenter precession (1pN) terms\n\
  --outerperi  Include outer pericenter precession (1pN) terms\n\
  --cross_lim  Include cross terms\n\
  --ignore_gsl Ignore errors from GSL (i.e. when a binary merges, the inner\n\
               binary can decouple from the outer s.t. the optimal timestep\n\
               for the inner is beyond machine-tolerance away from the outer)\n\
  --print_gsl  Print status to cour and state to cerr when GSL_ERROR flagged\n\
\n\
MODIFIED BUILD (Arakida+) -- options for comparing against our results.\n\
See README.md for the mapping onto the equations of the two papers.\n\
\n\
Units and output:\n\
  --units code   G = 1, Msun = AU = yr = 1 (the system of our Python runs).\n\
                 MUST be given before/independently of the other options.\n\
  --clight X     speed of light in the current unit system (use a large value\n\
                 for the Newtonian limit).  Default in code units: 1e6.\n\
  --raw          machine-precision dump of the full state (see header line).\n\
  --outer_frame  reference plane = outer orbit plane, not the invariant plane,\n\
                 so the printed Omega_1 is the Omega_1 of the papers.\n\
  --restricted   freeze the outer orbit (the restricted problem we integrate).\n\
\n\
Individual cross-term grades (all on by default under --cross_lim):\n\
  --no_rn4       drop grade (1,4)    [LR Eq. (4.6) etc.; canonical, keep it]\n\
  --no_rn2p5     drop grade (1,5/2)  [LR Eq. (4.4), de Sitter; canonical]\n\
  --no_rn1p5     drop grade (1,3/2)  [LR Eq. (4.5); our result: K == 0]\n\
  --lr_rn1p5     use LR's published mass prefactor for (1,3/2)  (this port\n\
                 has m*sqrt(M) where LR Eq. (4.5) has m*m3/sqrt(M))\n\
  --lr_rn2p5     use LR's published mass prefactor for (1,5/2)  (this port\n\
                 has the test-particle limit 3 M^{3/2}/2)\n\
  --x_rn3p5 kst  add grade (1,7/2) with our/KST value K = -K_star\n\
  --x_rn3p5 lr   add grade (1,7/2) with LR (A7)-(A11) reconstructed as a\n\
                 Hamiltonian flow (same prefactor, LR angular structure)\n\
  --x_rn3p5 lrraw  LR (A7),(A9),(A11) transcribed verbatim (non-canonical;\n\
                 implies --cross_lim)\n\
  --brown        add the (0,9/2) Brown Hamiltonian C eps21 F21 (Newtonian)\n\
  --audit        evaluate the rhs once at t=0, one grade at a time, and print\n\
                 Gdot, iotadot, Omegadot_1 and the normalisation-free ratios\n\
                 rho_iota, rho_Omega used in Paper I; then exit.\n\
\n\
Shorthands:\n\
  --ours         --cross_lim, (1,4) and (1,5/2) on, (1,3/2) OFF, (1,7/2)=KST\n\
  --lr_full      --cross_lim, everything LR publishes, incl. (A7)-(A11)\n";

// MODIFIED (Arakida+): machine-precision dump of the full state, in the current
// unit system, for direct comparison against our Python runs.  Everything the
// papers plot (jbar_z, e0, omega_1, Omega_1, iota) is a function of these, so
// the derived quantities are computed on the Python side rather than here --
// in particular the running average over one outer period that turns these
// osculating elements into the MEAN elements the papers compare.
//
//  1:t  2-4:j1  5-7:e1  8-10:j2  11-13:e2  14:a1
//  15:jz=j1.n2  16:iota  17:omega_1(mutual node)  18:Omega_1  19:e1  20:e2
int RAW_OUTPUT = 0;
int AUDIT_ONLY = 0;
void print_state_raw(double t_next, kozai_struct *kozai, ostream& stream = cerr){
  stream << setprecision(17);
  double *y = kozai->get_y();
  stream << (t_next/YEAR);
  for(int i=0; i<13; i++) stream << "  " << y[i];
  stream << "  " << kozai->get_jz()
         << "  " << kozai->get_inc()
         << "  " << kozai->get_omega1_mutual()
         << "  " << kozai->get_Omega1_acc()
         << "  " << kozai->get_ecc1()
         << "  " << kozai->get_ecc2() << endl;
}

// Function for printing the state of a triple; placed here for easy
// modification
void print_state(double t_next, kozai_struct *kozai, ostream& stream = cerr){

  if (RAW_OUTPUT){ print_state_raw(t_next, kozai, stream); return; }

  streamsize ss = stream.precision();

  stream << setprecision(14) << (t_next/YEAR) << "  "
    << setprecision(8) 
    << kozai->get_a1()/AU << "  "
    << kozai->get_a2()/AU << "  " 
    << kozai->get_ecc1() << "  "
    << kozai->get_ecc2() << "  " 
    << setprecision(14) 
    << kozai->get_inc()/DEG << "  " 
    << kozai->get_inc1()/DEG << "  " 
    << kozai->get_inc2()/DEG << "  " 
    << setprecision(8) 
    << kozai->get_g1()/DEG << "  " 
    << kozai->get_g2()/DEG << "  "
    << kozai->get_Omega1()/DEG << "  "
    << kozai->get_Omega2()/DEG << endl;
}

//=====================================================================
// MODIFIED (Arakida+): --audit
//
// Integrating a triple mixes every grade together, and the small grades are
// invisible in max(e1).  The audit instead evaluates the right-hand side ONCE,
// at the initial condition, with one grade switched on at a time, and reports
//
//     Gdot = L0 d|j1|/dt ,   iotadot ,   Omegadot_1
//
// together with the two NORMALISATION-FREE ratios used in Paper I, Sec. on the
// grade (1,7/2):
//
//     rho_iota  = G0 iotadot   / Gdot ,      rho_Omega = G0 Omegadot_1 / Gdot .
//
// Because they are ratios they do not depend on the overall prefactor or on
// any sign convention, so they isolate the ANGULAR STRUCTURE -- which is where
// LR (A7),(A9),(A11) and KST Eq. (35) disagree.  Paper I quotes, at
// iota = 1.1, omega_1 = 0.7, omega_2 = 1.3:
//
//     KST Eq. (35) :  rho_iota = 1.46406 ,  rho_Omega = 2.59997
//     LR  (A7)...  :  rho_iota = 2.74027 ,  rho_Omega = 5.39947
//
// The two ratios differ by different factors (1.872 vs 2.077), so the
// disagreement cannot be absorbed into a single overall constant.
//=====================================================================
static inline double wrap_pi(double x){
  while(x >  PI) x -= TWOPI;
  while(x < -PI) x += TWOPI;
  return x;
}

// Elements the papers use, read off the state vector currently in `kozai'.
struct audit_el { double G0, iota, Omega1, omega1; };
static audit_el elements_of(kozai_struct *kozai){
  audit_el e;
  vec j1 = kozai->get_j1();
  e.G0     = kozai->get_L1_no_a()*sqrt(kozai->get_a1())*abs(j1);
  e.iota   = kozai->get_inc();
  e.Omega1 = kozai->get_Omega1_acc();
  e.omega1 = kozai->get_omega1_mutual();
  return e;
}

// One row of the audit.  The rates are obtained by a centred difference along
// the flow, which is both simpler and safer than hand-differentiating the
// element definitions -- and it automatically picks up the fact that LR's
// variables are NOT contact elements, i.e. that a1 (hence L0, hence G0) moves.
static void audit_one(kozai_struct *kozai, const char *label){
  double *y = kozai->get_y();
  double f[DIMENSION], y0[DIMENSION];
  for(int i=0;i<DIMENSION;i++) y0[i] = y[i];
  rhs(0., y0, f, (void *) kozai);

  // step scale: j1,e1,j2,e2 are O(1) by construction, a1 is not
  double sc = 0.;
  for(int i=0;i<12;i++) sc = max(sc, fabs(f[i]));
  sc = max(sc, fabs(f[12])/max(fabs(y0[12]), 1e-300));
  if(sc == 0.){
    cout << setw(30) << left << label << right << scientific << setprecision(8)
         << setw(16) << 0. << setw(16) << 0. << setw(16) << 0. << setw(16) << 0.
         << setw(12) << "-" << setw(12) << "-" << endl;
    return;
  }
  double h = 1e-7/sc;

  for(int i=0;i<DIMENSION;i++) y[i] = y0[i] + h*f[i];
  audit_el ep = elements_of(kozai);
  for(int i=0;i<DIMENSION;i++) y[i] = y0[i] - h*f[i];
  audit_el em = elements_of(kozai);
  for(int i=0;i<DIMENSION;i++) y[i] = y0[i];
  audit_el e0 = elements_of(kozai);

  double Gdot = (ep.G0     - em.G0)/(2.*h);
  double idot = (ep.iota   - em.iota)/(2.*h);
  double Odot = wrap_pi(ep.Omega1 - em.Omega1)/(2.*h);
  double wdot = wrap_pi(ep.omega1 - em.omega1)/(2.*h);

  cout << setw(30) << left << label << right << scientific << setprecision(8)
       << setw(16) << Gdot << setw(16) << idot
       << setw(16) << Odot << setw(16) << wdot;
  if (fabs(Gdot) > 1e-300*fabs(e0.G0))
    cout << fixed << setprecision(5) << setw(12) << e0.G0*idot/Gdot
                                     << setw(12) << e0.G0*Odot/Gdot;
  else cout << setw(12) << "-" << setw(12) << "-";
  cout << endl;
}

static void run_audit(kozai_struct *kozai){
  cout << "\n=== term audit: rhs evaluated once at t = 0 ===\n";
  cout << "  iota = " << kozai->get_inc() << " rad,  omega_1 = "
       << kozai->get_omega1_mutual() << " rad,  e1 = " << kozai->get_ecc1()
       << ",  e2 = " << kozai->get_ecc2() << "\n\n";
  cout << setw(30) << left << "grade / block" << right
       << setw(16) << "Gdot" << setw(16) << "iotadot"
       << setw(16) << "Omegadot_1" << setw(16) << "omegadot_1"
       << setw(12) << "rho_iota" << setw(12) << "rho_Omega" << endl;
  cout << string(118,'-') << endl;

  // remember what the user asked for
  bool q=kozai->get_quadrupole(), o=kozai->get_octupole(),
       p=kozai->get_pericenter(), P=kozai->get_outer_pericenter(),
       X=kozai->get_cross_lim(),  r=kozai->get_radiation(),
       b=kozai->get_brown();
  int  x72=kozai->get_x_rn3p5();

  #define ONLY(qq,oo,pp,PP,XX,r4,r25,r15,r35,bb) \
    kozai->set_quadrupole(qq); kozai->set_octupole(oo); \
    kozai->set_pericenter(pp); kozai->set_outer_pericenter(PP); \
    kozai->set_cross_lim(XX);  kozai->set_radiation(false); \
    kozai->set_x_rn4(r4); kozai->set_x_rn2p5(r25); kozai->set_x_rn1p5(r15); \
    kozai->set_x_rn3p5(r35); kozai->set_brown(bb);

  ONLY(true ,false,false,false,false,false,false,false,0,false)
  audit_one(kozai, "(0,3)   Newtonian quadrupole");
  ONLY(false,true ,false,false,false,false,false,false,0,false)
  audit_one(kozai, "(0,4)   Newtonian octupole");
  ONLY(false,false,false,false,false,false,false,false,0,true )
  audit_one(kozai, "(0,9/2) Brown  [ours]");
  ONLY(false,false,true ,false,false,false,false,false,0,false)
  audit_one(kozai, "(1,0)   inner 1PN pericentre");
  ONLY(false,false,false,false,true ,true ,false,false,0,false)
  audit_one(kozai, "(1,4)   Rn4  [LR, canonical]");
  ONLY(false,false,false,false,true ,false,true ,false,0,false)
  audit_one(kozai, "(1,5/2) Rn2p5 de Sitter, LR(4.4)");
  kozai->set_lr_rn2p5_exact(true);
  audit_one(kozai, "(1,5/2) same, LR mass prefactor");
  kozai->set_lr_rn2p5_exact(false);
  ONLY(false,false,false,false,true ,false,false,true ,0,false)
  audit_one(kozai, "(1,3/2) Rn1p5 LR(4.5) DISPUTED");
  kozai->set_lr_rn1p5_exact(true);
  audit_one(kozai, "(1,3/2) same, LR mass prefactor");
  kozai->set_lr_rn1p5_exact(false);
  ONLY(false,false,false,false,false,false,false,false,1,false)
  audit_one(kozai, "(1,7/2) Rn3p5  ours = -K_star");
  ONLY(false,false,false,false,false,false,false,false,2,false)
  audit_one(kozai, "(1,7/2) Rn3p5  LR reconstructed");
  ONLY(false,false,false,false,true ,false,false,false,3,false)
  audit_one(kozai, "(1,7/2) Rn3p5  LR (A7)(A9)(A11)");
  #undef ONLY

  cout << "\n  K^(1,3/2) should be identically zero (Paper I): the row above is\n"
          "  the spurious libration term, and --no_rn1p5 removes it.\n";
  cout << "  The two (1,7/2) rows share every prefactor; compare rho_iota and\n"
          "  rho_Omega, which are independent of the overall normalisation.\n\n";

  // restore
  kozai->set_quadrupole(q); kozai->set_octupole(o); kozai->set_pericenter(p);
  kozai->set_outer_pericenter(P); kozai->set_cross_lim(X);
  kozai->set_radiation(r); kozai->set_brown(b); kozai->set_x_rn3p5(x72);
  kozai->set_x_rn4(true); kozai->set_x_rn2p5(true); kozai->set_x_rn1p5(true);
}

int main(int argc, char **argv){

  double eps_abs = 1.e-12;
  double eps_rel = 1.e-12;
  double h = 1.e-10;
  double tmin=0, tmax=3.2e8*YEAR, delta_t=1e5*YEAR;
  bool collision=false;
  bool failure=false;
  bool IGNORE_GSL_ERRORS=false;
  bool PRINT_GSL=false;
  bool TENhz_ecc=false;
  bool TENhz_circ=false;
  double t=tmin;
  double t_insp;
  double max_e=0.;
  int status;

  //MODIFIED (Arakida+): --units / --clight must take effect BEFORE any option
  //is parsed, because set_parameters() multiplies the masses by MSUN, the
  //semi-major axes by AU and the times by YEAR as it reads them.  So pre-scan.
  int unit_mode = 0;
  {
    double clight = -1.;
    for(int i=1; i<argc; i++){
      if(strcmp(argv[i],"--units")==0 && i+1<argc)
        unit_mode = (strcmp(argv[i+1],"code")==0) ? 1 : 0;
      if(strcmp(argv[i],"--clight")==0 && i+1<argc)
        clight = atof(argv[i+1]);
    }
    set_units(unit_mode, clight);
    if(unit_mode == 1){
      //sensible code-unit defaults: our Paper II configuration
      tmax = 1e5; delta_t = 1e2;
    }
  }

  //First, initialize the Kozai structure
  kozai_struct *kozai = new kozai_struct;

  //MODIFIED: in code units the built-in default triple (a BH triple in AU/Msun)
  //is meaningless, so fall back on the configuration of our Paper II:
  //  m = (1.5, 0.5, 7.0), a0 = 1, a3 = 14, e3 = 0.6, jz = 0.4, e0 = e0*,
  //  omega_1(0) = 0.30 rad, Omega_1(0) = 0, varpi_3 = 0.
  if (unit_mode == 1){
    kozai->set_m1(1.5); kozai->set_m2(0.5); kozai->set_m3(7.0);
    kozai->set_a1(1.0); kozai->set_a2(14.0);
    kozai->set_ecc1(0.24451925); kozai->set_ecc2(0.6);
    //jz = l0 cos(iota) = 0.4  ->  iota = acos(0.4/sqrt(1-e0^2))
    kozai->set_inc(acos(0.4/sqrt(1.-0.24451925*0.24451925)));
    kozai->set_g1(0.30); kozai->set_g2(0.0);
    kozai->set_Omega1(0.0); kozai->set_Omega2(PI);
    kozai->set_outer_frame(true);
    kozai->set_restricted(true);
  }

  set_parameters(argc, argv, kozai, tmax, delta_t, IGNORE_GSL_ERRORS, PRINT_GSL);
  kozai->initialize();

  //Then, set up the GSL integrator (we're using an 8th-order Runge Kutta)
  const gsl_odeiv2_step_type *type_ptr = gsl_odeiv2_step_rk8pd;

  gsl_odeiv2_step *step_ptr = gsl_odeiv2_step_alloc(type_ptr, DIMENSION);
  gsl_odeiv2_control *control_ptr = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
  gsl_odeiv2_evolve *evolve_ptr = gsl_odeiv2_evolve_alloc(DIMENSION);
  gsl_odeiv2_system my_system;

  my_system.function = rhs;
  my_system.dimension = DIMENSION;
  my_system.params = kozai;

  if (AUDIT_ONLY){ print_header_and_initial_state(kozai); run_audit(kozai); delete kozai; return 0; }

  print_header_and_initial_state(kozai);
  print_state(0.,kozai);


  //Then actually integrate the system forward over each timestep
  long step = 0;
  long delta_output = 1;
  bool output_every_nsteps = false;

  for (double t_next = delta_t; t_next < tmax; t_next += delta_t){

    // if printing every x timesteps, the while loop does the whole
    // integration
    if (delta_t <= 0){
      delta_output = pow(10,-delta_t);
      delta_t = 1e100;
      t_next = tmax;
      output_every_nsteps = true;
    }

    //Then integrate the system forward for this timestep
    while (t < t_next){
      status = gsl_odeiv2_evolve_apply (evolve_ptr, control_ptr, step_ptr,
                                   &my_system, &t, t_next, &h, kozai->get_y());

      //Check for collisions
      if(kozai->collision()){
        collision = true;
        break;
      }

      //Record the maximum ecc of the inner binary
      if (kozai->get_ecc1() > max_e)
        max_e = kozai->get_ecc1();
      
      //check for hitting the LIGO band
      if(kozai->gwave_freq() > 10 && TENhz_ecc == false){
        TENhz_ecc = true;
        cout << "LIGO band (eccentric): ";
        print_state(t,kozai,cout);
      }
      
      if(2*kozai->get_forb() > 10){
        TENhz_circ = true;
        cout << "LIGO band (circular): ";
        print_state(t,kozai,cout);
        break;
      }

      //print status if PRINT_GSL flagged
      if ((status != GSL_SUCCESS && PRINT_GSL) ){
        cout << "GSL_ERROR: (time, status) = " << t / YEAR << "," << gsl_strerror(status) << endl;
        print_state(t,kozai);
      }
      //check for some type of failure/NaN
      if ((status != GSL_SUCCESS && !IGNORE_GSL_ERRORS) || isnan(kozai->get_ecc1()) || kozai->get_ecc1() >= 1.){
        failure = true;
        break;
      }

      //otherwise, print the step (if we're printing on the integrator
      //timescale)
      if (output_every_nsteps && step % delta_output == 0)
        print_state(t,kozai);
        /*
        cout << "t = " << t/YEAR << endl;
        for(int z = 0;z <= 12; z++) cout << kozai->get_y()[z] << " ";
        cout << endl;
        for(int z = 0;z <= 12; z++) cout << evolve_ptr->yerr[z] << " " ;
        cout << endl << endl;
        */
      step += 1;
    }

    //Finally, print the state of the system
    print_state(t,kozai);

    if (collision == true){
      cout << "Collision!" << endl;
      break;
    }

    if (TENhz_circ == true){
      cout << "LIGO band!" << endl;
      break;
    }

    if (failure == true){
      cout << "NaNs (if near merger, try --ignore_gsl)!" << endl;
      break;
    }
  }

  streamsize ss = cout.precision();
  cout << setprecision(10);
  cout << "Maximum Eccentricity: " << max_e << endl;
  cout << setprecision(ss);

  /* all done; free up the gsl_odeiv stuff */
  gsl_odeiv2_evolve_free (evolve_ptr);
  gsl_odeiv2_control_free (control_ptr);
  gsl_odeiv2_step_free (step_ptr);
  delete kozai;
}

void print_header_and_initial_state(kozai_struct* kozai){
  cout << "Triple setup:\n"  
    << "  m1 = " << kozai->get_m1()/MSUN << endl
    << "  m2 = " << kozai->get_m2()/MSUN << endl
    << "  m3 = " << kozai->get_m3()/MSUN << endl
    << "  a1 = " << kozai->get_a1()/AU << endl
    << "  a2 = " << kozai->get_a2()/AU << endl
    << "  e1 = " << kozai->get_ecc1() << endl
    << "  e2 = " << kozai->get_ecc2() << endl
    << "  g1 = " << kozai->get_g1()/DEG << endl
    << "  g2 = " << kozai->get_g2()/DEG << endl
    << "  Omega1 = " << kozai->get_Omega1()/DEG << endl
    << "  Omega2 = " << kozai->get_Omega2()/DEG << endl
    << "  inc = " << kozai->get_inc()/DEG << endl
    << "  inc 1 = " << kozai->get_inc1()/DEG << endl
    << "  inc 2 = " << kozai->get_inc2()/DEG << endl
    << "  radius 1 = " << kozai->get_r1()/RSUN << endl
    << "  radius 2 = " << kozai->get_r2()/RSUN << endl;
  
  //Print the header with units and column names, and the initial state of the
  //system
  // cerr << "#1:t[yr]  #2:a1[AU]  #3:a2[AU] #4:e1  #5:e2  #6:inc[deg] #7:i1  #8:i2  #9:g1  #10:g2" << endl;
  if (RAW_OUTPUT)
    cerr << "#1:t  #2-4:j1  #5-7:e1  #8-10:j2  #11-13:e2  #14:a1  "
            "#15:jz  #16:iota[rad]  #17:omega1[rad]  #18:Omega1[rad]  #19:e1  #20:e2" << endl;
  else
    cerr << "#1:t[yr]  #2:a1[AU]  #3:a2[AU] #4:e1  #5:e2  #6:inc[deg] #7:i1  #8:i2  #9:g1  #10:g2  #11:h1  #12:h2" << endl;
}


int rhs(double t, const double y[], double f[], void *kozai_ptr){

  kozai_struct *kozai = (kozai_struct *) kozai_ptr;
  //First, convert the y array into vectors
  //
  //Note: you can't use the vectors in the kozai structure; for adaptive
  //timstepping, GSL sometimes uses a temporary y array, which won't be
  //updated in the kozai class
  vec j1 = vec(y[0],y[1],y[2]);
  vec e1 = vec(y[3],y[4],y[5]);
  vec j2 = vec(y[6],y[7],y[8]);
  vec e2 = vec(y[9],y[10],y[11]);

  //extract semi-major axes, masses and angular momenta
  double a1 = y[12];
  double a2 = kozai->get_a2();
  double m1 = kozai->get_m1();
  double m2 = kozai->get_m2();
  double m3 = kozai->get_m3();
  double L1 = kozai->get_L1_no_a()*sqrt(a1);
  double L2 = kozai->get_L2_no_a()*sqrt(a2);

  //pre-compute a bunch of the vectorial quantities that are used repeatedly
  double j1n = abs(j1);
  double j2n = abs(j2);
  vec n1 = j1/j1n;
  vec n2 = j2/j2n;

  double e1n = abs(e1);
  double e2n = abs(e2);
  vec u1 = e1/e1n;
  vec u2 = e2/e2n;

  double j1n2 = j1*n2;
  double j1u2 = j1*u2;
  double e1n2 = e1*n2;
  double e1u2 = e1*u2;

  vec j1xn2 = j1^n2;
  vec j1xu2 = j1^u2;
  vec e1xn2 = e1^n2;
  vec e1xu2 = e1^u2;

  double L1L2 = L1 / L2;

  //Compute the (local) quadrupole timescale
  double tsec = (sqrt((m1+m2)/(G*pow(a1,3.)))*pow(a2,3.)*pow(1.-sqr(e2n),1.5) /m3);

  vec dj1dt=0., de1dt=0., dj2dt=0., de2dt=0.;
  double dadt  = 0;

  //Add the quadrupole-order secular evolution equations
  if(kozai->get_quadrupole()){
    double quad_coef = 0.75 / tsec;
    dj1dt += quad_coef*(j1n2*j1xn2 - 5.*e1n2*e1xn2);
    de1dt += quad_coef*((j1n2*e1xn2) + (2.*(j1^e1)) - (5.*(e1n2*j1xn2)));
    dj2dt += quad_coef*L1L2*((-j1n2*j1xn2) + (5*e1n2*e1xn2));
    de2dt += quad_coef*(L1L2/j2n)*((j1n2*(e2^j1)) - (5.*e1n2*(e2^e1)) 
         - ((0.5 - 3.*sqr(e1n) + 12.5*sqr(e1n2) - 2.5*sqr(j1n2))*(n2^e2)));
  }

  //Add the octupole-order secular evolution equations
  if(kozai->get_octupole() == true){
    double octo_coef = -1.171875*(fabs(m1-m2)/(m1+m2)*(a1/a2)*(e2n/sqr(j2n)))/tsec;
    dj1dt += octo_coef*((((2.*(e1u2*j1n2+e1n2*j1u2)*j1)
        + 2.*(j1u2*j1n2-7.*e1u2*e1n2)*e1)^n2) + ((((2.*e1n2*j1n2)*j1)
        + (1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2) + sqr(j1n2))*e1)^u2));
    de1dt += octo_coef*(((((2.*e1n2*j1n2)*e1) + (1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2)   
        + sqr(j1n2))*j1)^u2) + (((2.*(e1u2*j1n2+e1n2*j1u2)*e1) 
        + (2.*(j1n2*j1u2-7*e1n2*e1u2)*j1))^n2) + (3.2*e1u2*(j1^e1)));
    dj2dt += octo_coef*L1L2*((2*(((e1n2*j1u2)*n2) + ((e1u2*j1n2)*n2)
        + (e1n2*j1n2)*u2)^j1) + ((((2*j1u2*j1n2)*n2) - ((14.*e1u2*e1n2)*n2) 
        + ((1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2) + sqr(j1n2))*u2))^e1));
    de2dt += octo_coef*(L1L2/j2n)*(2*(((e1n2*((j1*e2)*u2))
        + (j1n2*((e1*e2)*u2)) + (((sqr(j2n)/e2n)*e1n2*j1n2)*n2))^j1) + ((((2.*e2n*j1u2*j1n2)*u2)
        - ((14.*e1n2*(e1*e2))*u2) + (sqr(j2n)/e2n)*((1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2)
        + sqr(j1n2))*n2))^e1) - (((2.*(0.2-1.6*sqr(e1n))*e1u2*e2) + ((14.*e1n2*j1u2*j1n2)*e2)
        + (7.*e1u2*(1.6*sqr(e1n) - 0.2 - 7.*sqr(e1n2) + sqr(j1n2))*e2))^n2));
  }

  //Add the pericenter precession of the inner binary
  if (kozai->get_pericenter() == true)
    de1dt += (3./(c*c*a1*sqr(j1n)))*(pow(G*(m1+m2)/a1, 1.5) * (n1^e1));

  //Add the pericenter precession of the outer binary
  if (kozai->get_outer_pericenter() == true)
    de2dt += (3./(c*c*a2*sqr(j2n)))*(pow(G*(m1+m2+m3)/a2, 1.5) * (n2^e2));

  //Add 1pn-quadrupole cross terms
  if (kozai->get_cross_lim() == true){   
    // precompute various quantities
    vec v1 = n1^u1;
    vec v2 = n2^u2;
    double p2    = a2 * sqr(j2n);
    double m = m1+m2;
    double eta = m1*m2/sqr(m);
    double mtot = m + m3;
        
    // Define trig funcs wrt to inc, i1
    // Avoid evaluating angles; use vectors; Naoz et al (2013a) Appendix A
    double cinc    = n1*n2;
    double sinc    = sqrt(1.-sqr(cinc));
    double s2inc   = 2. * sinc * cinc;
    double sincsq  = sqr(sinc);
    double cincsq  = sqr(cinc);
    double c2inc   = 1 - 2*sqr(sinc);
    double cinc1   = n1[2];
    // double cinc2   = n2[2];
    double sinc1   = sqrt(1.-sqr(cinc1));
    double cscinc1 = 1 / sinc1;
    // double sinc2   = sqrt(1.-sqr(cinc2));
    double cotinc1     = cinc1 / sinc1;
    double cotinc1sinc = cinc1 * sinc / sinc1;
    double cscinc1sinc = sinc / sinc1;

    // Define trig functions wrt to g1, g2; Goes bad as sinc -> 0
    double cg1 = (n2*v1) / sinc;
    double sg1 = (n2*u1) / sinc;
    double cg2 = (n1*v2) / sinc;
    double sg2 = (n1*u2) / sinc;
    double c2g1 = 1 - 2*sqr(sg1);
    double s2g1 = 2 * sg1 * cg1;
    double c2g2 = 1 - 2*sqr(sg2);
    double s2g2 = 2 * sg2 * cg2;

    double cinc2 = n2[2];
    double sinc2 = sqrt(1.-sqr(cinc2));
    // 1PN-quadrupole cross terms; Lim Rodriguez (2019)
    // Inner binary
    double de1dtcross=0, dp1dtcross=0, dg1dtcross=0, dh1dtcross=0, di1dtcross=0;

    //=================================================================
    // MODIFIED (Arakida+): each grade is now behind its own switch so
    // that the disputed pieces can be isolated.  See README.md.
    //   x_rn4   -> blocks (1,4)      : LR Eqs. (4.6), (A12)-... ; canonical  [OK]
    //   x_rn2p5 -> block  (1,5/2)    : LR Eq.  (4.4) de Sitter  ; canonical  [OK]
    //   x_rn1p5 -> block  (1,3/2)    : LR Eq.  (4.5) libration  ; SPURIOUS
    //=================================================================
    if (kozai->get_x_rn4()) {
    /* 
    // Original: Just IDPNB and ID6, (nr, nm) = (1.5, -0.5) 
    double factor = pow(a1,1.5)*pow(a2,-4)*pow(c,-2)*pow(G,1.5)*pow(j2n,-5)*pow(m,-0.5)*pow(mtot,2);
    de1dtcross += factor * (-15*e1n*j1n*(8*c2g1*cinc*s2g2 - 2*s2g1*(c2g2*(3 + c2inc) + 6*sincsq))*pow(e2n,2))/32.;
    dp1dtcross += factor * (-15*a1*j1n*(-4*c2g1*cinc*s2g2 + s2g1*(c2g2*(3 + c2inc) + 6*sincsq))*pow(e1n,2)*pow(e2n,2))/8.;
    di1dtcross += factor * (-3*sinc*pow(e2n,2)*pow(j1n,-1)*(-5*(-3 + c2g2)*cinc*s2g1*pow(e1n,2) + s2g2*(5 + 5*c2g1*pow(e1n,2) - 3*pow(j1n,2))))/8.;
    dg1dtcross += factor * (-3*pow(j1n,-1)*(10*cotinc1sinc*s2g1*s2g2*pow(e1n,2)*(-16 + 3*pow(j2n,2)) + 20*cinc*s2g1*s2g2*pow(j1n,2)*(-16 + 3*pow(j2n,2)) - cotinc1*s2inc*(5 - 5*c2g1*pow(e1n,2) - 3*pow(j1n,2))*(4 - 10*pow(j2n,2) + c2g2*(-16 + 3*pow(j2n,2))) + (-3 + 5*c2g1)*c2inc*pow(j1n,2)*(4 - 10*pow(j2n,2) + c2g2*(-16 + 3*pow(j2n,2))) + (1 + 5*c2g1)*pow(j1n,2)*(-4 + 10*pow(j2n,2) + c2g2*(-48 + 9*pow(j2n,2)))))/64.;
    dh1dtcross += factor * (3*cscinc1*sinc*pow(e2n,2)*pow(j1n,-1)*(-5*s2g1*s2g2*pow(e1n,2) + (-3 + c2g2)*cinc*(5 - 5*c2g1*pow(e1n,2) - 3*pow(j1n,2)))*pow(j2n,-5))/8.;
   */ 

    // Modified: IDPNB, ID6 and DN6 contributions to (nr, nm) = (1.5, -0.5) 
    // IDPNB + ID6 contributions
    double factor = pow(a1,1.5)*pow(a2,-4)*pow(c,-2)*pow(G,1.5)*pow(j2n,-5)*pow(m,-0.5)*pow(mtot,2);
    de1dtcross += factor * (15*e1n*j1n*(-4*c2g1*cinc*s2g2 + s2g1*((-3 + c2g2)*c2inc + 6*pow(cg2,2)))*pow(e2n,2))/16.;
    dp1dtcross += factor * (-15*a1*j1n*(-4*c2g1*cinc*s2g2 + s2g1*((-3 + c2g2)*c2inc + 6*pow(cg2,2)))*pow(e1n,2)*pow(e2n,2))/8.;
    di1dtcross += factor * (3*sinc*pow(e2n,2)*pow(j1n,-1)*(5*(-3 + c2g2)*cinc*s2g1*pow(e1n,2) + s2g2*(-5 - 5*c2g1*pow(e1n,2) + 3*pow(j1n,2))))/8.;
    dg1dtcross += factor * (3*pow(e2n,2)*pow(j1n,-1)*(10*cotinc1*s2g1*s2g2*sinc*pow(e1n,2) + 3*pow(j1n,2) + 9*c2inc*pow(j1n,2) + 20*cinc*s2g1*s2g2*pow(j1n,2) + 6*c2g2*sincsq*pow(j1n,2) + (-3 + c2g2)*cotinc1*s2inc*(-5 + 3*pow(j1n,2)) + 5*c2g1*(6*pow(cg2,2)*pow(j1n,2) + (-3 + c2g2)*(cotinc1*s2inc*pow(e1n,2) + c2inc*pow(j1n,2)))))/16.;
    dh1dtcross += factor * (-3*cscinc1sinc*pow(j1n,-1)*(5*s2g1*s2g2*pow(e1n,2)*pow(e2n,2) + (-3 + c2g2)*cinc*pow(e2n,2)*(-5 + 5*c2g1*pow(e1n,2) + 3*pow(j1n,2))))/8.;

    // DN6 contributions
    de1dtcross += factor * (15*e1n*j1n*(4*c2g1*cinc*s2g2*pow(e2n,2) - s2g1*(12*sincsq + c2g2*(3 + c2inc)*pow(e2n,2))))/32.;
    dp1dtcross += factor * (15*a1*j1n*pow(e1n,2)*(12*s2g1*sincsq + (c2g2*(3 + c2inc)*s2g1 - 4*c2g1*cinc*s2g2)*pow(e2n,2)))/16.;
    di1dtcross += factor * (-3*sinc*pow(j1n,-1)*(5*cinc*s2g1*pow(e1n,2)*(-6 + c2g2*pow(e2n,2)) + s2g2*pow(e2n,2)*(-5 - 5*c2g1*pow(e1n,2) + 3*pow(j1n,2))))/16.;
    dg1dtcross += factor * (-3*pow(j1n,-1)*(10*cotinc1*s2g1*s2g2*sinc*pow(e1n,2)*pow(e2n,2) + 6*sincsq*(10*c2g1 + c2g2*pow(e2n,2))*pow(j1n,2) + (6 + 18*c2inc + 5*(c2g1*c2g2*(3 + c2inc) + 4*cinc*s2g1*s2g2)*pow(e2n,2))*pow(j1n,2) + cotinc1*s2inc*(-6 + c2g2*pow(e2n,2))*(-5 + 5*c2g1*pow(e1n,2) + 3*pow(j1n,2))))/32.;
    dh1dtcross += factor * (3*cscinc1sinc*pow(j1n,-1)*(5*s2g1*s2g2*pow(e1n,2)*pow(e2n,2) + cinc*(-6 + c2g2*pow(e2n,2))*(-5 + 5*c2g1*pow(e1n,2) + 3*pow(j1n,2))))/16.;

    // Cross terms due to Period Corrections
    double quad_coef = (-3 * G * mtot * (-3 + sqr(j2n)) * (-5 + 2 * sqr(j2n)) / (2 * sqr(c) * p2 * sqr(j2n)) ) * 0.75 / tsec;
    dj1dt += quad_coef*(j1n2*j1xn2 - 5.*e1n2*e1xn2);
    de1dt += quad_coef*((j1n2*e1xn2) + (2.*(j1^e1)) - (5.*(e1n2*j1xn2)));
    }   // end x_rn4

    // (nr, nm) = (0, 0)   == Rn2p5 DN6 == LR Eq. (4.4), de Sitter precession
    //
    // NOTE (Arakida+): the prefactor hard-wired here is 3 M^{3/2}/2, which is the
    // TEST-PARTICLE limit of LR's published A = G^{3/2}(4m + 3 m3) m3 /(2 M^{1/2}
    // c^2 a2^{5/2} l2^2).  The two agree only as m -> 0.  --lr_rn2p5 restores LR's
    // published mass dependence.
    if (kozai->get_x_rn2p5()) {
    double dS_pref = 1.5*pow(mtot,1.5);
    if (kozai->get_lr_rn2p5_exact()) dS_pref = (4.*m + 3.*m3)*m3/(2.*sqrt(mtot));
    dh1dtcross += (dS_pref*cscinc1sinc*pow(a2,-2.5)*pow(c,-2)*pow(G,1.5)*pow(j2n,-2));
    dg1dtcross += (dS_pref*(cinc - cotinc1sinc)*pow(a2,-2.5)*pow(c,-2)*pow(G,1.5)*pow(j2n,-2));
    }   // end x_rn2p5

    // (nr, nm) = (-1, 1)  == Rn1p5 IDNQ == LR Eq. (4.5)/(A59), the "libration" term
    //
    // The angle factor below is exactly 4*AA2 with
    //   AA2 = A1 A3 - A2 A4 = cos(i) cos2w1 cos2w2 + (1+cos^2 i)/2 sin2w1 sin2w2,
    // and (1 + l2 - 2 l2^2)/(1 + l2) = (1 - l2)(1 + 2 l2)/(1 + l2) is LR's outer
    // eccentricity factor -- which our Paper I identifies as <cos 2F> taken with a
    // UNIFORM-F measure instead of the time measure; the correct value of this
    // whole grade is K^(1,3/2) == 0.  Off whenever --no_rn1p5 is given.
    //
    // NOTE: the mass prefactor here is m*sqrt(M), whereas LR Eq. (4.5) has
    // m*m3/sqrt(M).  --lr_rn1p5 applies the ratio m3/M to recover LR's published
    // expression.
    if (kozai->get_x_rn1p5()) {
    double lib = (-15*m*(4*c2g1*c2g2*cinc + (3 + c2inc)*s2g1*s2g2)*pow(a1,-1)*pow(a2,-1.5)*pow(c,-2)*pow(e1n,2)*pow(G,1.5)*pow(j1n,-3)*pow(j2n,-3)*(1 + j2n - 2*pow(j2n,2))*pow(1 + j2n,-1)*pow(mtot,0.5))/16.;
    if (kozai->get_lr_rn1p5_exact()) lib *= m3/mtot;
    dg1dtcross += lib;
    }   // end x_rn1p5

    // grade (1,7/2) == Rn3p5 DN6 == LR Eqs. (A7)-(A11), TRANSCRIBED VERBATIM.
    //
    // --x_rn3p5 lrraw.  Unlike the two Hamiltonian variants further down, this
    // reproduces LR's published rates exactly as they stand -- which is the
    // point: they are NOT the flow of any Hamiltonian (Paper I shows all three
    // mixed-partial identities fail), so no reconstructed K can match all of
    // them at once.  Only the three cells LR gives are implemented; LR's
    // appendix has no e1 or omega_1 cell at this grade, so those stay zero.
    //
    //   dp1/dt      = (3 S p1^2/4)  [ cos i (5cos^2 i - 1) sin w1 cos w2
    //                                + (3 - 7 cos^2 i) cos w1 sin w2 ]
    //   diota1/dt   = (3 S p1/16) sin i [ (16 - 10 cos^2 i) sin w1 cos w2
    //                                     + 14 cos i cos w1 sin w2 ]
    //   dOmega1/dt  = (3 S p1/8) csc i1 sin i [ 7 cos w1 cos w2
    //                                           + 17 cos i sin w1 sin w2 ]
    //   S = e1 e2 G^{3/2} M^{3/2} l2^3 sqrt(1-4 eta) / (c^2 p2^{7/2} l1^2)
    if (kozai->get_x_rn3p5() == 3){
    double p1  = a1*sqr(j1n);
    double dm  = fabs(m1-m2)/m;                  // sqrt(1 - 4 eta)
    double S3  = 3.*e1n*e2n*pow(G,1.5)*pow(mtot,1.5)*pow(j2n,3.)*dm
                 /(sqr(c)*pow(p2,3.5)*sqr(j1n));
    dp1dtcross += 0.25*S3*sqr(p1)*( cinc*(5.*cincsq - 1.)*sg1*cg2
                                  + (3. - 7.*cincsq)*cg1*sg2 );
    di1dtcross += (S3*p1/16.)*sinc*( (16. - 10.*cincsq)*sg1*cg2
                                   + 14.*cinc*cg1*sg2 );
    dh1dtcross += (S3*p1/8.)*cscinc1sinc*( 7.*cg1*cg2 + 17.*cinc*sg1*sg2 );
    }   // end x_rn3p5 == lrraw

    de1dt += de1dtcross*u1 + (e1n*(dg1dtcross + dh1dtcross*cinc1))*v1 + (e1n*(di1dtcross*sg1-dh1dtcross*cg1*sinc1))*n1;
    dj1dt += (j1n*(-di1dtcross*sg1 + dh1dtcross*cg1*sinc1))*u1 + (-j1n*(di1dtcross*cg1 + dh1dtcross*sg1*sinc1))*v1 + ((-e1n/j1n)*de1dtcross)*n1;
    dadt  += (dp1dtcross + 2. * a1 * e1n * de1dtcross) / sqr(j1n);
    }

  //=====================================================================
  // MODIFIED (Arakida+): terms that are NOT in Lim & Rodriguez's code.
  //
  // Both are added as genuine Hamiltonian flows through the Milankovitch
  // (vectorial Lagrange) equations
  //
  //    dj1/dt = -(1/L1) [ j1 x dPhi/dj1 + e1 x dPhi/de1 ]
  //    de1/dt = -(1/L1) [ e1 x dPhi/dj1 + j1 x dPhi/de1 ]
  //
  // with Phi the doubly-averaged Hamiltonian in the same sign convention as
  // the quadrupole block above.  (Calibration: writing the code's own
  // quadrupole block in this form gives
  //    Phi_quad = kappa [ (j1.n2)^2 - 5 (e1.n2)^2 + 2 e1^2 ],
  //    kappa = -3 L1 / (8 tsec) = 6 A_pref,  A_pref = -mu0 m3 a0^2/(16 a3^3 eta3^3),
  // which is exactly the normalisation C = 6 A_pref of our Paper II.  The two
  // blocks below therefore share the papers' normalisation with no free factor.)
  //
  // As with every other cross term in this code, the back-reaction on the outer
  // orbit is not included, so these blocks touch j1/e1 only.
  //=====================================================================

  // ---- grade (1,7/2) : Rn3p5, absent from this port ------------------
  //
  // K^(1,7/2) = s * N * Xi / l2^5 ,  with e1 the eccentricity vector of x1 - x2,
  //   N   = 3 G^2 (m1-m2) m1 m2 m3 sqrt(M3) a1^{3/2} / (c^2 M0^{3/2} a2^{7/2})
  //
  //   --x_rn3p5 kst  (s = -1):  K^(1,7/2) = -N Xi / l2^5  (KST Eq. (35), second
  //        term; Paper II Eqs. (45)-(47)).  The e1 of this code is that of
  //        x_light - x_heavy (octupole block), so the code applies
  //        +N' |m1-m2| Xi / l2^5 with N' = N/(m1-m2); see "Orientation" below.
  //        Xi     = (j1 x e1).(j2 x e2) - 2 (j1.j2)(e1.e2)
  //        <=> B    = cos(i) cos w1 cos w2 + (2cos^2 i - 1) sin w1 sin w2
  //   --x_rn3p5 lr   (s = +1/8, LR Eqs. (A7)-(A11) reconstructed):
  //        Xi_LR  = -2 (j1.j2)(e1.e2) + (j1 x e1).(j2 x e2) [5 cos^2 i - 3]
  //        <=> B_LR = (5cos^2 i - 1) cos i cos w1 cos w2 + (7 cos^2 i - 3) sin w1 sin w2
  //
  // The two share all eight prefactor structures (G^{3/2}M^{3/2}, e1 e2, p1^2,
  // l2^3, p2^{-7/2}, l1^{-2}, c^{-2} and the mass-difference factor); only the
  // angular structure and the numerical coefficient (-1 vs +1/8) differ, which
  // is precisely the disagreement documented in Paper I.
  if (kozai->get_x_rn3p5() == 1 || kozai->get_x_rn3p5() == 2){
    double mm   = m1 + m2;
    double mmtot= mm + m3;
    // Orientation: the octupole block above uses |m1-m2|, i.e. it treats e1 as the
    // eccentricity vector of x_light - x_heavy.  Physically K^(1,7/2) = -N Xi with
    // N ~ (m1-m2) and e1 of x1 - x2; written with the code's e1 this is +N' |m1-m2| Xi.
    // X72LEGACY=1 restores the sign of the submitted versions (correct only for m1 < m2).
    static const bool legacy72 = (getenv("X72LEGACY") && atoi(getenv("X72LEGACY")) == 1);
    double dmf  = legacy72 ? (m1-m2) : -fabs(m1-m2);
    double Nbase = 3.*sqr(G)*m1*m2*m3*sqrt(mmtot)*pow(a1,1.5)
                  /(sqr(c)*pow(mm,1.5)*pow(a2,3.5));
    double Nfac  = Nbase*(m1-m2);    // LR reconstruction: unchanged
    double inv5 = 1./pow(j2n,5.);
    vec P = j2^e2;                 // j2 x e2
    vec Q = j1^e1;                 // j1 x e1
    double R    = Q*P;
    double jj   = j1*j2;
    double ee12 = e1*e2;
    vec gj, ge;                    // dPhi/dj1 , dPhi/de1
    if (kozai->get_x_rn3p5() == 1){          // ours / KST
      double s = -1.0*Nbase*dmf*inv5;
      gj = s*((e1^P) - (2.*ee12)*j2);
      ge = s*((P^j1) - (2.*jj)*e2);
    } else {                                  // LR (A7)-(A11), reconstructed
      double s  = 0.125*Nfac*inv5;
      double cI = jj/(j1n*j2n);
      double fI = 5.*sqr(cI) - 3.;
      vec dcdj1 = (1./(j1n*j2n))*j2 - (cI/sqr(j1n))*j1;
      gj = s*(((-2.*ee12)*j2) + fI*(e1^P) + (10.*cI*R)*dcdj1);
      ge = s*(((-2.*jj)*e2)   + fI*(P^j1));
    }
    dj1dt += (-1./L1)*((j1^gj) + (e1^ge));
    de1dt += (-1./L1)*((e1^gj) + (j1^ge));
  }

  // ---- grade (0,9/2) : the Brown Hamiltonian (Newtonian, single-averaging) --
  //
  //   Phi_Brown = C eps21 F21,
  //   F21 = (3/8) jz (1 - jz^2 + 24 e1^2 - 15 ez^2),  jz = j1.n2,  ez = e1.n2,
  //   C   = 6 A_pref = -3 L1/(8 tsec),
  //   eps21 = (n3/n0)(m3/M3)(1 + 2 e3^2/3)/eta3^3 .
  //
  // This is the omega_2-INDEPENDENT part of our K^(0,9/2) (Luo, Katz & Dong
  // 2016; Tremaine 2023; Lei & Grishin 2025).  Being axisymmetric it conserves
  // jz exactly; its job in the comparison is to move the DeltaOmega_1 = -pi
  // resonant surface, which our Paper II locates at e0 = 0.24451925 without it
  // and e0 = 0.28512069 with it (a3 = 14, jz = 0.4).  The 2 omega_2 harmonic of
  // K^(0,9/2) -- the part that actually drives the libration -- has no closed
  // form and is not included here; see README.md.
  if (kozai->get_brown()){
    double mm   = m1 + m2;
    double mmtot= mm + m3;
    double Cq   = -3.*L1/(8.*tsec);                     // = 6 A_pref
    double n0   = sqrt(G*mm/pow(a1,3.));
    double n3   = sqrt(G*mmtot/pow(a2,3.));
    double eps21= (n3/n0)*(m3/mmtot)*(1. + 2.*(e2*e2)/3.)/pow(j2n,3.);
    double z    = j1*n2;
    double w    = e1*n2;
    double S    = 1. - sqr(z) + 24.*(e1*e1) - 15.*sqr(w);
    double pre  = 0.375*Cq*eps21;
    vec gj = (pre*(S - 2.*sqr(z)))*n2;
    vec ge = (pre*z)*(48.*e1 - (30.*w)*n2);
    dj1dt += (-1./L1)*((j1^gj) + (e1^ge));
    de1dt += (-1./L1)*((e1^gj) + (j1^ge));
  }

  //Add gravitational-wave emission for the inner binary
  if (kozai->get_radiation() == true){
    double dedt  = -(c304o15*e1n*(G3 * m1*m2*(m1+m2)) /
         (c5 * pow(a1,4.) * pow(j1n,5.))) * (1. + c121o304*sqr(e1n));
    dadt += -c64o5 * ((G3 * m1*m2*(m1+m2)) / (c5
        * pow(a1,3.) * pow(j1n,7.))) * (1. + c73o24*sqr(e1n) + c37o96*pow(e1n,4.)); 

    de1dt += dedt * u1;
    dj1dt += (-e1n / j1n) * dedt * n1;
  }

  //MODIFIED: --restricted freezes the outer orbit, which is the restricted
  //problem our Python single-/double-averaged runs integrate (the outer orbit is
  //an external, fixed perturber; L0 and the outer elements are constants).
  if (kozai->get_restricted()){ dj2dt = 0.; de2dt = 0.; }

  //Finally, copy the derivatives into the output array
  for(int i=0; i<3; i++){
    f[i] = dj1dt[i];
    f[i+3] = de1dt[i];
    f[i+6] = dj2dt[i];
    f[i+9] = de2dt[i];
  }
  f[12] = dadt;

  return GSL_SUCCESS;
}

// Integrate the peters equation to find the inspiral time for a binary
// Not actually used at the moment, but probably useful to have around...
double peters_t(kozai_struct *kozai){
  double beta = (64./5.) * G3c5 * kozai->get_m1()*kozai->get_m2()*(kozai->get_m1()+kozai->get_m2());
  double e = kozai->get_ecc1();
  double c0 = kozai->get_a1() * (1-sqr(e)) * pow(e,-c12o19) * pow(1+(c121o304)*sqr(e),-c870o2299);

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  gsl_function F;
  F.function = &peters_integral;
  double integral, error;

  gsl_integration_qags(&F,0,e,1e-7,1e-8,1000, w,&integral, &error);
  gsl_integration_workspace_free(w);

  return integral * c12o19 * sqr(c0)*sqr(c0) / beta;
}

double peters_integral(double e, void *params){
  return pow(e,c29o19) * pow(1+(c121o304)*sqr(e),c1181o2299) / pow(1-sqr(e),1.5);
}

void set_parameters(int argc, char **argv, kozai_struct *kozai, double &t_end, double &t_step, bool &IGNORE_GSL_ERRORS, bool &PRINT_GSL){
  int c;
  int opterr = 0;
  int option_index = 0;
  double radius;

  const struct option longopts[] =
  {
    {"m1",    1,  0, 'm'},
    {"m2",    1,  0, 'M'},
    {"m3",    1,  0, 'n'},
    {"a1",    1,  0, 'a'},
    {"a2",    1,  0, 'A'},
    {"e1",    1,  0, 'e'},
    {"e2",    1,  0, 'E'},
    {"g1",    1,  0, 'g'},
    {"g2",    1,  0, 'G'},
    {"omega1",1,  0, 'l'},
    {"omega2",1,  0, 'L'},
    {"inc"   ,1,  0, 'i'},
    {"rad1"   ,1,  0, 'b'},
    {"rad2"   ,1,  0, 'B'},
    {"quad"  ,0,  0, 'q'},
    {"oct"   ,0,  0, 'o'},
    {"peri"  ,0,  0, 'p'},
    {"outerperi",0,0,'P'},
    {"cross_lim",0,0,'X'},
    {"rad"   ,0,  0, 'r'},
    {"print_gsl",0,0,'f'},
    {"time"   ,1,  0, 'd'},
    {"dt"   ,1,  0, 'D'},
    {"ignore_gsl"   ,0,  0, 'I'},
    {"help"   ,0,  0, 'h'},
    //----- MODIFIED (Arakida+) -------------------------------------------
    {"units"      ,1, 0, 300},   // cgs (default) | code  (G = 1)
    {"clight"     ,1, 0, 301},   // speed of light in the current unit system
    {"no_rn4"     ,0, 0, 310},   // drop grade (1,4)
    {"no_rn2p5"   ,0, 0, 311},   // drop grade (1,5/2)  (de Sitter, LR 4.4)
    {"no_rn1p5"   ,0, 0, 312},   // drop grade (1,3/2)  (LR 4.5)  <-- our correction
    {"lr_rn1p5"   ,0, 0, 313},   // use LR's PUBLISHED mass prefactor for (1,3/2)
    {"lr_rn2p5"   ,0, 0, 314},   // use LR's PUBLISHED mass prefactor for (1,5/2)
    {"x_rn3p5"    ,1, 0, 315},   // kst | lr : add grade (1,7/2)
    {"brown"      ,0, 0, 316},   // add the (0,9/2) Brown Hamiltonian
    {"restricted" ,0, 0, 317},   // freeze the outer orbit
    {"outer_frame",0, 0, 318},   // reference plane = outer orbit plane
    {"raw"        ,0, 0, 319},   // machine-precision state dump
    {"ours"       ,0, 0, 320},   // shorthand: LR corrected as in our papers
    {"lr_full"    ,0, 0, 321},   // shorthand: LR as published, incl. (1,7/2)
    {"audit"      ,0, 0, 322},   // evaluate each grade once at t=0 and stop
    {0,0,0,0},
  };

  while ((c = getopt_long (argc, argv, 
          "m:M:n:a:A:e:E:g:G:l:L:i:c:C:t:T:u:U:qopPXsSrfId:D:b:B:h",
          longopts,&option_index)) != -1)
  switch (c)
    {
    case 'm':
      kozai->set_m1(atof(optarg)*MSUN); break;
    case 'M':
      kozai->set_m2(atof(optarg)*MSUN); break;
    case 'n':
      kozai->set_m3(atof(optarg)*MSUN); break;
    case 'a':
      kozai->set_a1(atof(optarg)*AU); break;
    case 'A':
      kozai->set_a2(atof(optarg)*AU); break;
    case 'e':
      kozai->set_ecc1(atof(optarg)); break;
    case 'E':
      kozai->set_ecc2(atof(optarg)); break;
    case 'g':
      kozai->set_g1(atof(optarg)*DEG); break;
    case 'G':
      kozai->set_g2(atof(optarg)*DEG); break;
    case 'l':
      kozai->set_Omega1(atof(optarg)*DEG); break;
    case 'L':
      kozai->set_Omega2(atof(optarg)*DEG); break;
    case 'i':
      kozai->set_inc(atof(optarg)*DEG); break;
    case 'b':
      kozai->set_r1(atof(optarg)*RSUN); break;
    case 'B':
      kozai->set_r2(atof(optarg)*RSUN); break;
    case 'q':
      kozai->set_quadrupole(true); break;
    case 'o':
      kozai->set_octupole(true); break;
    case 'p':
      kozai->set_pericenter(true); break;
    case 'P':
      kozai->set_outer_pericenter(true); break;
    case 'X':
      kozai->set_cross_lim(true); break;
    case 'r':
      kozai->set_radiation(true); break;
    case 'f':
      PRINT_GSL = true; break;
    case 'I':
      IGNORE_GSL_ERRORS = true; break;
    case 'd':
      t_end = atof(optarg)*YEAR; break;
    case 'D':
      t_step = atof(optarg);
      if (t_step <= 0) break;
      t_step*=YEAR;
      break;
    //----- MODIFIED (Arakida+) -------------------------------------------
    case 300: break;   // --units  : already handled in the pre-scan in main()
    case 301: break;   // --clight : likewise
    case 310: kozai->set_x_rn4(false);   break;
    case 311: kozai->set_x_rn2p5(false); break;
    case 312: kozai->set_x_rn1p5(false); break;
    case 313: kozai->set_lr_rn1p5_exact(true); break;
    case 314: kozai->set_lr_rn2p5_exact(true); break;
    case 315:
      if      (strcmp(optarg,"kst")==0 || strcmp(optarg,"ours")==0) kozai->set_x_rn3p5(1);
      else if (strcmp(optarg,"lr") ==0)                             kozai->set_x_rn3p5(2);
      else if (strcmp(optarg,"lrraw")==0){ kozai->set_x_rn3p5(3); kozai->set_cross_lim(true); }
      else { cout << "--x_rn3p5 takes kst|lr|lrraw\n"; exit(1); }
      break;
    case 316: kozai->set_brown(true);       break;
    case 317: kozai->set_restricted(true);  break;
    case 318: kozai->set_outer_frame(true); break;
    case 319: RAW_OUTPUT = 1;               break;
    case 320:  // --ours : LR's cross terms with our Paper I corrections applied
      kozai->set_cross_lim(true);
      kozai->set_x_rn4(true);
      kozai->set_x_rn2p5(true); kozai->set_lr_rn2p5_exact(true);
      kozai->set_x_rn1p5(false);        // K^(1,3/2) == 0
      kozai->set_x_rn3p5(1);            // K^(1,7/2) = -K_star
      break;
    case 321:  // --lr_full : LR as published, with (A7)-(A11) restored
      kozai->set_cross_lim(true);
      kozai->set_x_rn4(true);
      kozai->set_x_rn2p5(true);  kozai->set_lr_rn2p5_exact(true);
      kozai->set_x_rn1p5(true);  kozai->set_lr_rn1p5_exact(true);
      kozai->set_x_rn3p5(2);
      break;
    case 322: AUDIT_ONLY = 1; break;
    case 'h':
      cout << help;
      exit(1);
    case '?':
      cout << help;
      exit(1);
    default:
      exit(1);
      }

}
