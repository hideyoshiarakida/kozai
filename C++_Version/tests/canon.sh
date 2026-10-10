#!/bin/bash
# Smoke test of the canonical cross terms of Paper II (option --canon) and of the
# environment switches CANON13, CANON14, CANON14Z, CANON14Y (factors multiplying
# K^(1,3), K^(1,4), c_Z Z and c_Y Y; 0 switches a term off).
#
# A. LR's librating system (LR Fig. 3; Paper II Sec. VII A, Fig. 5 top row,
#    blue trajectory): m = (30, 20, 2e7) Msun, a1 = 0.10 AU, a2 = 209.84 AU,
#    e2 = 0.8, l_z = -0.6593, initial e1 = 0.40, 2000 yr, output every 0.05 yr,
#    --quad --peri --outerperi --canon.  Range of e1 and maximum of l_z (= jz)
#    over the output points, for the full canonical model, c_Y Y only, and the
#    electric term K^(1,3) only.  Paper II quotes e1 in [0.226, 0.487],
#    max l_z = -0.6453 (full); e1 in [0.194, 0.511], max l_z = -0.6386 (c_Y Y
#    only); e1 in [0.318, 0.400] with l_z constant to 1e-5 (electric only).
#
# B. LR's Fig. 8 system, orientation A of Paper I (Paper II Sec. VII B,
#    Table VIII): a1 = 0.0713 AU, e1 = 0.808, a2 = 144.694 AU, iota = 37.935 deg,
#    --quad --oct --peri --outerperi, no GW emission.  "Maximum Eccentricity"
#    printed by the code.  Table VIII is for 3e5 yr; e1,max is reached in the
#    first 1e4 yr, so this test integrates 1e4 yr with output every 1e3 yr
#    (about 1 s per run).  With FULL=1 the runs cover 3e5 yr with output every
#    1e4 yr (about 30 s per run), as in the paper.
#
# Every run is made twice: with the corrected code (sign of the (1,7/2) term
# for m1 > m2 fixed in 2026-10, see README.md) and with X72LEGACY=1, which
# restores the submitted code.  Each number is compared (i) with the value
# recorded when this test was set up (tolerance TOL, default 1e-8, to allow
# for other compilers or libraries) and (ii) with the value quoted in Paper II
# (to the digits quoted there): in part A for both codes (the correction does
# not change the quoted digits), in part B (Table VIII) for X72LEGACY=1 only,
# since the correction changes Table VIII in the 4th digit.  The X72LEGACY=1
# values of part A, and of the --canon rows of part B with FULL=1, are
# identical to all printed digits to the research runs from which Paper II
# Sec. VII was made (the "no 3BpN" and "Paper I" rows with FULL=1 agree with
# the research table, which has 6 digits).
#
# Usage: ./canon.sh   (needs ../kozai; run `make` first)

set -u
cd "$(dirname "$0")"
EXE=../kozai
[ -x "$EXE" ] || { echo "canon.sh: $EXE not found; run make first"; exit 1; }
TOL=${TOL:-1e-8}
FULL=${FULL:-0}
fail=0

# check NAME VALUE RECORDED PAPER NDIG : |VALUE-RECORDED| <= TOL and VALUE
# rounded to NDIG decimals equals PAPER (PAPER = "-" skips the second check)
check() {
  local line st
  line=$(awk -v v="$2" -v r="$3" -v p="$4" -v n="$5" -v tol="$TOL" 'BEGIN{
        d = v - r; if (d < 0) d = -d;
        ok1 = (v != "") && (d <= tol);
        ok2 = 1; if (p != "-") ok2 = (sprintf("%." n "f", v) == sprintf("%." n "f", p));
        printf "%16.10f  recorded %16.10f  paper %8s  %s", v, r, p, (ok1 && ok2) ? "ok" : "MISMATCH";
        exit !(ok1 && ok2) }')
  st=$?
  printf '%-36s %s\n' "$1" "$line"
  [ $st -eq 0 ] || fail=1
}

# ---------------------------------------------------------------- part A
LIB="--m1 30 --m2 20 --m3 2e7 --a1 0.1 --a2 209.84 --e1 0.4 --e2 0.8 --g1 90 --g2 282.27
     --omega1 192.5 --omega2 12.5 --inc 136.00127638618591
     --quad --peri --outerperi --canon --time 2000 --dt 0.05 --raw"
# inc = acos(l_z / sqrt(1 - e1^2)) in degrees, l_z = -0.6593, e1 = 0.4

libstats() {   # e1 min, e1 max, l_z max over the --raw output
  env "$@" $EXE $LIB 2>&1 >/dev/null | awk '!/^#/{e=$19; z=$15;
    if (n==0 || e<emin) emin=e; if (n==0 || e>emax) emax=e; if (n==0 || z>zmax) zmax=z;
    if (n==0 || z<zmin) zmin=z; n++}
    END{printf "%.16g %.16g %.16g %.16g\n", emin, emax, zmax, zmax-zmin}'
}

echo "== Paper II, Sec. VII A: librating system, 2000 yr =="
echo "-- corrected code (recorded values; paper digits still agree)"
read a b z dz <<< "$(libstats)"
check "full: e1 min"                        "$a" 0.2256392264 0.226 3
check "full: e1 max"                        "$b" 0.4868826069 0.487 3
check "full: max l_z"                       "$z" -0.6453205444 -0.6453 4
read a b z dz <<< "$(libstats CANON13=0 CANON14Z=0)"
check "c_Y Y only: e1 min"                  "$a" 0.1940174164 0.194 3
check "c_Y Y only: e1 max"                  "$b" 0.5108847169 0.511 3
check "c_Y Y only: max l_z"                 "$z" -0.6386087568 -0.6386 4
read a b z dz <<< "$(libstats CANON14=0)"
check "electric only: e1 min"               "$a" 0.3179513265 0.318 3
check "electric only: e1 max"               "$b" 0.4000583281 0.400 3
check "electric only: range of l_z"         "$dz" 3.1002e-06 - 0
awk -v d="$dz" 'BEGIN{exit !(d < 1e-5)}' || { echo "   range of l_z not < 1e-5"; fail=1; }
echo "-- X72LEGACY=1 (submitted code; values of the research runs)"
read a b z dz <<< "$(libstats X72LEGACY=1)"
check "full: e1 min"                        "$a" 0.2255926272694566 0.226 3
check "full: e1 max"                        "$b" 0.4869301644058858 0.487 3
check "full: max l_z"                       "$z" -0.6453083850058789 -0.6453 4
read a b z dz <<< "$(libstats X72LEGACY=1 CANON13=0 CANON14Z=0)"
check "c_Y Y only: e1 min"                  "$a" 0.1939742134719959 0.194 3
check "c_Y Y only: e1 max"                  "$b" 0.5109001184766297 0.511 3
check "c_Y Y only: max l_z"                 "$z" -0.6385957053507977 -0.6386 4
read a b z dz <<< "$(libstats X72LEGACY=1 CANON14=0)"
check "electric only: e1 min"               "$a" 0.3180174712500569 0.318 3
check "electric only: e1 max"               "$b" 0.4000006553968322 0.400 3

# ---------------------------------------------------------------- part B
if [ "$FULL" = 1 ]; then T=3e5; DT=1e4; else T=1e4; DT=1e3; fi
F8="--m1 30 --m2 20 --m3 2e7 --a1 0.0713 --a2 144.694 --e1 0.808 --e2 0.8 --g1 82.076975
    --g2 159.369837 --omega1 127.304623 --omega2 307.304623 --inc 37.935
    --quad --oct --peri --outerperi --time $T --dt $DT"
# angles of orientation A: numpy default_rng(1016).uniform(0, 360, 3) -> g1, g2, Omega1;
# Omega2 = Omega1 + 180

emax() { env "${@:2}" $EXE $F8 $1 2>/dev/null | awk '/Maximum Eccentricity/{print $NF}'; }

echo
echo "== Paper II, Sec. VII B, Table VIII (orientation A): e1,max over $T yr =="
# corrected code (N = new, recorded) and X72LEGACY=1 (O = old, recorded; P = Table VIII)
if [ "$FULL" = 1 ]; then
  N=(0.8103972854 0.8862220092 0.9091603514 0.8125709596 0.8123999171 0.9067373599 0.9091632121 0.8123994368)
  O=(0.8103972854 0.8863644115 0.9092392631 0.8115367549 0.811458831  0.9068239385 0.9092400949 0.8114586755)
else
  N=(0.810395085  0.8862220047 0.9091579712 0.8125641186 0.8123990757 0.9067380862 0.9091608312 0.8123990134)
  O=(0.810395085  0.8863644097 0.9092388688 0.8115111369 0.8114579693 0.9068234241 0.9092437052 0.8114570293)
fi
P=(0.8104 0.8864 0.9092 0.8115 0.8115 0.9068 - -)
LBL=("no 3BpN" "Paper I (--ours)" "this work (--canon)" "  electric only (CANON14=0)"
     "  c_Z Z only (CANON13=0 CANON14Y=0)" "  c_Y Y only (CANON13=0 CANON14Z=0)"
     "  K^(1,4) only (CANON13=0)" "  without c_Y Y (CANON14Y=0)")
ARG=("" --ours --canon --canon --canon --canon --canon --canon)
ENV=("" "" "" "CANON14=0" "CANON13=0 CANON14Y=0" "CANON13=0 CANON14Z=0" "CANON13=0" "CANON14Y=0")
echo "-- corrected code (recorded values)"
for k in 0 1 2 3 4 5 6 7; do
  check "${LBL[$k]}" "$(emax "${ARG[$k]}" ${ENV[$k]})" ${N[$k]} - 0
done
echo "-- X72LEGACY=1 (submitted code; compared with Table VIII of the submitted Paper II)"
for k in 0 1 2 3 4 5 6 7; do
  check "${LBL[$k]}" "$(emax "${ARG[$k]}" X72LEGACY=1 ${ENV[$k]})" ${O[$k]} ${P[$k]} 4
done

echo
if [ $fail -eq 0 ]; then echo "canon: PASS"; exit 0; else echo "canon: FAIL"; exit 1; fi
