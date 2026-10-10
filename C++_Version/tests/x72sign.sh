#!/bin/bash
# Sign of the grade-(1,7/2) term (--x_rn3p5 kst, used by --ours and --canon)
# relative to the Newtonian octupole (correction of 2026-10, see README.md).
#
# The octupole block of the code uses |m1-m2|, i.e. its e1 is the eccentricity
# vector of x_light - x_heavy, while K^(1,7/2) = -N Xi / l2^5 with N ~ (m1-m2)
# is written for e1 of x1 - x2.  The submitted versions used (m1-m2) with the
# code's e1, so that for m1 > m2 the (1,7/2) term had the wrong sign relative to
# the octupole.  This test runs --audit at
#
#   --units code --m1 30 --m2 20 --m3 40 --a1 1 --a2 12 --e1 0.4 --e2 0.35
#   --g1 40 --g2 70 --omega1 0 --omega2 180 --inc 65 --clight 300
#   --quad --oct --x_rn3p5 kst
#
# and checks the ratio Gdot(1,7/2) / Gdot(0,4) of the audit rows
# "(1,7/2) Rn3p5  ours = -K_star" and "(0,4) Newtonian octupole":
#   m1 = 30, m2 = 20:                  -3.4913e-3   (relative tolerance 1e-4)
#   m1 = 20, m2 = 30 (mirrored):       -3.4913e-3   (old and new code agree)
#   m1 = 30, m2 = 20, X72LEGACY=1:     +3.4913e-3   (sign of the submitted code)
#
# Usage: ./x72sign.sh   (needs ../kozai; run `make` first)

set -u
cd "$(dirname "$0")"
EXE=../kozai
[ -x "$EXE" ] || { echo "x72sign.sh: $EXE not found; run make first"; exit 1; }
REF=-3.4913e-3
RTOL=1e-4
fail=0

ratio() {   # ratio <legacy 0|1> <m1> <m2>
  X72LEGACY=$1 $EXE --units code --m1 $2 --m2 $3 --m3 40 --a1 1 --a2 12 --e1 0.4 \
    --e2 0.35 --g1 40 --g2 70 --omega1 0 --omega2 180 --inc 65 --clight 300 \
    --quad --oct --x_rn3p5 kst --audit 2>/dev/null |
  awk 'index($0, "(0,4)") == 1 && index($0, "Newtonian octupole") {g4=$(NF-5)}
       index($0, "(1,7/2) Rn3p5  ours = -K_star") == 1 {g72=$(NF-5)}
       END{if (g4 == "" || g72 == "") print "nan"; else printf "%.6e %.8e %.8e\n", g72/g4, g72, g4}'
}

check() {   # check <label> <legacy> <m1> <m2> <expected>
  local out r
  out=$(ratio $2 $3 $4); r=${out%% *}
  if awk -v r="$r" -v e="$5" -v t="$RTOL" 'BEGIN{d=(r-e)/e; if (d<0) d=-d; exit !(r!="nan" && d<=t)}'; then
    st=ok
  else
    st=MISMATCH; fail=1
  fi
  printf '%-36s Gdot(1,7/2)/Gdot(0,4) = %13s  expected %11s  %s\n' "$1" "$r" "$5" "$st"
}

echo "== sign of the (1,7/2) term relative to the octupole (kozai --audit) =="
check "m1=30, m2=20"                 0 30 20 $REF
check "m1=20, m2=30 (mirrored)"      0 20 30 $REF
check "m1=30, m2=20, X72LEGACY=1"    1 30 20 3.4913e-3
check "m1=20, m2=30, X72LEGACY=1"    1 20 30 $REF

if [ $fail -eq 0 ]; then echo "x72sign: PASS"; exit 0; else echo "x72sign: FAIL"; exit 1; fi
