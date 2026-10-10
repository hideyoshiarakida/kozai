#!/bin/bash
# Paper I, Table VII: normalization-free ratios of the (1,7/2) rates,
#
#     rho_iota  = G_1 iotadot   / Gdot_1 ,   rho_Omega = G_1 Omegadot_1 / Gdot_1 ,
#
# at iota = 1.1, omega_1 = 0.7, omega_2 = 1.3 (rad), for
#   - the KST/our Hamiltonian K^(1,7/2) = -K_star     (--x_rn3p5 kst)
#   - LR's rates (A7), (A9), (A11), verbatim         (--x_rn3p5 lrraw)
#
# The ratios are produced by the --audit option of ../kozai: the right-hand side
# is evaluated once at the initial condition with one grade switched on at a
# time, and the rates of G_1, iota and Omega_1 are obtained by a centred
# difference along the flow (see run_audit() in kozai.C).  The two (1,7/2) rows
# of the audit table are compared with the values printed in Table VII.
#
# Angles are given in degrees: 63.0253... = 1.1 rad, 40.1070... = 0.7 rad,
# 74.4845... = 1.3 rad; with --outer_frame, omega_2 is set by --g2 and
# --omega2 180 (LR's node convention Omega_2 = Omega_1 + pi).  The masses and
# semi-major axes do not enter the ratios (Paper I, Sec. on the (1,7/2) term).
#
# The Paper II executable prints one more row, "(1,3)+(1,4) canonical [II]"
# (the --canon terms evaluated alone); it does not enter this test.
#
# Usage: ./audit.sh   (from any directory; needs ../kozai, run `make` first)

set -u
cd "$(dirname "$0")"
EXE=../kozai
[ -x "$EXE" ] || { echo "audit.sh: $EXE not found; run make first"; exit 1; }

out=$($EXE --units code --outer_frame --restricted --audit \
      --m1 1.5 --m2 0.5 --m3 7.0 --a1 1.0 --a2 14.0 --e1 0.4 --e2 0.6 \
      --inc 63.02535746439056 --g1 40.10704565915762 --g2 74.48451336700701 \
      --omega1 0 --omega2 180 --clight 300 2>/dev/null)

echo "== Paper I, Table VII: (1,7/2) rate ratios (kozai --audit) =="
echo "$out" | sed -n '/=== term audit/,/^ *$/p' | sed '/^ *$/d'
echo "$out" | grep -E '^grade|^\(1,7/2\)'

kst=$(echo "$out" | grep -F '(1,7/2) Rn3p5  ours = -K_star' | awk '{print $(NF-1), $NF}')
lr=$(echo "$out"  | grep -F '(1,7/2) Rn3p5  LR (A7)(A9)(A11)' | awk '{print $(NF-1), $NF}')

echo
printf '%-26s %10s %10s   %s\n' "" "rho_iota" "rho_Omega" "Table VII"
printf '%-26s %10s %10s   %s\n' "KST Eq. (35) (ours)" $kst "1.46406 2.59997"
printf '%-26s %10s %10s   %s\n' "LR (A7), (A9), (A11)" $lr "2.74027 5.39947"
echo "$kst $lr" | awk '{printf "%-26s %10.3f %10.3f   %s\n", "ratio LR/KST", $3/$1, $4/$2, "1.872 2.077"}'

if [ "$kst" = "1.46406 2.59997" ] && [ "$lr" = "2.74027 5.39947" ]; then
  echo "audit: PASS"
  exit 0
else
  echo "audit: FAIL (expected 1.46406 2.59997 / 2.74027 5.39947)"
  exit 1
fi
