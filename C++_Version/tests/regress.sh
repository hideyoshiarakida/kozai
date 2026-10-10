#!/bin/bash
# Regression test against the unmodified upstream code.
#
#  1. clones https://github.com/halstonblim/kozai (commit 1393251, branch
#     Sept2019_crossterms) into tests/upstream (reused if already there);
#  2. checks that ../upstream.patch applied with `patch -p1` to that commit
#     reproduces ../kozai.C, ../kozai.h and ../vector.h byte for byte;
#  3. builds the unmodified code as ../kozai_orig;
#  4. runs three option sets that use none of the new options with ../kozai
#     and ../kozai_orig and requires the integrated trajectory (stderr) to be
#     bit-identical;
#  5. rebuilds the Paper I version of the code from ../kozai.C, ../kozai.h by
#     applying ../PaperI_to_PaperII.patch in reverse, and checks that without
#     --canon the Paper I and Paper II executables give bit-identical output
#     (stdout and stderr) for five option sets, including the options of
#     Paper I (--ours, --brown, --x_rn3p5, --units code, --raw).
#
# Known, intended difference: upstream sets the default collision radius to
# 2*1476*m/MSUN, i.e. 2Gm/c^2 with 1476 in metres while all lengths are in cm
# (100 times too small); this version uses 2Gm/c^2.  Only the "radius" lines
# of stdout can therefore differ; they are shown but not counted as failures.
#
# Needs git and network access (or an existing tests/upstream clone).
# Usage: ./regress.sh   (run `make` first; `make clean` removes the clone)

set -u
cd "$(dirname "$0")"
URL=https://github.com/halstonblim/kozai.git
COMMIT=1393251748fca1942f16c8639ed81b54286615c7
CXX=${CXX:-g++}
CXXFLAGS=${CXXFLAGS:--O3 -ffp-contract=off}
GSL_CFLAGS=$(gsl-config --cflags 2>/dev/null)
GSL_LIBS=$(gsl-config --libs 2>/dev/null || echo -lgsl -lgslcblas -lm)
EXE=../kozai
[ -x "$EXE" ] || { echo "regress.sh: $EXE not found; run make first"; exit 1; }
fail=0

echo "== regression test against upstream commit ${COMMIT:0:7} =="

# 1. upstream
if [ ! -d upstream/.git ]; then
  git clone -q "$URL" upstream || { echo "regress: FAIL (cannot clone $URL)"; exit 1; }
fi
git -C upstream checkout -q "$COMMIT" 2>/dev/null || {
  git -C upstream fetch -q origin "$COMMIT" && git -C upstream checkout -q "$COMMIT"; } ||
  { echo "regress: FAIL (commit $COMMIT not found)"; exit 1; }
echo "upstream: $(git -C upstream log -1 --format='%h %ad %s' --date=short)"

# 2. the patch reproduces our sources
rm -rf work && mkdir -p work/patched
cp -r upstream/C++_Version work/patched/
if (cd work/patched && patch -s -p1 < ../../../upstream.patch); then
  for f in kozai.C kozai.h vector.h; do
    if cmp -s "work/patched/C++_Version/$f" "../$f"; then
      echo "upstream.patch -> $f : identical"
    else
      echo "upstream.patch -> $f : DIFFERENT"; fail=1
    fi
  done
else
  echo "upstream.patch does not apply"; fail=1
fi

# 3. build the unmodified code
$CXX $CXXFLAGS $GSL_CFLAGS -o ../kozai_orig upstream/C++_Version/kozai.C $GSL_LIBS ||
  { echo "regress: FAIL (cannot build kozai_orig)"; exit 1; }

# 4. trajectories
OPTS=("--quad --oct --peri --outerperi --cross_lim --time 3e7 --dt 1e6"
      "--quad --oct --time 1e7 --dt 5e5"
      "--quad --oct --peri --outerperi --rad --cross_lim --time 3e7 --dt 1e6")
for o in "${OPTS[@]}"; do
  ../kozai_orig $o > work/orig.out 2> work/orig.err
  $EXE          $o > work/new.out  2> work/new.err
  n=$(wc -l < work/new.err)
  if cmp -s work/orig.err work/new.err; then
    echo "stderr IDENTICAL ($n lines): $o"
  else
    echo "stderr DIFFERENT ($n lines): $o"; fail=1
  fi
  if diff -q <(grep -v 'radius' work/orig.out) <(grep -v 'radius' work/new.out) > /dev/null; then
    echo "   stdout identical apart from the radius lines" \
         "(radius 1 [RSUN]: upstream $(awk '/radius 1/{print $NF}' work/orig.out)," \
         "this version $(awk '/radius 1/{print $NF}' work/new.out))"
  else
    echo "   stdout differs beyond the radius lines"; fail=1
  fi
done

# 5. Paper I behaviour unchanged when --canon is not given
mkdir -p work/paperI
cp ../kozai.C ../kozai.h ../vector.h work/paperI/
if (cd work/paperI && patch -s -R -p1 < ../../../PaperI_to_PaperII.patch) &&
   $CXX $CXXFLAGS $GSL_CFLAGS -o work/paperI/kozai_paperI work/paperI/kozai.C $GSL_LIBS; then
  echo "Paper I version rebuilt from PaperI_to_PaperII.patch (reverse)"
  F8="--m1 30 --m2 20 --m3 2e7 --a1 0.0713 --a2 144.694 --e1 0.808 --e2 0.8 --g1 82.076975
      --g2 159.369837 --omega1 127.304623 --omega2 307.304623 --inc 37.935"
  LIB="--m1 30 --m2 20 --m3 2e7 --a1 0.1 --a2 209.84 --e1 0.4 --e2 0.8 --g1 90 --g2 282.27
       --omega1 192.5 --omega2 12.5 --inc 136.00128"
  OPTS1=("--quad --oct --peri --outerperi --cross_lim --time 3e7 --dt 1e6"
         "--quad --oct --peri --outerperi --ours --time 3e7 --dt 1e6"
         "$F8 --quad --oct --peri --outerperi --ours --time 1e4 --dt 100 --raw"
         "$LIB --quad --peri --outerperi --cross_lim --time 2000 --dt 1 --raw"
         "--units code --quad --brown --x_rn3p5 lrraw --clight 300 --time 2000 --dt 10 --raw")
  for o in "${OPTS1[@]}"; do
    work/paperI/kozai_paperI $o > work/p1.out 2> work/p1.err
    $EXE                     $o > work/p2.out 2> work/p2.err
    n=$(wc -l < work/p2.err)
    if cmp -s work/p1.err work/p2.err && cmp -s work/p1.out work/p2.out; then
      echo "Paper I vs Paper II IDENTICAL ($n lines): $(echo $o | cut -c1-90)..."
    else
      echo "Paper I vs Paper II DIFFERENT ($n lines): $o"; fail=1
    fi
  done
else
  echo "cannot rebuild the Paper I version from PaperI_to_PaperII.patch"; fail=1
fi

rm -rf work
if [ $fail -eq 0 ]; then echo "regress: PASS"; exit 0; else echo "regress: FAIL"; exit 1; fi
