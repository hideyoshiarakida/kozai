# Paper I — C++ secular three-body code

Modified version of the secular three-body code of C. L. Rodriguez with the
relativistic cross terms of Lim & Rodriguez (LR; Phys. Rev. D 102, 064033, 2020), used for
H. Arakida, *Post-Newtonian secular dynamics of hierarchical triples. I* (Paper I):

* Table VII (Sec. VIII, "Tests"): the normalization-free ratios ρ_ι, ρ_Ω of the (1,7/2)
  rates — reproduced by `tests/audit.sh`;
* Sec. IX (ZLK evolution, LR's Figs. 3–8): the integrations of the models "no 3BpN",
  "LR", "LR without Eq. (4.5)" and "this work" (options in the table below; the drivers
  and figure scripts are in the code release of Papers I and II, `PaperI/python/zlk`).

## Provenance

This folder is the branch `arakida-crossterms` of a fork of C. L. Rodriguez's repository
<https://github.com/carlrodriguez/kozai> (branch `crossterms`, commit `e7f5203`, whose tree is
identical to the upstream commit below); the history of the fork records the changes commit
by commit.

Upstream: <https://github.com/halstonblim/kozai>, branch `Sept2019_crossterms`, commit
`1393251748fca1942f16c8639ed81b54286615c7` (2020-09-28), directory `C++_Version/`.
`vector.h` is unchanged; `kozai.C` and `kozai.h` are modified. `upstream.patch` is the
`diff -ruN` of the three files from the upstream commit to this version; applied to a
clone of that commit with

    git clone https://github.com/halstonblim/kozai.git && cd kozai
    git checkout 1393251748fca1942f16c8639ed81b54286615c7
    patch -p1 < upstream.patch          # (path to this file)

it reproduces `kozai.C`, `kozai.h` and `vector.h` of this folder byte for byte (checked by
`tests/regress.sh`).

The upstream README says: *"You are welcome to use either code in any way you like, though
we ask that if you use any results from it in a paper that you cite Rodriguez and Antonini
(2018, https://arxiv.org/abs/1805.08212) and Antonini et al (2017,
https://arxiv.org/abs/1711.07142), and Lim and Rodriguez (2020,
https://arxiv.org/abs/2001.03654)"*. Please cite these papers, and Paper I for the
modifications.

## Changes with respect to upstream

All changes are marked `MODIFIED (Arakida+)` in the source. Without the new options the
code is bit-for-bit the upstream code (`tests/regress.sh`), except for the collision
radius (item 2).

1. **Unit system.** `G`, `c`, `YEAR`, `AU`, `RSUN`, `MSUN` are no longer `const`.
   `--units code` switches to G = 1, M_sun = AU = yr = 1 (the units of the Python runs of
   the papers); `--clight` sets the speed of light (in code units default 1e6, i.e. the
   Newtonian limit). Both are read before all other options. In code units the default
   triple is the configuration of Paper II, Sec. VI (m = 1.5, 0.5, 7; a1 = 1, a2 = 14,
   e2 = 0.6, j_z = 0.4, e1 = 0.24451925, ω1 = 0.30 rad), with `--outer_frame` and
   `--restricted`.
2. **Collision radius (bug fix).** Upstream sets the default radius to
   `2*1476.*m/MSUN`, i.e. 2Gm/c² with 1476 in metres while all lengths are in cm (100
   times too small). This version uses `2.*G*m/(c*c)`. Only the printed radius changes in
   the regression runs; it would matter only for collisions.
3. **`get_L1()` (bug fix).** Upstream reads `y[18]` (an 18-dimensional state vector with
   spins) while `DIMENSION` is 13; now `y[12]` (= a1). `get_L1()` is not called in the
   upstream code, so this did not affect any result.
4. **Cross terms grade by grade.** The cross-term block of `--cross_lim` is split into its
   grades, each with its own switch: (1,4) [LR Eq. (4.6) etc.], (1,5/2) [LR Eq. (4.4), de
   Sitter], (1,3/2) [LR Eq. (4.5), shown in Paper I to vanish, K^(1,3/2) = 0]. With
   `--cross_lim` alone all are on and the output is upstream's.
   `--lr_rn1p5` / `--lr_rn2p5` replace the mass prefactors hard-wired in the port by
   LR's published ones (m m3/√M instead of m √M for (1,3/2); (4m+3m3) m3/(2√M) instead of
   the test-particle limit 3M^{3/2}/2 for (1,5/2)).
5. **Grade (1,7/2)**, absent from the port: `--x_rn3p5 kst` (K^(1,7/2) = −N Ξ/ℓ2⁵ with
   N ∝ (m1 − m2) and e1 the eccentricity vector of x1 − x2, the magnetic-quadrupole term of
   KST Eq. (35), as a Hamiltonian flow through the Milankovitch equations; sign corrected
   for m1 > m2 in 2026-10, see below), `--x_rn3p5 lr` (LR
   (A7)–(A11) reconstructed as a Hamiltonian), `--x_rn3p5 lrraw` (LR (A7), (A9), (A11)
   transcribed verbatim; non-canonical; LR give no e1 and ω1 rates at this grade, so these
   are zero).
6. **Brown term** (`--brown`): the ω2-independent part of K^(0,9/2) (Newtonian,
   quadrupole squared, discussed in Paper II, Sec. VI), C ε21 F21, as a Hamiltonian flow.
7. **Frames and the restricted problem.** `--outer_frame` (reference plane = outer orbit
   plane, so that the printed Ω1 is the Ω1 of the papers and ι1 = ι), `--restricted`
   (outer orbit frozen).
8. **Output.** `--raw` (full state at 17 digits); `--audit` (rates of each grade at t = 0
   and the ratios ρ_ι, ρ_Ω of Table VII; then exit). New accessors `get_jz`,
   `get_omega1_mutual`, `get_Omega1_acc` (atan2 version of `get_Omega1`, which uses acos
   and loses digits near 0 and π; `get_Omega1` is kept for the ordinary output), `get_n0`,
   `get_n3`, `get_eps21`.
9. **Shorthands** `--ours` and `--lr_full` (below).

Not changed: the port contains only the grades (1,4), (1,5/2) and (1,3/2) of LR's
`domega1dtcross` (LR's Rn3 and Rn2 pieces are absent upstream); as in upstream, the cross
terms act on the inner orbit only (no back-reaction on the outer orbit).

## Correction of the (1,7/2) sign for m1 > m2 (2026-10)

**What was wrong.** The octupole block of the (upstream) code is written with |m1 − m2|,
i.e. its e1 is the eccentricity vector of x_light − x_heavy. The (1,7/2) block
(`--x_rn3p5 kst`, used by `--ours` and by `--canon` of Paper II) used (m1 − m2) with the same
e1, whereas K^(1,7/2) = −N Ξ/ℓ2⁵ with N ∝ (m1 − m2) holds for e1 of x1 − x2. For m1 > m2 the
(1,7/2) term therefore entered with the wrong sign relative to the octupole (for m1 < m2 it
was correct). At the state of `tests/x72sign.sh` the physical ratio
Ġ1(1,7/2)/Ġ1(0,4) is −3.4913e-3; the submitted code gave +3.4913e-3. The block now uses
s N = +N′|m1 − m2| (N′ = N/(m1 − m2)) with the code's e1; the LR reconstruction
(`--x_rn3p5 lr`) and the verbatim LR rates (`lrraw`) are unchanged.

**Effect.** The normalization-free ratios of Table VII are unchanged (`tests/audit.sh`).
In the integrations of Sec. IX (m1 = 30 > m2 = 20) the model "this work" (`--ours`)
changes: Table IX of Paper I and Table VIII of Paper II change in the 4th digit (e.g.
orientation A, e1,max over 3e5 yr: 0.8864 → 0.8862 for `--ours`, 0.9092 → 0.9092
[0.90924 → 0.90916] for `--canon`); the merger time of the Paper I model in LR's Fig. 8
system changes from 1.433e7 to 1.450e7 yr; the population statistics (orientation
sample) are unchanged.

**Submitted behaviour.** The environment variable `X72LEGACY=1` restores the sign of the
submitted version (e.g. `X72LEGACY=1 ./kozai ... --ours`); it is read once.

## Options

Upstream options (`--m1 --m2 --m3 --a1 --a2 --e1 --e2 --g1 --g2 --omega1 --omega2 --inc
--rad1 --rad2 --time --dt --quad --oct --peri --outerperi --cross_lim --rad --ignore_gsl
--print_gsl --help`) are unchanged; angles are in degrees. New options:

| option | meaning |
|---|---|
| `--units code` | G = 1, M_sun = AU = yr = 1; code-unit default triple (see above) |
| `--clight X` | speed of light in the current units |
| `--raw` | stderr: `t, j1(3), e1(3), j2(3), e2(3), a1, jz, iota, omega1, Omega1, e1, e2` at 17 digits |
| `--outer_frame` | reference plane = outer orbit plane (then ω2 is set by `--omega2 180 --g2 <ω2>`) |
| `--restricted` | freeze the outer orbit |
| `--no_rn4` | drop grade (1,4) |
| `--no_rn2p5` | drop grade (1,5/2) (de Sitter, LR Eq. (4.4)) |
| `--no_rn1p5` | drop grade (1,3/2) (LR Eq. (4.5)) |
| `--lr_rn1p5`, `--lr_rn2p5` | LR's published mass prefactors for (1,3/2), (1,5/2) |
| `--x_rn3p5 kst\|lr\|lrraw` | add grade (1,7/2): KST/ours, LR reconstructed, LR verbatim (`lrraw` implies `--cross_lim`); `X72LEGACY=1` in the environment restores the submitted sign of `kst` |
| `--brown` | add the Brown term (ω2-independent part of K^(0,9/2)) |
| `--audit` | evaluate the rates of each grade once at t = 0, print Ġ1, ι̇, Ω̇1, ω̇1, ρ_ι, ρ_Ω; exit |
| `--ours` | = `--cross_lim`, (1,4) on, (1,5/2) on with LR's mass prefactor, (1,3/2) off, (1,7/2) = KST |
| `--lr_full` | = `--cross_lim`, all LR terms with LR's published prefactors, (1,7/2) = LR reconstructed |

Models of Paper I, Sec. IX (all with `--quad --peri --outerperi`, and `--oct`, `--rad`
where stated): "no 3BpN": no further option; "LR": `--cross_lim`; "LR without Eq. (4.5)":
`--cross_lim --no_rn1p5`; "this work": `--ours`.

Output: summary to stdout (setup, "Maximum Eccentricity"), trajectory to stderr (as
upstream).

## Build and test

Requires a C++11 compiler and GSL (tested with g++ 13.3 and GSL 2.7.1 on Linux x86-64).

    make            # builds ./kozai (uses gsl-config if present, else -lgsl -lgslcblas -lm)
    make test       # tests/audit.sh, tests/x72sign.sh and tests/regress.sh
    make clean      # removes kozai, kozai_orig and tests/upstream

The default flags are `-O3 -ffp-contract=off` (the latter prevents fused multiply-adds,
so that bit-for-bit comparisons are portable). `CXX` and `CXXFLAGS` can be overridden.

| test | what | runtime | expected output |
|---|---|---|---|
| `tests/audit.sh` | Paper I **Table VII**: `kozai --audit` at ι = 1.1, ω1 = 0.7, ω2 = 1.3 rad; compares the (1,7/2) rows | < 1 s | ρ_ι, ρ_Ω = 1.46406, 2.59997 (KST, `--x_rn3p5 kst`) and 2.74027, 5.39947 (LR (A7), (A9), (A11), `--x_rn3p5 lrraw`); ratios 1.872, 2.077; `audit: PASS` |
| `tests/x72sign.sh` | sign of the (1,7/2) term: `--audit` at `--units code --m1 30 --m2 20 --m3 40 --a1 1 --a2 12 --e1 0.4 --e2 0.35 --g1 40 --g2 70 --omega1 0 --omega2 180 --inc 65 --clight 300 --quad --oct --x_rn3p5 kst`, ratio of the Ġ1 of the rows (1,7/2) and (0,4) | < 1 s | −3.4913e-3 (rel. tol. 1e-4) for m1, m2 = 30, 20 and 20, 30; +3.4913e-3 with `X72LEGACY=1` for 30, 20; `x72sign: PASS` |
| `tests/regress.sh` | clones upstream at commit 1393251 into `tests/upstream`, checks that `upstream.patch` reproduces the sources, builds `kozai_orig`, runs 3 option sets without the new options with both executables | 5 s + clone | stderr (trajectory) bit-identical for `--quad --oct --peri --outerperi --cross_lim`, `--quad --oct`, `--quad --oct --peri --outerperi --rad --cross_lim`; stdout identical apart from the radius lines (item 2); `regress: PASS` |

`tests/regress.sh` needs `git` and network access to github.com (or an existing clone in
`tests/upstream`). The output of `make test` from a clean checkout is in
`tests/expected_output.txt`.

How the ratios are produced: `--audit` evaluates the right-hand side once at the initial
condition with one grade switched on at a time and obtains Ġ1, ι̇ and Ω̇1 by a centred
difference of the elements along the flow (`run_audit()` in `kozai.C`). The ratios are
independent of the overall normalization, the sign conventions and the masses.

## Files

| file | content |
|---|---|
| `kozai.C`, `kozai.h` | the modified code |
| `vector.h` | 3-vector class (unchanged from upstream) |
| `upstream.patch` | `diff -ruN` upstream → this version (apply with `patch -p1` in the upstream clone) |
| `Makefile` | `make`, `make test`, `make clean` |
| `tests/audit.sh`, `tests/x72sign.sh`, `tests/regress.sh` | tests (relative paths only) |
| `tests/expected_output.txt` | log of `make test` |
