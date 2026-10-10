# Paper II — C++ secular three-body code with the canonical cross terms

The secular three-body code of Paper I (a modified version of C. L. Rodriguez's code with
the cross terms of Lim & Rodriguez, LR, Phys. Rev. D 102, 064033, 2020) with the canonical
1PN × quadrupole terms of H. Arakida, *Post-Newtonian secular dynamics of hierarchical
triples. II* (Paper II) added as the option `--canon`. It was used for the orbit
integrations of Paper II, Sec. VII (Figs. 5–8, Tables VII and VIII); the drivers and figure
scripts are in the code release of Papers I and II, `PaperII/python/orbit`.

Without `--canon` the code behaves exactly as the Paper I version (tag `paperI-v1`),
bit for bit (checked by `tests/regress.sh`). Everything said in `C++_Version/README.md` at tag `paperI-v1`
about the other options and the changes with respect to upstream applies here; the
essentials are repeated below so that this folder is self-contained.

## Provenance

This folder is the branch `arakida-crossterms` of a fork of C. L. Rodriguez's repository
<https://github.com/carlrodriguez/kozai> (branch `crossterms`, commit `e7f5203`, whose tree is
identical to the upstream commit below); the history of the fork records the changes commit
by commit, and the Paper I version is tagged `paperI-v1`.

Upstream: <https://github.com/halstonblim/kozai>, branch `Sept2019_crossterms`, commit
`1393251748fca1942f16c8639ed81b54286615c7` (2020-09-28), directory `C++_Version/`.

* `upstream.patch`: `diff -ruN` of `kozai.C`, `kozai.h`, `vector.h` from the upstream
  commit to this version; `patch -p1 < upstream.patch` in a clone of the upstream commit
  reproduces the three files of this folder byte for byte.
* `PaperI_to_PaperII.patch`: `diff -ruN` from the Paper I version (tag `paperI-v1`, `C++_Version/kozai.C`,
  `kozai.h`) to this version; apply with `patch -p1` in `C++_Version/` of a checkout of `paperI-v1`
  (`patch -R -p1` in a copy of this folder gives back the Paper I version).
* `vector.h` is the upstream file, unchanged.

The upstream README says: *"You are welcome to use either code in any way you like, though
we ask that if you use any results from it in a paper that you cite Rodriguez and Antonini
(2018, https://arxiv.org/abs/1805.08212) and Antonini et al (2017,
https://arxiv.org/abs/1711.07142), and Lim and Rodriguez (2020,
https://arxiv.org/abs/2001.03654)"*. Please cite these papers, and Papers I and II for the
modifications.

## The option `--canon`

`--canon` is a shorthand for the model "this work" of Paper II, Table VII:

| setting | term | Paper II |
|---|---|---|
| `--cross_lim` with LR's (1,4) block off (`--no_rn4`) | LR's (1,4) terms replaced by ours | Sec. V C |
| (1,5/2) on with LR's published mass prefactor (`--lr_rn2p5`) | de Sitter term, prefactor ∝ (4m + 3m3) m3 | Eq. (38) |
| (1,3/2) off (`--no_rn1p5`) | LR Eq. (4.5); K^(1,3/2) = 0 | Paper I |
| (1,7/2) = KST (`--x_rn3p5 kst`) | magnetic quadrupole K^(1,7/2) = K_MQ (sign for m1 > m2 corrected in 2026-10, below) | Eqs. (45)–(47) |
| new: K^(1,3) | electric 1PN × quadrupole term, F0, F2 | Eqs. (42)–(44), vector form Eq. (D4) |
| new: K^(1,4)_quad = U4 [c_Z Z + c_Y Y] | next-to-leading term (quadrupole type) | Eqs. (50)–(53), Z of Eq. (40) |

The new terms are coded in `rhs()` in `kozai.C` (block "Paper II canonical grades (1,3) and
(1,4)") as a Hamiltonian K(**j**, **e**) = K^(1,3) + K^(1,4)_quad of the inner orbit, with

* K^(1,3) = U3 [F0 (1 − 3 j_z²/j²) + (F2/e1²)(e1² − e1² j_z²/j² − 2 e_z²)],
  U3 = G² m1 m2 m3 a1 /(c² a2³ ℓ2³) — Eq. (D4) with F0, F2 of Eqs. (43), (44);
  F2/e1² is evaluated with a series for e1² < 1e-4 (regular at e1 = 0);
* K^(1,4)_quad = U4 [c_Z Z + c_Y Y], U4 = G² m1 m2 m3 a1² /(c² a2⁴ ℓ2⁵),
  Z = −2 + 12 e1² + 6 (**j**·ẑ2)² − 30 (**e**·ẑ2)² [Eq. (40)],
  Y = 10[(**e**·û2)² − (**e**·v̂2)²] − 2[(**j**·û2)² − (**j**·v̂2)²] [Eq. (51)],
  c_Z, c_Y of Eqs. (52), (53); the octupole-type part of K^(1,4) (∝ m1 − m2) is not included;
* the flow is given by the Milankovitch equations, Eq. (D5),
  d**j**/dt = −(**j** × ∇_j K + **e** × ∇_e K)/L1, d**e**/dt = −(**j** × ∇_e K + **e** × ∇_j K)/L1,
  with the gradients taken by central differences (step 1e-6); as for all cross terms of
  the code, the terms act on the inner orbit only.

## Correction of the (1,7/2) sign for m1 > m2 (2026-10)

**What was wrong.** The octupole block of the (upstream) code is written with |m1 − m2|,
i.e. its e1 is the eccentricity vector of x_light − x_heavy. The (1,7/2) block
(`--x_rn3p5 kst`, used by `--ours` and `--canon`) used (m1 − m2) with the same e1, whereas
K^(1,7/2) = −N Ξ/ℓ2⁵ with N ∝ (m1 − m2) holds for e1 of x1 − x2. For m1 > m2 the (1,7/2)
term therefore entered with the wrong sign relative to the octupole (for m1 < m2 it was
correct). At the state of `tests/x72sign.sh` the physical ratio Ġ1(1,7/2)/Ġ1(0,4) is
−3.4913e-3; the submitted code gave +3.4913e-3. The block now uses s N = +N′|m1 − m2|
(N′ = N/(m1 − m2)) with the code's e1; the LR reconstruction (`--x_rn3p5 lr`), the verbatim
LR rates (`lrraw`) and the canonical (1,3), (1,4) terms are unchanged. The same correction
is in the Paper I version, so `PaperI_to_PaperII.patch` is unchanged in substance.

**Effect.** The ratios of Paper I, Table VII are unchanged (`tests/audit.sh`). For LR's
triples (m1 = 30 > m2 = 20) the models "Paper I" (`--ours`) and "this work" (`--canon`)
change: Table VIII of Paper II and Table IX of Paper I change in the 4th digit
(orientation A, e1,max over 3e5 yr, submitted → corrected: Paper I 0.88636 → 0.88622; this
work 0.90924 → 0.90916; electric only 0.81154 → 0.81257; c_Z Z only 0.81146 → 0.81240;
c_Y Y only 0.90682 → 0.90674); the merger time of the Paper I model in LR's Fig. 8 system
changes from 1.433e7 to 1.450e7 yr; the population statistics (orientation sample) are
unchanged. The numbers quoted in Sec. VII A for the librating system are unchanged at the
digits quoted.

**Submitted behaviour.** The environment variable `X72LEGACY=1` restores the sign of the
submitted version (e.g. `X72LEGACY=1 ./kozai ... --canon`); it is read once.

### Environment switches

For the diagnostics of Paper II, Sec. VII ("switched the canonical terms on one at a time")
the terms are multiplied by factors read from the environment (default 1; read once, at
the first evaluation):

| variable | multiplies | Paper II |
|---|---|---|
| `CANON13` | K^(1,3) | electric term, Eqs. (42)–(44) |
| `CANON14` | K^(1,4)_quad (both parts) | Eq. (50) |
| `CANON14Z` | c_Z Z (ZLK-like part) | Eqs. (40), (52) |
| `CANON14Y` | c_Y Y (ω2-dependent part) | Eqs. (51), (53) |

The rows of Table VIII and the statements of Sec. VII A use: electric term only:
`CANON14=0`; c_Z Z only: `CANON13=0 CANON14Y=0`; c_Y Y only: `CANON13=0 CANON14Z=0`
(the de Sitter and magnetic-quadrupole terms stay on). Example:

    CANON13=0 CANON14Z=0 ./kozai --m1 30 --m2 20 --m3 2e7 --a1 0.0713 --a2 144.694 \
        --e1 0.808 --e2 0.8 --g1 82.076975 --g2 159.369837 --omega1 127.304623 \
        --omega2 307.304623 --inc 37.935 --quad --oct --peri --outerperi --canon \
        --time 3e5 --dt 1e4

## Other options

As in the Paper I version (see `C++_Version/README.md` at tag `paperI-v1` and `./kozai --help`): upstream
options unchanged (angles in degrees); `--units code`, `--clight`, `--raw`, `--outer_frame`,
`--restricted`, `--no_rn4`, `--no_rn2p5`, `--no_rn1p5`, `--lr_rn1p5`, `--lr_rn2p5`,
`--x_rn3p5 kst|lr|lrraw`, `--brown`, `--audit`, and the shorthands `--ours` (model
"Paper I") and `--lr_full`. The audit table of `--audit` has one more row here,
"(1,3)+(1,4) canonical [II]". Changes with respect to upstream: unit system made variable;
collision radius fixed to 2Gm/c² (upstream: 2·1476·m/MSUN, i.e. metres mixed with cm);
out-of-bounds read in `get_L1()` fixed; LR's cross terms split by grade; (1,7/2) term and
Brown term added; outer-plane frame, restricted problem, raw output, audit; and `--canon`.

Models of Paper II, Table VII (all with `--quad --peri --outerperi`, plus `--oct`, `--rad`
where stated): "no 3BpN": none; "LR": `--cross_lim`; "Paper I": `--ours`; "this work":
`--canon`.

## Build and test

Requires a C++11 compiler and GSL (tested with g++ 13.3 and GSL 2.7.1 on Linux x86-64).

    make            # builds ./kozai (uses gsl-config if present, else -lgsl -lgslcblas -lm)
    make test       # tests/audit.sh, tests/x72sign.sh, tests/regress.sh, tests/canon.sh
    make clean      # removes kozai, kozai_orig and tests/upstream

Default flags `-O3 -ffp-contract=off` (no fused multiply-adds, so that the bit-for-bit
comparisons are portable); `CXX`, `CXXFLAGS` can be overridden.

| test | what | runtime | expected output |
|---|---|---|---|
| `tests/audit.sh` | Paper I Table VII with this executable (`--audit`) | < 1 s | ρ_ι, ρ_Ω = 1.46406, 2.59997 (KST) and 2.74027, 5.39947 (LR (A7), (A9), (A11)); `audit: PASS` |
| `tests/x72sign.sh` | sign of the (1,7/2) term: `--audit` at `--units code --m1 30 --m2 20 --m3 40 --a1 1 --a2 12 --e1 0.4 --e2 0.35 --g1 40 --g2 70 --omega1 0 --omega2 180 --inc 65 --clight 300 --quad --oct --x_rn3p5 kst`, ratio of the Ġ1 of the rows (1,7/2) and (0,4) | < 1 s | −3.4913e-3 (rel. tol. 1e-4) for m1, m2 = 30, 20 and 20, 30; +3.4913e-3 with `X72LEGACY=1` for 30, 20; `x72sign: PASS` |
| `tests/regress.sh` | (1) clones upstream (commit 1393251) into `tests/upstream`, checks `upstream.patch`, builds `kozai_orig`, compares 3 option sets without new options; (2) rebuilds the Paper I version by reversing `PaperI_to_PaperII.patch` and compares it with this executable for 5 option sets without `--canon` (incl. `--ours`, `--brown`, `--x_rn3p5 lrraw`, `--units code`, `--raw`) | 5 s + clone | upstream: stderr bit-identical, stdout identical apart from the radius lines; Paper I vs Paper II: stdout and stderr bit-identical; `regress: PASS` |
| `tests/canon.sh` | smoke test of `--canon` and of the `CANON*` switches (below) | 35 s (`FULL=1`: 6 min) | `canon: PASS` |

`tests/canon.sh` runs everything twice, with the corrected code and with `X72LEGACY=1`
(the submitted code), and checks each number against the value recorded when the test was
set up (tolerance `TOL`, default 1e-8) and against the value quoted in Paper II:

* **Sec. VII A, Fig. 5 (top), blue trajectory** — LR's librating system, initial
  e1 = 0.40, 2000 yr. Corrected code: full model e1 ∈ [0.2256392264, 0.4868826069], max
  l_z = −0.6453205444; c_Y Y only [0.1940174164, 0.5108847169], −0.6386087568; electric term
  only [0.3179513265, 0.4000583281], l_z constant to 3.1e-6. `X72LEGACY=1`:
  [0.2255926273, 0.4869301644], −0.6453083850; [0.1939742135, 0.5109001185], −0.6385957054;
  [0.3180174713, 0.4000006554] — identical to all digits to the research runs behind the
  paper. Both agree with the paper at the digits quoted ([0.226, 0.487], −0.6453;
  [0.194, 0.511], −0.6386; [0.318, 0.400], l_z constant to 1e-5).
* **Sec. VII B, Table VIII, orientation A** — LR's Fig. 8 system with octupole, no GW;
  e1,max over the first 1e4 yr (the maximum over 3e5 yr is reached there to 4 digits).
  Corrected code: no 3BpN 0.8103950850, Paper I 0.8862220047, this work 0.9091579712,
  electric only 0.8125641186, c_Z Z only 0.8123990757, c_Y Y only 0.9067380862,
  `CANON13=0` 0.9091608312, `CANON14Y=0` 0.8123990134. `X72LEGACY=1`: 0.8103950850,
  0.8863644097, 0.9092388688, 0.8115111369, 0.8114579693, 0.9068234241, 0.9092437052,
  0.8114570293, compared with Table VIII of the submitted paper (0.8104, 0.8864, 0.9092,
  0.8115, 0.8115, 0.9068). With `FULL=1` the runs cover 3e5 yr, as in Table VIII:
  corrected 0.8103972854, 0.8862220092, 0.9091603514, 0.8125709596, 0.8123999171,
  0.9067373599, 0.9091632121, 0.8123994368; with `X72LEGACY=1` the `--canon` rows
  (0.9092392631, 0.8115367549, 0.8114588310, 0.9068239385, 0.9092400949) are identical to
  the research runs, and no 3BpN (0.8103972854) and Paper I (0.8863644115) agree with the
  6-digit research table.

`tests/regress.sh` needs `git` and network access to github.com (or an existing clone in
`tests/upstream`). The output of `make test` from a clean checkout is in
`tests/expected_output.txt`.

## Files

| file | content |
|---|---|
| `kozai.C`, `kozai.h` | the code (Paper I version + `--canon`) |
| `vector.h` | 3-vector class (unchanged from upstream) |
| `upstream.patch` | `diff -ruN` upstream → this version |
| `PaperI_to_PaperII.patch` | `diff -ruN` Paper I version → this version |
| `Makefile` | `make`, `make test`, `make clean` |
| `tests/audit.sh`, `tests/x72sign.sh`, `tests/regress.sh`, `tests/canon.sh` | tests (relative paths only) |
| `tests/expected_output.txt` | log of `make test` |
