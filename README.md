# Fork with the post-Newtonian cross terms of Arakida (Papers I and II)

This is a fork of Carl Rodriguez's Kozai code, <https://github.com/carlrodriguez/kozai>
(branch `crossterms`), which contains the three-body relativistic cross terms of
Halston Lim's version, <https://github.com/halstonblim/kozai> (branch `Sept2019_crossterms`;
Lim & Rodriguez 2020). The branch `arakida-crossterms` adds, in `C++_Version/`, the
modifications used in

* H. Arakida, *Post-Newtonian secular dynamics of hierarchical triples. I*, submitted to
  Phys. Rev. D, arXiv:2609.36990 (Paper I; tag `paperI-v1`);
* H. Arakida, *Post-Newtonian secular dynamics of hierarchical triples. II*, submitted to
  Phys. Rev. D (Paper II; the canonical cross terms, option `--canon`; tag `paperII-v1`).

See `C++_Version/README.md` for the changes, the new options and the tests
(`make && make test` in `C++_Version/`). The Fortran version is not modified.

**Citation.** If you use this code in a paper, please cite the papers requested by the
original authors below — Rodriguez & Antonini (2018), Antonini et al. (2017) and
Lim & Rodriguez (2020) — and, for the modifications, Papers I and II.

The README of the original repository follows unchanged.

---

# Kozai Codes with orbit-averaged Spin Terms and three-body relativistic effects

Two codes for integrating the secular equations of motions for the three-body problem, including the quadrupole and octupole terms (see Tremaine+2009, Petrovich2015 and Liu+2015).

Also includes relativistic terms for the inner binary, including:  
----1pN Pericenter precession (ala Ford+2000)  
----Gravitational Wave emission (ala Blaes+2002)  
----Spin-Orbit coupling (ala me, following Barker and O'Connell 1975)  
----Spin-Spin coupling (same)  
Also includes cross terms from relativistic three-body effects (Lim and Rodriguez 2020)
You are welcome to use either code in any way you like, though we ask that if you use any results from it in a paper that you cite Rodriguez and Antonini (2018, https://arxiv.org/abs/1805.08212) and Antonini et al (2017, https://arxiv.org/abs/1711.07142), and Lim and Rodriguez (2020, https://arxiv.org/abs/2001.03654)

