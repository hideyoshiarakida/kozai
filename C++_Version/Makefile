# Secular three-body code (modified version of C. L. Rodriguez's code with the
# cross terms of Lim & Rodriguez 2020), Paper I version.
#
#   make          build ./kozai
#   make test     run the tests in tests/ (needs git and network access for
#                 the regression test against the upstream code)
#   make clean    remove the executables and the files made by the tests
#
# -ffp-contract=off keeps the arithmetic free of fused multiply-adds on every
# platform, so that the bit-for-bit comparisons in tests/ are meaningful.

CXX      ?= g++
CXXFLAGS ?= -O3 -ffp-contract=off
GSL_CFLAGS := $(shell gsl-config --cflags 2>/dev/null)
GSL_LIBS   := $(shell gsl-config --libs 2>/dev/null || echo -lgsl -lgslcblas -lm)

kozai: kozai.C kozai.h vector.h
	$(CXX) $(CXXFLAGS) $(GSL_CFLAGS) -o $@ kozai.C $(GSL_LIBS)

test: kozai
	@cd tests && ./audit.sh && ./x72sign.sh && ./regress.sh

clean:
	rm -f kozai kozai_orig
	rm -rf tests/upstream tests/work

.PHONY: test clean
