#include "matrix.h"
#include "quantumstate.h"
#include "schrodinger.h"

#include <complex>
#include <cmath>
#include <iostream>

typedef std::complex<double> cd;

double hBar = 1.054571817e-34;
Matrix I({{cd{1, 0}, cd{0, 0}}, {cd{0, 0}, cd{1, 0}}});

Matrix findEvolution(double step) {
	Matrix O({{cd{1, 0}, cd{0, 1}}, {cd{0, 1}, cd{1, 0}}});

	O = (((O * cd{step, 0}) * cd{0, 1}) * cd{1 / hBar, 0}) * cd{-1, 0};

	O = I + O;

	return O;
}

void evolve(QuantumState state, Matrix evolution) {
	state.set(evolution * state.get());
}
