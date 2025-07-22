#ifndef SCHRODINGER_H
#define SCHRODINGER_H

#include "hilbert.h"
#include "matrix.h"
#include "quantumstate.h"

Matrix split(StateVector state);

Matrix outerProduct(StateVector state);

Matrix identityPad(Matrix A, size_t size);

std::vector<Matrix> qrDecompose(Matrix A);

bool offDiagonalValues(Matrix A, std::complex<double> tolerance);

Matrix algorithmQR(Matrix A);

std::vector<Matrix> algorithmSVD(StateVector state);

std::vector<StateVector> extractQubitStates(StateVector state);

Matrix findEvolution(double step);

void evolve(QuantumState state, Matrix evolution);

#endif
