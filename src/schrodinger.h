#ifndef SCHRODINGER_H
#define SCHRODINGER_H

#include "hilbert.h"
#include "matrix.h"
#include "quantumstate.h"

Matrix split(StateVector state);

void tempNormalize(StateVector vector);

double tempNorm(StateVector vector);

std::complex<double> tempInnerProduct(StateVector vector);

Matrix outerProduct(StateVector state);

Matrix identityPad(Matrix A, size_t size);

std::vector<Matrix> qrDecompose(Matrix A);

void SVD(StateVector state);

Matrix findEvolution(double step);

void evolve(QuantumState state, Matrix evolution);

#endif
