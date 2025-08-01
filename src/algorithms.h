#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include "matrix.h"
#include "hilbert.h"

std::vector<Matrix> qrDecompose(Matrix A);
bool offDiagonalValues(Matrix A, std::complex<double> tolerance);
Matrix algorithmQR(Matrix A);
std::vector<Matrix> algorithmSVD(StateVector state);

#endif
