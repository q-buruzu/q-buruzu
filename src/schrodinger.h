#ifndef SCHRODINGER_H
#define SCHRODINGER_H

#include "quantumstate.h"
#include "matrix.h"

Matrix findEvolution(double step);
void evolve(QuantumState state, Matrix evolution);

#endif
