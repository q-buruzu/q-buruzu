#ifndef SCHRODINGER_H
#define SCHRODINGER_H

#include "matrix.h"
#include "quantumstate.h"

Matrix findEvolution(double step);
void evolve(QuantumState state, Matrix evolution);

#endif
