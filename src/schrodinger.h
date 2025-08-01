#ifndef SCHRODINGER_H
#define SCHRODINGER_H

#include "hilbert.h"
#include "matrix.h"
#include "quantumstate.h"

std::vector<StateVector> extractQubitStates(StateVector state);
Matrix findEvolution(QuantumState system, double step);
void evolve(QuantumState state, Matrix evolution);

#endif
