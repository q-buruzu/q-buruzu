#include "algorithms.h"
#include "error_utils.h"
#include "hilbert.h"
#include "matrix.h"
#include "quantumstate.h"
#include "random_generator.h"
#include "schrodinger.h"

#include <complex>
#include <cmath>
#include <iostream>

typedef std::complex<double> cd;

double hBar = 1.054571817e-34;
Matrix I({{{1, 0}, {0, 0}}, {{0, 0}, {1, 0}}});
Matrix X({{{0, 0}, {1, 0}}, {{1, 0}, {0, 0}}});
Matrix Y({{{0, 0}, {0, -1}}, {{0, 1}, {0, 0}}});
Matrix Z({{{1, 0}, {0, 0}}, {{0, 0}, {-1, 0}}});

std::vector<StateVector> extractQubitStates(StateVector state) {
        size_t nQubits = static_cast<size_t>(std::log2(state.size()));
        std::vector<StateVector> individualQubits;
        individualQubits.reserve(nQubits);

        StateVector copyState(state);
        HilbertSpace space;

        for (size_t i = 0; i < nQubits; ++i) {
                size_t qubitsLeft = nQubits - i;

                if (qubitsLeft == 1) {
                        StateVector lastQubit(space.normalize(copyState.get()));
                        individualQubits.push_back(lastQubit);
                        break;
                } else {
                        std::vector<Matrix> USV;
                        USV = algorithmSVD(copyState);
                        Matrix& U(USV[0]);
                        Matrix& sigma(USV[1]);
                        Matrix& V(USV[2]);

                        StateVector firstQubitState(2);
                        firstQubitState[0] = U[0][0];
                        firstQubitState[1] = U[1][0];
                        individualQubits.push_back(firstQubitState);

                        size_t width = V.getRows();
                        StateVector remainingQubitStates(width);

                        for (size_t j = 0; j < width; ++j) {
                                remainingQubitStates[j] = sigma[0][0] * V[j][0];
                        }

                        remainingQubitStates = space.normalize(remainingQubitStates);
                        copyState.set(remainingQubitStates.get());
                }
        }

        return individualQubits;
}

Matrix findHamiltonian(std::vector<StateVector> qubits) {
        std::vector<Matrix> pauliStrings;
        std::vector<size_t> chosenQubits;
        std::vector<Matrix> gates({I, X, Y, Z});
        size_t otherQubit;
        bool alreadyChosen = false;
	HilbertSpace space;

        for (size_t i = 0; i < qubits.size(); ++i) {
                for (size_t j = 0; j < chosenQubits.size(); ++j) {
                        if (i == chosenQubits[j]) {
                                alreadyChosen = true;
                                break;
                        }
                }

                if (!randChooseQubit() || alreadyChosen) {
                        alreadyChosen = false;
                        continue;
                }

                otherQubit = (i + randDeviationChoose() + qubits.size()) % qubits.size();

                pauliStrings.push_back(space.outerProduct((gates[randChooseGate()] * qubits[i]), (gates[randChooseGate()]) * qubits[otherQubit]));

                chosenQubits.push_back(i);
                chosenQubits.push_back(otherQubit);
        }

        Matrix resultMatrix(0, 0);

        if (pauliStrings.size() != 0) {
                for (size_t i = 0; i < pauliStrings.size(); ++i) {
                        resultMatrix = resultMatrix + pauliStrings[i];
                }
        }

        // handle nonchosen qubits

        return resultMatrix;
}

Matrix findEvolution(QuantumState system, double step) {
        // Matrix O(findHamiltonian(extractQubitStates(system.get())));
        Matrix O({{{1, 0}, {0, 1}}, {{0, 1}, {1, 0}}});

        O = (((O * cd{step, 0}) * cd{0, 1}) * cd{1 / hBar, 0}) * cd{-1, 0};
	O = O.identity() + O;

        return O;
}

void evolve(QuantumState state, Matrix evolution) {
        state.set(evolution * state.get());
}
