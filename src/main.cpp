#include "algorithms.h"
#include "hilbert.h"
#include "io_handler.h"
#include "main_utils.h"
#include "matrix.h"
#include "schrodinger.h"

#include <iostream>
#include <thread>

extern QuantumState* blud;

StateVector tensorProduct(StateVector state1, StateVector state2) {
        size_t dim1 = state1.size();
        size_t dim2 = state2.size();

        StateVector resultVector(dim1 * dim2);

        for (size_t i = 0; i < dim1; ++i) {
                for (size_t j = 0; j < dim2; ++j) {
                        resultVector[(i * dim2) + j] = state1[i] * state2[j];
                }
        }

        return resultVector;
}


void printState(const StateVector &state) {
        auto data = state.get();

        for (auto &thing : data) {
                std::cout << "(" << thing.real() << (thing.imag() >= 0 ? "+" : "") << thing.imag() << "i) ";
                std::cout << "\n";
        }

        std::cout << "\n";
}

void printStateList(const std::vector<StateVector> &list) {
        for (size_t i = 0; i < list.size(); ++i) {
                std::cout << "qubit " << i << ": ";
                printState(list[i]);
        }
}

int main() {
        StateVector state1({{1, 0}, {0, 0}});
	state1.print();

        StateVector state2({{0, 0}, {1, 0}});
	state2.print();

        StateVector yurt(tensorProduct(state1, state2));

        std::cout << "\n";
        printState(yurt);
        std::cout << "\n";

        std::vector<StateVector> bart;
        bart = extractQubitStates(yurt.get());

        printStateList(bart);


/*
	Matrix A({{{0, 0}, {0, 0}}, {{0, 0}, {1, 0}}});
	std::cout << "A:\n";
	A.print();

	std::vector<Matrix> QR = qrDecompose(A);

	std::cout << "Q:\n";
	QR[0].print();

	std::cout << "R:\n";
	QR[1].print();
*/

/*
        initialize();

        std::thread listener(inputListener);
        listener.detach();

        mainLoop();

        delete blud;
*/

        return 0;
}
