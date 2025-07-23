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

	for (auto &fart : data) {
		std::cout << "(" << fart.real() << (fart.imag() >= 0 ? "+" : "") << fart.imag() << "i) ";
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
/*
	StateVector state1({{1, 1}, {1, 0}});
	StateVector state2({{2, 1}, {1, 3}});
	StateVector state3({{1, 0}, {0, 2}});

	StateVector yurt(tensorProduct(tensorProduct(state1, state2), state3));

	printState(yurt);
	std::cout << "\n\n";

	std::vector<StateVector> bart;
	bart = extractQubitStates(yurt);

	printStateList(bart);


/*	Matrix A({{{0, 1}, {1, 1}}, {{2, 0}, {1, 3}}});
	Matrix B({{{0, 1}, {1, 1}}, {{2, 0}, {1, 3}}});

	Matrix C(kroneckerProduct(A, B));
	C.print();
*/

        initialize();

        std::thread listener(inputListener);
        listener.detach();

        mainLoop();

        delete blud;

        return 0;
}
