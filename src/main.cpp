#include "hilbert.h"
#include <iostream>
#include "main_utils.h"
#include "matrix.h"
#include "schrodinger.h"

/* StateVector tensorProduct(StateVector state1, StateVector state2) {
	size_t dim1 = state1.size();
	size_t dim2 = state2.size();

	Matrix tensorProductResult;

	for (size_t i = 0; i < dim1; ++i) {
		
	}
}
*/

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
	QuantumState yurt(3);

	printState(yurt.get());

	std::vector<StateVector> bart;
	bart = extractQubitStates(yurt.get());

	printStateList(bart);


/*	Matrix A({{{0, 1}, {1, 1}}, {{2, 0}, {1, 3}}});
	Matrix B({{{0, 1}, {1, 1}}, {{2, 0}, {1, 3}}});

	Matrix C(kroneckerProduct(A, B));
	C.print();
*/
	return 0;
}
