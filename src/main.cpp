#include "quantumstate.h"

#include <iostream>

int main() {
	size_t qubits;
	std::cout << "enter number of qubits blud\n";
	std::cin >> qubits;
	std::cout << "\n";

	QuantumState blud(qubits);
	blud.measure();
	return 0;
}
