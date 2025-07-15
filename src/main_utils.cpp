#include "main_utils.h"
#include "matrix.h"
#include "quantumstate.h"
#include "schrodinger.h"

#include <iostream>
#include <new>

double simTime;
double timeStep;
double maxTime;

size_t qubits;
QuantumState* blud;
Matrix evolution(2, 2);

void initialize(){
        simTime = 0;

        std::cout << "Enter time step length:\n";
        std::cin >> timeStep;
        std::cout << "\n";

        std::cout << "Enter simulation length:\n";
        std::cin >> maxTime;
        std::cout << "\n";

        std::cout << "Enter number of qubits:\n";
        std::cin >> qubits;
        std::cout << "\n";

        blud = new QuantumState(qubits);
        evolution = findEvolution(timeStep);
}

void mainLoop() {
        while (simTime < maxTime) {
                evolve(*blud, evolution);
                simTime += timeStep;
        }
}

void measure() {
        blud->measure();
}
