#include "io_handler.h"
#include "main_utils.h"
#include "matrix.h"
#include "quantumstate.h"
#include "schrodinger.h"

#include <atomic>
#include <chrono>
#include <iostream>
#include <new>
#include <string>
#include <thread>

double simTime;
double timeStep;
double maxTime;

std::atomic<bool> running;

size_t qubits;
QuantumState* blud;
Matrix evolution(2, 2);

void initialize(){
        simTime = 0;

        std::cout << "Enter number of qubits:\n";
        std::cin >> qubits;
        std::cout << "\n";

        std::cout << "Enter time step length:\n";
        std::cin >> timeStep;
        std::cout << "\n";

        std::cout << "Enter simulation length:\n";
        std::cin >> maxTime;
        std::cout << "\n";

        blud = new QuantumState(qubits);
        evolution = findEvolution(*blud, timeStep);

	running = true;
}

double timeOperation() {
	using clock = std::chrono::steady_clock;
	auto start = clock::now();

	inputHandler();
	evolve(*blud, evolution);
	simTime += timeStep;

	auto end = clock::now();
	std::chrono::duration<double> elapsed = end - start;

	return (std::max((double) 0, (timeStep - elapsed.count())));
}

void mainLoop() {
	double delay = timeOperation();

        while (simTime < maxTime) {
		std::this_thread::sleep_for(std::chrono::duration<double>(delay));

		std::cout << simTime << "\n";

		inputHandler();

		while(!running) {
			inputHandler();
		}

		evolve(*blud, evolution);
                simTime += timeStep;
        }
}

void measure() {
        blud->measure();
}
