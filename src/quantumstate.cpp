#include "error_utils.h"
#include "hilbert.h"
#include "matrix.h"
#include "quantumstate.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <random>
#include <string>
#include <vector>

QuantumState::QuantumState(size_t qubits)
	: state((size_t) pow(2, qubits)) {
		dimension = pow(2, qubits);
		set(randGauss(0.0, 1.0));
}

void QuantumState::set(StateVector otherState) {
	state.set(otherState.get());
}

void QuantumState::set(std::vector<std::complex<double>> otherState) {
	state.set(otherState);
}

std::vector<std::complex<double>> QuantumState::randGauss(double mean, double stddev) {
	std::random_device rd;
	std::mt19937 gen(rd());

	std::normal_distribution<double> distribution(mean, stddev);

	std::vector<std::complex<double>> randomVector(dimension);

	for (size_t i = 0; i < dimension; ++i) {
		randomVector[i] = {distribution(gen), distribution(gen)};
	}

	return randomVector;
}

int QuantumState::randChoose() {
	std::random_device rd;
	std::mt19937 gen(rd());

	std::vector<double> probabilities = state.convert();

	std::discrete_distribution<> distribution(probabilities.begin(), probabilities.end());

	return distribution(gen);
}

void QuantumState::normalize() {
	std::vector<std::complex<double>> normalizedState(dimension);

	double factor = 1 / space.norm(state);

	for (size_t i = 0; i < dimension; ++i) {
		normalizedState[i] = state[i] * factor;
	}

	state.set(normalizedState);
}

std::string QuantumState::printKet(int value, size_t length) {
	std::string ket = "|";

	for (size_t i = 0; i < length; ++i) {
		ket += ((value >> (length - i - 1)) & 1) ? '1' : '0';
	}

	return ket + ">";
}

void QuantumState::probAmplitudes() {
	normalize();

	for (size_t i = 0; i < dimension; ++i) {
		state[i] *= std::conj(state[i]);
	}
}

void QuantumState::measure() {
	probAmplitudes();
	int decimalValue = randChoose();

	std::string binaryValue = printKet(decimalValue, std::log2(dimension));
	std::cout << binaryValue << "\n";

	for (size_t i = 0; i < dimension; ++i) {
		state[i] = {0, 0};
	}

	state[decimalValue] = {1, 0};
}
