#include "error_utils.h"
#include "hilbert.h"
#include "matrix.h"

#include <bitset>
#include <cmath>
#include <complex>
#include <random>
#include <string>
#include <vector>

size_t qubits = 1;

class QuantumState {
	public:
		QuantumState(size_t qubits) {
			dimension = pow(2, qubits);
			set(randGauss());
		}

		void set(StateVector otherState) {
			state.set(otherState.get());
		}

		void set(std::vector<std::complex<double>> otherState) {
			state.set(otherState);
		}

		std::vector<std::complex<double>> randGauss(double mean = 0, double stddev = 1) {
			std::random_device rd;
			std::mt19937 gen(rd());

			std::normal_distribution<double> distribution(mean, stddev);

			std::vector<std::complex<double>> randomVector(dimension);

			for (size_t i = 0; i < dimension; ++i) {
				randomVector[i] = {distribution(gen), distribution(gen)};
			}

			return randomVector;
		}

		int randChoose() {
			std::random_device rd;
			std::mt19937 gen(rd());

			std::discrete_distribution<> distribution(state.convert());

			return distribution(gen);
		}

		void normalize() {
			std::vector<std::complex<double>> normalizedState(dimension);

			double factor = 1 / space.norm(state);

			for (size_t i = 0; i < dimension; ++i) {
				normalizedState[i] = state[i] * factor;
			}

			state.set(normalizedState);
		}

		void probAmplitudes() {
			normalize();

			for (size_t i = 0; i < dimension; ++i) {
				state[i] *= std::conj(state[i]);
			}
		}

		void measure() {
			probAmplitudes();

			int decimalValue = randChoose()
			std::string binaryValue = std::bitset<(int)std::log2(dimension)>(decimalValue).to_string();

			std::cout << binaryValue << "\n";

			for (size_t i = 0; i < dimension; ++i) {
				state[i] = {0, 0};
			}

			state[decimalValue] = {1, 0};
		}

	private:
		HilbertSpace space;
		StateVector state;
		size_t dimension;
};
