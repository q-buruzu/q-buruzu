#ifndef QUANTUMSTATE_H
#define QUANTUMSTATE_H

#include "error_utils.h"
#include "hilbert.h"
#include "matrix.h"

#include <complex>
#include <vector>

class QuantumState {
	public:
		QuantumState(size_t qubits);
		void set(StateVector otherState);
		void set(std::vector<std::complex<double>> otherState);
		std::vector<std::complex<double>> randGauss(double mean, double stddev) const;
		int randChoose() const;
		void normalize();
		std::string printKet(int value, size_t length) const;
		void probAmplitudes();
		void measure();

	private:
		HilbertSpace space;
		StateVector state;
		size_t dimension;
};

#endif
