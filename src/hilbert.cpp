#include "error_utils.h"
#include "hilbert.h"

#include <vector>
#include <iostream>
#include <stdexcept>
#include <complex>
#include <cmath>
#include <string>

StateVector::StateVector(std::vector<std::complex<double>> initVector) {
	vector = initVector;
}

StateVector::StateVector(const StateVector& otherVector) {
	vector = otherVector.vector;
}

StateVector::StateVector(int elements) {
	vector = std::vector<std::complex<double>>(elements);
}

void StateVector::set(std::vector<std::complex<double>> otherVector) {
	vector = otherVector;
}

std::vector<std::complex<double>> StateVector::get() {
	return vector;
}

std::complex<double>& StateVector::operator[](int i) {
	return vector[i];
}

size_t StateVector::size() const {
	return vector.size();
}

void StateVector::print() {
	for (size_t i = 0; i < vector.size(); ++i) {
		std::cout << vector[i] << "\t";
	}

	std::cout << "\n";
}

StateVector HilbertSpace::add(StateVector vector1, StateVector vector2) {
	throwError(sameDimensions(vector1, vector2), "VECTORS MUST BE SAME DIMENSIONS");

	StateVector resultVector(vector1.size());

	for (size_t i = 0; i < vector1.size(); ++i) {
		resultVector[i] = vector1[i] + vector2[i];
	}

	return resultVector;
}

StateVector HilbertSpace::scalarMultiply(std::complex<double> scalar, StateVector vector1) {
	StateVector resultVector(vector1.size());

	for (size_t i = 0; i < vector1.size(); ++i) {
		resultVector[i] = scalar * vector1[i];
	}

	return resultVector;
}

std::complex<double> HilbertSpace::innerProduct(StateVector vector1, StateVector vector2) {
	throwError(sameDimensions(vector1, vector2), "VECTORS MUST BE SAME DIMENSIONS");

	std::complex<double> innerProductValue = {0, 0};

	for (size_t i = 0; i < vector1.size(); ++i) {
		innerProductValue += std::conj(vector1[i]) * vector2[i];
	}

	return innerProductValue;
}

double HilbertSpace::norm(StateVector vector1) {
	return std::sqrt(std::real(innerProduct(vector1, vector1)));
}

bool HilbertSpace::isCauchyConvergent(std::vector<StateVector> vectorSequence, double epsilon = 1e-12) {
	for (size_t i = 0; i < vectorSequence.size(); ++i) {
		for (size_t j = i + 1; j < vectorSequence.size(); ++j) {
			StateVector difference(vectorSequence[i].size());

			for (size_t k = 0; k < difference.size(); ++k) {
				difference[k] = vectorSequence[i][k] - vectorSequence[j][k];

				if (HilbertSpace::norm(difference) > epsilon) return false;
			}
		}
	}

	return true;
}
