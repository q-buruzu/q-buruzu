#include "error_utils.h"
#include "hilbert.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

StateVector::StateVector(const std::vector<std::complex<double>>& otherVector) {
	vector = otherVector;
}

StateVector::StateVector(const StateVector& otherVector) {
	vector = otherVector.vector;
}

StateVector::StateVector(size_t elements) {
	vector = std::vector<std::complex<double>>(elements);
}

void StateVector::set(const std::vector<std::complex<double>>& otherVector) {
	vector = otherVector;
}

std::vector<std::complex<double>> StateVector::get() const {
	return vector;
}

std::complex<double>& StateVector::operator[](int i) {
	return vector[i];
}

const std::complex<double>& StateVector::operator[](int i) const {
	return vector[i];
}

size_t StateVector::size() const {
	return vector.size();
}

void StateVector::resize(size_t size) {
	vector.resize(size);
}

StateVector StateVector::conjugate() const {
	StateVector resultVector(vector.size());

	for (size_t i = 0; i < vector.size(); ++i) {
		resultVector[i] = std::conj(vector[i]);
	}

	return resultVector;
}

std::vector<double> StateVector::convert() const {
	std::vector<double> resultVector(vector.size());

	for (size_t i = 0; i < vector.size(); ++i) {
		resultVector[i] = vector[i].real();
	}

	return resultVector;
}

void StateVector::print() const {
	for (size_t i = 0; i < vector.size(); ++i) {
		std::cout << vector[i] << "\t";
	}

	std::cout << "\n";
}

StateVector HilbertSpace::add(const StateVector& vector, const StateVector& otherVector) const {
	throwError(sameDimensions(vector, otherVector), "VECTORS MUST BE SAME DIMENSIONS");

	StateVector resultVector(vector.size());

	for (size_t i = 0; i < vector.size(); ++i) {
		resultVector[i] = vector[i] + otherVector[i];
	}

	return resultVector;
}

StateVector HilbertSpace::scalarMultiply(std::complex<double> scalar, const StateVector& vector) const {
	StateVector resultVector(vector.size());

	for (size_t i = 0; i < vector.size(); ++i) {
		resultVector[i] = scalar * vector[i];
	}

	return resultVector;
}

std::complex<double> HilbertSpace::innerProduct(const StateVector& vector, const StateVector& otherVector) const {
	throwError(sameDimensions(vector, otherVector), "VECTORS MUST BE SAME DIMENSIONS");

	std::complex<double> innerProductValue = {0, 0};

	for (size_t i = 0; i < vector.size(); ++i) {
		innerProductValue += std::conj(vector[i]) * otherVector[i];
	}

	return innerProductValue;
}

double HilbertSpace::norm(const StateVector& vector) const {
	return std::sqrt(std::real(innerProduct(vector, vector)));
}

StateVector HilbertSpace::normalize(const StateVector& vector) const {
        std::vector<std::complex<double>> resultVector(vector.size());

        double factor = 1 / norm(vector);

        for (size_t i = 0; i < vector.size(); ++i) {
                resultVector[i] = vector[i] * factor;
        }

	return resultVector;
}

bool HilbertSpace::isCauchyConvergent(const std::vector<StateVector>& sequence, double epsilon = 1e-12) const {
	for (size_t i = 0; i < sequence.size(); ++i) {
		for (size_t j = i + 1; j < sequence.size(); ++j) {
			StateVector difference(sequence[i].size());

			for (size_t k = 0; k < difference.size(); ++k) {
				difference[k] = sequence[i][k] - sequence[j][k];
			}

			if (HilbertSpace::norm(difference) > epsilon) {
				return false;
			}
		}
	}

	return true;
}
