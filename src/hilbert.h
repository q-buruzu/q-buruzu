#ifndef HILBERT_H
#define HILBERT_H

#include <complex>
#include <vector>

class StateVector {
        public:
                StateVector(const std::vector<std::complex<double>>& initVector);
                StateVector(const StateVector& otherVector);
		StateVector(size_t elements);

                void set(const std::vector<std::complex<double>>& otherVector);
                std::vector<std::complex<double>> get() const;
		std::complex<double>& operator[](int i);
		const std::complex<double>& operator[](int i) const;
		size_t size() const;
		std::vector<double> convert() const;
                void print() const;

        private:
                std::vector<std::complex<double>> vector;
};

class HilbertSpace {
        public:
		friend class StateVector;

		StateVector add(const StateVector& vector, const StateVector& otherVector);
		StateVector scalarMultiply(std::complex<double> scalar, const StateVector& vector);
		std::complex<double> innerProduct(const StateVector& vector, const StateVector& otherVector);
		double norm(const StateVector& vector);
		bool isCauchyConvergent(const std::vector<StateVector>& vectorSequence, double epsilon);
};

#endif
