#ifndef HILBERT_H
#define HILBERT_H

#include <vector>
#include <complex>


class StateVector {
        public:
                StateVector(std::vector<std::complex<double>> initVector);
                StateVector(const StateVector& otherVector);
		StateVector(int elements);

                void set(std::vector<std::complex<double>> otherVector);
                std::vector<std::complex<double>> get();
		std::complex<double>& operator[](int i);
		size_t size() const;
                void print();

        private:
                std::vector<std::complex<double>> vector;
};

class HilbertSpace {
        public:
		friend class StateVector;

		StateVector add(StateVector vector1, StateVector vector2);
		StateVector scalarMultiply(std::complex<double> scalar, StateVector vector1);
		std::complex<double> innerProduct(StateVector vector1, StateVector vector2);
		double norm(StateVector vector1);
		bool isCauchyConvergent(std::vector<StateVector> vectorSequence, double epsilon);
};

#endif
