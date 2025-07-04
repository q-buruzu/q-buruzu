#include <vector>
#include <iostream>
#include <stdexcept>
#include <complex>
#include <string>

class HilbertSpace {
        public:
		HilbertSpace(StateVector vector) {
			stateVector = vector;
		}

		StateVector operator+(StateVector otherVector) {
			throwError(sameDimensions(otherVector), "VECTORS MUST BE SAME SIZE");

			StateVector resultVector(stateVector.size());

			for (int i = 0; i < stateVector.size(); ++i) {
				resultVector[i] = stateVector[i] + otherVector[i];
			}

			return resultVector;
		}

		StateVector operator*(std::complex<double> scalar) {
                        StateVector resultVector(stateVector.size());

                        for (int i = 0; i < stateVector.size(); ++i) {
                                resultVector[i] = scalar * stateVector[i];
                        }

                        return resultVector;
                }

		std::complex<double> innerProduct(StateVector otherVector) {
                        throwError(sameDimensions(otherVector), "VECTORS MUST BE SAME SIZE");

                        std::complex<double> innerProductValue;

                        for (int i = 0; i < stateVector.size(); ++i) {
                                innerProductValue += stateVector[i] * otherVector[i];
                        }

                        return innerProductValue;
                }

	private:
		StateVector stateVector;

		bool sameDimensions(StateVector otherVector) {
			return (stateVector.size() == otherVector.size());
		}

		void throwError(bool condition, std::string message) {
			if (!condition) {
				throw std::invalid_argument(message);
			}
		}
};

class StateVector {
        public:
                StateVector(std::vector<std::complex<double>> initVector) {
                        vector = initVector;
                }

                StateVector(StateVector otherVector) {
                        vector = otherVector.vector;
                }

		StateVector(int elements) {
			vector = std::vector<std::complex<double>>(elements);
		}

                void set(std::vector<std::complex<double>> otherVector) {
                        vector = otherVector;
                }

                std::vector<std::complex<double>> get() {
                        return vector;
                }

		std::complex<double>& operator[](int i) {
			return vector[i];
		}

		int size() {
			return vector.size();
		}

                void print() {
                        for (int i = 0; i < vector.size(); ++i) {
                                std::cout << vector[i] << "\t";
                        }

                        std::cout << "\n";
                }

        private:
                std::vector<std::complex<double>> vector;
};
