#include "matrix.h"
#include "hilbert.h"

#include <stdexcept>
#include <string>

void throwError(bool condition, std::string message) {
	if (!condition) {
		throw std::invalid_argument(message);
	}
}

bool sameDimensions(const StateVector& currentVector, const StateVector& otherVector) {
	return (currentVector.size() == otherVector.size());
}

bool sameDimensions(const Matrix& currentMatrix, const Matrix& otherMatrix) {
        return (currentMatrix.getColumns() == otherMatrix.getColumns() && currentMatrix.getRows() == otherMatrix.getRows());
}

bool multiplyApplicable(const Matrix& currentMatrix, const Matrix& otherMatrix) {
        return (currentMatrix.getColumns() == otherMatrix.getRows());
}

bool isSquare(const Matrix& currentMatrix) {
        return (currentMatrix.getRows() == currentMatrix.getColumns());
}

bool isZeroColumn(const Matrix& currentMatrix, int column) {
        const auto& copyMatrix = currentMatrix.get();

	for (size_t i = 0; i < currentMatrix.getRows(); ++i) {
                if (copyMatrix[i][column] != std::complex<double>{0.0, 0.0}) {
                        return false;
                }
        }

        return true;
}
