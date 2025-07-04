#ifndef ERROR_UTILS_H
#define ERROR_UTILS_H

#include "matrix.h"
#include "hilbert.h"

#include <string>

void throwError(bool condition, std::string message);
bool sameDimensions(const StateVector& currentVector, const StateVector& otherVector);
bool sameDimensions(const Matrix& currentMatrix, const Matrix& otherMatrix);
bool multiplyApplicable(const Matrix& currentMatrix, const Matrix& otherMatrix);
bool isSquare(const Matrix& currentMatrix);
bool isZeroColumn(const Matrix& currentMatrix, int column);

#endif
