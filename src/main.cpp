#include "matrix.h"
#include "hilbert.h"

int main() {
	HilbertSpace space;

	Matrix A({{{0, 0}, {0, 1}}, {{1, 0}, {1, 1}}});
	A = A * std::complex<double>{1, 2};
	A.print();

	StateVector b({{0, 1}, {0, 0}, {1, 0}, {8, 8}});
	b = space.scalarMultiply({0.25, 0.25}, b);
	b.print();
	return 0;
}
