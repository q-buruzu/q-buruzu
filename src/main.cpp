#include "matrix.h"

int main() {
	Matrix A({{{0, 0}, {0, 1}}, {{1, 0}, {1, 1}}});
	A.print();

	Matrix B({{3,1}, {2,9}});
	B.set(A.invert());
	B.print();
	return 0;
}
