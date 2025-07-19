#include "main_utils.h"
#include "matrix.h"
#include "schrodinger.h"

int main() {
	Matrix A({{{1, 0}, {0, 1}, {0, 0}}, {{2, 1}, {0, 1}, {3, 0}}, {{1, 0}, {2, 2}, {3, 1}}});

	A = identityPad(A, 5);
	A.print();

	/*
	std::vector<Matrix> thing;

	thing = qrDecompose(A);

	A.print();

	thing[0].print();
	thing[1].print();

	Matrix B(thing[0] * thing[1]);
	B.print();
	*/

	return 0;
}
