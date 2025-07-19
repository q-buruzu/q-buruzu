#include "main_utils.h"
#include "matrix.h"
#include "schrodinger.h"

int main() {
	Matrix A({{{1, 0}, {0, 1}, {0, 0}}, {{2, 1}, {0, 1}, {3, 0}}, {{1, 0}, {2, 2}, {3, 1}}});

	std::vector<Matrix> thing;

	thing = qrDecompose(A);

	Matrix B(thing[0] * thing[1]);
	A.print();
	thing[0].print();
	thing[1].print();
	B.print();

	return 0;
}
