#include "main_utils.h"
#include "matrix.h"
#include "schrodinger.h"

int main() {
	Matrix A({{{1, 1}, {3, 1}, {0, 0}}, {{2, 1}, {0, 1}, {2, 0}}, {{1, 4}, {2, 2}, {1, 1}}});

	std::vector<Matrix> thing;
	thing = qrDecompose(A);

	Matrix Q(thing[0]);
	Matrix R(thing[1]);

	Matrix QR(thing[0] * thing[1]);
	QR.roundValues();

	A.print();
	Q.print();
	R.print();
	QR.print();

	return 0;
}
