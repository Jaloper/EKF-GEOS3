#include "..\include\PoleMatrix.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_x.hpp"

Matrix PoleMatrix (double xp,double yp) {
	Matrix Ry = R_y(-xp);
	Matrix Rx = R_x(-yp);
	return Ry * Rx;
}