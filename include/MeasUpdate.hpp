#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include "..\include\matrix.hpp"
#include <cmath>

std::tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix x, Matrix z, Matrix g, Matrix s, Matrix G, Matrix P, int n);

#endif