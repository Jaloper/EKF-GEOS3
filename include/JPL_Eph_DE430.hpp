#ifndef _JPL_EPH_DE430_
#define _JPL_EPH_DE430_

#include "..\include\matrix.hpp"
#include <cmath>

std::tuple<Matrix&, Matrix&, Matrix&, Matrix&, Matrix&,
 Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&> JPL_Eph_DE430(double Mjd_TDB);

#endif