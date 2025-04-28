#ifndef _IERS_
#define _IERS_

#include "..\include\matrix.hpp"
#include <cmath>

std::tuple<double, double, double, double, double, double, double, double, double> 
IERS(
    Matrix eop, double Mjd_UTC, char interp
);

#endif