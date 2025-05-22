#ifndef _TIMEDIFF_
#define _TIMEDIFF_

#include "..\include\matrix.hpp"
#include <cmath>

std::tuple<double,double,double,double,double> timediff(double UT1_UTC, double TAI_UTC);

#endif