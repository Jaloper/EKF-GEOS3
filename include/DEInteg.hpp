#ifndef _DEINTEG_
#define _DEINTEG_

#include "..\include\matrix.hpp"
#include <cmath>

Matrix DEInteg(Matrix func(double t, Matrix y), double t,double tout,double relerr,double abserr,int n_eqn,Matrix& y);

#endif