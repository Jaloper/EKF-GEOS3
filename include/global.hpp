#ifndef _GLOBAL_
#define _GLOBAL_

#include "..\include\matrix.hpp"
#include <cmath>

extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;

void eop19620101(int c);

void GGM03S(int n);

void DE430Coeff(int i, int j);

#endif