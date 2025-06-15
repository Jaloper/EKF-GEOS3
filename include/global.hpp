#ifndef _GLOBAL_
#define _GLOBAL_

#include "..\include\matrix.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <string.h>

typedef struct{
	double Mjd_UTC, Mjd_TT;
	int n, m, sun, moon, planets;
} Param;

extern Param AuxParam;
extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;
extern Matrix obs;
extern int n_eqn;

void eop19620101(int c);

void GGM03S(int n);

void DE430Coeff(int i, int j);

void GEOS3(int n);

#endif