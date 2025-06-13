#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include "..\include\matrix.hpp"
#include <cmath>

Matrix TimeUpdate(Matrix P, Matrix Phi); // I can´t define a default Qdt because I don´t know its dimensions.
Matrix TimeUpdate(Matrix P, Matrix Phi, Matrix Qdt); 

#endif