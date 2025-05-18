#include "..\include\gast.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\gmst.hpp"

double gast (double Mjd_UT1){

return fmod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*M_PI );
}