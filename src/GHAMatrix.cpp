#include "..\include\GHAMatrix.hpp"
#include "..\include\R_z.hpp"
#include "..\include\gast.hpp"

 Matrix GHAMatrix(double Mjd_UT1){
	return R_z( gast(Mjd_UT1) );
 }
 