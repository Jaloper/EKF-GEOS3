#include "..\include\matrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\Frac.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\Position.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\sign_.hpp"
#include "..\include\timediff.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\IERS.hpp"
#include "..\include\global.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\LTC.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\gmst.hpp"
#include <iostream>

using namespace std;

int main() {
    eop19620101(21413); //c=21413
	GGM03S(181);
	DE430Coeff(2285,1020);
    //cout << "eopdata\n" << eopdata << "\n";
	double Mjd_UTC = 37668;
	auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC]= IERS(Mjd_UTC, 'l');
	
	cout<<"x_pole "<<x_pole<<endl;
cout<<"y_pole "<<y_pole<<endl;
cout<<"UT1_UTC "<<UT1_UTC<<endl;
cout<<"LOD "<<LOD<<endl;
cout<<"dpsi "<<dpsi<<endl;
cout<<"deps "<<deps<<endl;
cout<<"dx_pole "<<dx_pole<<endl;
cout<<"dy_pole "<<dy_pole<<endl;
cout<<"TAI_UTC "<<TAI_UTC<<endl;



    return 0;
}
