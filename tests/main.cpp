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
#include "..\include\gast.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\Accel.hpp"
#include "..\include\VarEqn.hpp"
#include "..\include\DEInteg.hpp"
#include <iostream>

using namespace std;

int main() {
    eop19620101(21413); //c=21413
	GGM03S(181);
	DE430Coeff(2285,1020);
    //cout << "eopdata\n" << eopdata << "\n";
 AuxParam.Mjd_UTC = 49746.1112847221;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;


    Matrix Y0_apr(6);

    Y0_apr(1) = 6221397.62857869;
	Y0_apr(2) = 2867713.77965738;
	Y0_apr(3) = 3006155.98509949;
    Y0_apr(4) = 4645.04725161807;
	Y0_apr(5) = -2752.21591588205;
	Y0_apr(6) = -7507.99940987033; 
	Y0_apr=transpose(Y0_apr);
	std::cout <<"Y0_apr\n"<< Y0_apr<<endl;

    
    Matrix Y = DEInteg(Accel,0,-134.999991953373,1e-13,1e-6,6,Y0_apr);

    std::cout <<"Y\n"<< Y<<endl;


    return 0;
}
