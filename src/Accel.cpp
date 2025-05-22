#include "..\include\Accel.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\global.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\IERS.hpp"
#include "..\include\timediff.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\Mjday_TDB.hpp"

Matrix Accel(double x, Matrix Y) {
	
    auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC]=
	IERS(AuxParam.Mjd_UTC + x/86400, 'l');
	
	auto[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = 
	timediff(UT1_UTC,TAI_UTC);
	
    double Mjd_UT1 = AuxParam.Mjd_UTC + x / 86400.0 + UT1_UTC / 86400.0;
    double Mjd_TT  = AuxParam.Mjd_UTC + x / 86400.0 + TT_UTC / 86400.0;
	
    Matrix P = PrecMatrix(SAT_Const::MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    double MJD_TDB = Mjday_TDB(Mjd_TT);
	auto [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun]= 
	JPL_Eph_DE430(MJD_TDB);
	
	r_Mercury=transpose(r_Mercury);
	r_Venus=transpose(r_Venus);
	r_Earth=transpose(r_Earth);
	r_Mars=transpose(r_Mars);
	r_Jupiter=transpose(r_Jupiter);
	r_Saturn=transpose(r_Saturn);
	r_Uranus=transpose(r_Uranus);
	r_Neptune=transpose(r_Neptune);
	r_Pluto=transpose(r_Pluto);
	r_Moon=transpose(r_Moon);
	r_Sun=transpose(r_Sun);
	
    Matrix r(3);
    for (int j = 1; j <= 3; ++j)
        r(1, j) = Y(j, 1);
    // Acceleration due to harmonic gravity field
    Matrix a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
	a=transpose(a);
    // Luni-solar perturbations
    if (AuxParam.sun)
        a = a + AccelPointMass(r, r_Sun, SAT_Const::GM_Sun);
    if (AuxParam.moon)
        a = a + AccelPointMass(r, r_Moon, SAT_Const::GM_Moon);
    // Planetary perturbations
    if (AuxParam.planets) {
        a = a + AccelPointMass(r, r_Mercury, SAT_Const::GM_Mercury);
        a = a + AccelPointMass(r, r_Venus, SAT_Const::GM_Venus);
        a = a + AccelPointMass(r, r_Mars, SAT_Const::GM_Mars);
        a = a + AccelPointMass(r, r_Jupiter, SAT_Const::GM_Jupiter);
        a = a + AccelPointMass(r, r_Saturn, SAT_Const::GM_Saturn);
        a = a + AccelPointMass(r, r_Uranus, SAT_Const::GM_Uranus);
        a = a + AccelPointMass(r, r_Neptune, SAT_Const::GM_Neptune);
        a = a + AccelPointMass(r, r_Pluto, SAT_Const::GM_Pluto);
    }
	a=transpose(a);
    Matrix dY(6, 1);
    for (int i = 1; i <= 3; ++i)
        dY(i, 1) = Y(i + 3, 1); 
    for (int i = 1; i <= 3; ++i)
        dY(i + 3, 1) = a(i, 1);
	return dY;
}

