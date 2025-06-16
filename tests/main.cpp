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
#include <iomanip>
using namespace std;

int main() {

	DE430Coeff(2285,1020);
	GGM03S(181);
	eop19620101(21413);
	GEOS3(46);
	

	double sigma_range = 92.5;          // [m]
	double sigma_az = 0.0224*SAT_Const::Rad; // [rad]
	double sigma_el = 0.0139*SAT_Const::Rad; // [rad]

// Kaena Point station
	double lat = SAT_Const::Rad*21.5748;     // [rad]
	double lon = SAT_Const::Rad*(-158.2706); // [rad]
	double alt = 300.20;                // [m]

	Matrix Rs = Position(lon, lat, alt);
		Rs=transpose(Rs);
	double Mjd1 = obs(1,1);
	double Mjd2 = obs(9,1);
	double Mjd3 = obs(18,1);

	Matrix r2(1,3);
		r2(1) = 6221397.62857869; r2(2) = 2867713.77965738;	r2(3) = 3006155.98509949;

	Matrix v2(1,3);
		v2(1) = 4645.04725161807; v2(2) = -2752.21591588205; v2(3) = -7507.99940987033;
		
	Matrix Y0_apr = union_vector(r2,v2);
	double Mjd0 = Mjday(1995,1,29,02,38,0);

	double Mjd_UTC = obs(9,1);

	AuxParam.Mjd_UTC = Mjd_UTC;
	AuxParam.n      = 20;
	AuxParam.m      = 20;
	AuxParam.sun     = 1;
	AuxParam.moon    = 1;
	AuxParam.planets = 1;

	n_eqn  = 6;
	Matrix Y = DEInteg(Accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,transpose(Y0_apr));
	Matrix P = zeros(6,6);
	for (int i=1;i<=3;i++) P(i, i) = 1e8;
    for (int i=4;i<=6;i++) P(i, i) = 1e3;

	Matrix LT = LTC(lon,lat);
	Matrix yPhi = zeros(1,42);
	Matrix Phi  = zeros(6,6);
	// Measurement loop
	double t = 0.0;
	double t_old,Mjd_TT,Mjd_UT1;
	for (int i=1;i<=46;i++){    
		// Previous step
		t_old = t;
		Matrix Y_old = transpose(Y);
		// Time increment and propagation
		Mjd_UTC = obs(i,1);                       // Modified Julian Date
		t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]
		auto[x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(Mjd_UTC,'l');
		auto[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
		Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
		Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
		AuxParam.Mjd_UTC = Mjd_UTC;
		AuxParam.Mjd_TT = Mjd_TT;
		for (int ii=1;ii<=6;ii++){
			yPhi(ii) = Y_old(ii);
			for (int j=1;j<=6;j++){  
				if (ii==j) 
					yPhi(6*j+ii) = 1; 
				else
					yPhi(6*j+ii) = 0;
				
			}
		}   
		yPhi=transpose(yPhi);
					Matrix yAux = DEInteg (VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);
		
		yPhi=transpose(yAux);
		// Extract state transition matrices
		for (int j=1;j<=6;j++){  
			for (int i=1;i<=6;i++){
			Phi(i, j) = yPhi(6 * j + i);
			}
		}
		Y_old=transpose(Y_old);
		Y = DEInteg (Accel,0,t-t_old,1e-13,1e-6,6,Y_old);
		 Y=transpose(Y);
		// Topocentric coordinates
		double theta = gmst(Mjd_UT1);                    // Earth rotation
		Matrix U = R_z(theta);
		Matrix r = extract_vector(Y,1,3);
			r=transpose(r);
		Matrix s = LT*(U*r-Rs);                          // Topocentric position [m]
		// Time update
		Matrix Px = TimeUpdate(P, Phi);
		P=transpose(Px);
			s=transpose(s);
		// Azimuth and partials
		auto [Azim, Elev, dAds, dEds] = AzElPa(s);     // Azimuth, Elevation
		Matrix dAdY = union_vector(dAds*LT*U,zeros(3));
		// Measurement update
		Y=transpose(Y);
		Matrix obsAux(1);
		obsAux(1,1)=obs(i,2);
		Matrix AzimAux(1);
		AzimAux(1,1)=Azim;
		Matrix sigma_azAux(1);
		sigma_azAux(1)=sigma_az;
		auto[K, Yaux, Paux] = MeasUpdate ( Y, obsAux, AzimAux, sigma_azAux, dAdY, P, 6 );
		Y=transpose(Yaux);
		P=Paux;
		// Elevation and partials
		r = transpose(extract_vector(Y,1,3));
		s = LT*(U*r-Rs);                          // Topocentric position [m]
		s=transpose(s);
		auto [Azimaux, Elevaux, dAdsaux, dEdsaux] = AzElPa(s);     // Azimuth, Elevation
		
		Matrix dEdY = union_vector(dEdsaux*LT*U,zeros(3));
		// Measurement update
		Y=transpose(Y);
		obsAux(1,1)=obs(i,3);
		Matrix ElevauxAux(1);
		ElevauxAux(1,1)=Elevaux;
		Matrix sigma_elAux(1);
		sigma_elAux(1)=sigma_el;
		auto [K2, Yaux2, Paux2] = MeasUpdate ( Y, obsAux, ElevauxAux, sigma_elAux, dEdY, P, 6 );
		Y=transpose(Yaux2);
		P=Paux2;
		// Range and partials
		r = transpose(extract_vector(Y,1,3));
		s = LT*(U*r-Rs);                          // Topocentric position [m]
		s=transpose(s);
		double Dist = norm(s); Matrix dDds = s/Dist;         // Range
		Matrix dDdY = union_vector(dDds*LT*U,zeros(3));
		
		// Measurement update
		Y=transpose(Y);
		obsAux(1,1)=obs(i,4);
		Matrix DistAux(1);
		DistAux(1,1)=Dist;
		Matrix sigma_rangeAux(1);
		sigma_rangeAux(1)=sigma_range;
		
		auto [K3, Yaux3, Paux3] = MeasUpdate ( Y, obsAux, DistAux, sigma_rangeAux, dDdY, P, 6 );
		Y=Yaux3;
		P=Paux3;
	}	
	auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(obs(46,1),'l');
	auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400;
	AuxParam.Mjd_UTC = Mjd_UTC;
	AuxParam.Mjd_TT = Mjd_TT;

	Matrix Y0 = DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);

	Matrix Y_true(6);
		Y_true(1) = 5753.173e3; 
		Y_true(2) = 2673.361e3;
		Y_true(3) = 3440.304e3;
		Y_true(4) = 4.324207e3;
		Y_true(5) = -1.924299e3;
		Y_true(6) = -5.728216e3;
	Y0=transpose(Y0);
	printf("\nError of Position Estimation\n");
	printf("dX%10.1f [m]\n",Y0(1)-Y_true(1));
	printf("dY%10.1f [m]\n",Y0(2)-Y_true(2));
	printf("dZ%10.1f [m]\n",Y0(3)-Y_true(3));
	printf("\nError of Velocity Estimation\n");
	printf("dVx%8.1f [m/s]\n",Y0(4)-Y_true(4));
	printf("dVy%8.1f [m/s]\n",Y0(5)-Y_true(5));
	printf("dVz%8.1f [m/s]\n",Y0(6)-Y_true(6));


}
