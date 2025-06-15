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
#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)
	#define EPSILON 1e-4
#define ASSERT_CLOSE(a, b) \
    _assert(fabs((a) - (b)) < EPSILON * fmax(fabs(a), fabs(b)))
int m_equals(Matrix A, Matrix B, double tolerance) {
    if (A.n_row != B.n_row || A.n_column != B.n_column)
        return 0;

    for (int i = 1; i <= A.n_row; i++) {
        for (int j = 1; j <= A.n_column; j++) {
            double a = A(i, j);
            double b = B(i, j);
            double abs_diff = fabs(a - b);
            double rel_diff = abs_diff / fmax(fabs(a), fabs(b));

            if (abs_diff > tolerance && rel_diff > tolerance) {
                printf("Differences at (%d, %d): A = %2.20lf, B = %2.20lf\n", i, j, a, b);
                return 0;
            }
        }
    }

    return 1;
}


int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_mult_01() {
	
    Matrix A(3, 3);
    A(1,1) = 0.9; A(1,2) = 0.1; A(1,3) = 0;
    A(2,1) = 0; A(2,2) = 0.95; A(2,3) = 0.05;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0.99;
    
    Matrix B(3, 3);
    B(1,1) = 0.1;  B(1,2) = 0.01; B(1,3) = 0.005;
    B(2,1) = 0.01;  B(2,2) = 0.2; B(2,3) = 0.01;
    B(3,1) = 0.005; B(3,2) = 0.01; B(3,3) = 0.15;
    
    Matrix C(3, 3);
    C(1,1) = 0.0910;  C(1,2) = 0.0290; C(1,3) = 0.0055;
    C(2,1) = 0.0097; C(2,2) = 0.1905; C(2,3) = 0.0170;
	C(3,1) = 0.0050; C(3,2) = 0.0099; C(3,3) = 0.1485;
    
    Matrix R = A * B;
	
    _assert(m_equals(C, R, 1e-4));
    return 0;
}

int m_div_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;
    
    Matrix B(2, 2);
    B(1,1) = 5; B(1,2) = 6;
    B(2,1) = 7; B(2,2) = 8;
    
    Matrix C(2, 2);
    C(1,1) = 3.0;  C(1,2) = -2.0;
    C(2,1) = 2.0;  C(2,2) = -1.0;
    
    Matrix R = A / B;
    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_assign_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;
    
    Matrix B=A;
    
    _assert(m_equals(A, B, 1e-8));
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_inv_01() {
    Matrix A(2, 2);
    A(1,1) = 4; A(1,2) = 7;
    A(2,1) = 2; A(2,2) = 6;
    
    Matrix B(2, 2);
    B(1,1) = 0.6;  B(1,2) = -0.7;
    B(2,1) = -0.2; B(2,2) = 0.4;
    
    Matrix R = inv(A);
    
    _assert(m_equals(B, R, 1e-10));
    return 0;
}

int m_eye_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 0;
    A(2,1) = 0; A(2,2) = 1;
    
    Matrix B=eye(2);
    
    _assert(m_equals(B, A, 1e-10));
    return 0;
}

int m_transpose_01() {
    Matrix A(2, 2);
    A(1,1) = 4; A(1,2) = 7;
    A(2,1) = 2; A(2,2) = 6;
    
    Matrix B(2, 2);
    B(1,1) = 4; B(1,2) = 2;
    B(2,1) = 7; B(2,2) = 6;
    
	Matrix R = transpose(A);
	
    _assert(m_equals(B, R, 1e-10));
    return 0;
}

int op_sum_01() {
    double d = 5.2;
	
	Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 0;
    A(2,1) = 3; A(2,2) = 4;
	
	Matrix B(2, 2);
    B(1,1) = 6.2; B(1,2) = 5.2;
    B(2,1) = 8.2; B(2,2) = 9.2;
	
	
	Matrix R = A + d;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int op_sub_01() {
    double d = -5.2;
	
	Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 0;
    A(2,1) = 3; A(2,2) = 4;
	
	Matrix B(2, 2);
    B(1,1) = 6.2; B(1,2) = 5.2;
    B(2,1) = 8.2; B(2,2) = 9.2;
	
	
	Matrix R = A - d;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int op_mult_01() {
    double d = 2;
	
	Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 0;
    A(2,1) = 3; A(2,2) = 4;
	
	Matrix B(2, 2);
    B(1,1) = 2; B(1,2) = 0;
    B(2,1) = 6; B(2,2) = 8;
	
	
	Matrix R = A * d;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int op_div_01() {
    double d = 2;
	
	Matrix A(2,2);
    A(1,1) = 1; A(1,2) = 0;
    A(2,1) = 3; A(2,2) = 4;
	
	Matrix B(2, 2);
    B(1,1) = 0.5; B(1,2) = 0;
    B(2,1) = 1.5; B(2,2) = 2;
	
	
	Matrix R = A / d;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_norm_01() {
	
	Matrix A(2,2);
    A(1,1) = 3; A(1,2) = 0;
    A(2,1) = 0; A(2,2) = 4;
	
	double result=5;
	double R = norm(A);
    ASSERT_CLOSE(result, R);
    
    return 0;
}

int m_dot_01() {
	
	Matrix A(1,2);
    A(1,1) = 4; A(1,2) = 7;
	
	double result=65;
	double R = dot(A,A);
    _assert(m_equals(result, R, 1e-10));
    
    return 0;
}

int m_cross_01() {
	
	Matrix A(1,3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	
	Matrix B(1,3);
	B(1,1) = 1; B(1,2) = 5; B(1,3) = 7;
	
	Matrix C(1,3);
	C(1,1) = -1; C(1,2) = -4; C(1,3) = 3;
	
    _assert(m_equals(cross(A,B), C, 1e-10));
    
    return 0;
}

int extract_vector_01() {
	
	Matrix A(1,3);
	A(1,1) = 1; A(1,2) = 5; A(1,3) = 7;
	
	Matrix B =  extract_vector(A,2,3); 
	
	Matrix C(1,2);
	C(1,1) = 5; C(1,2) = 7;
	
    _assert(m_equals(B, C, 1e-10));
    
    return 0;
}

int union_vector_01() {
    Matrix A(1,2);
    A(1,1) = 4; A(1,2) = 7;
    
    Matrix B(1,3);
    B(1,1) = 1; B(1,2) = 5; B(1,3) = 7;
    
    Matrix C(1,5);
    C(1,1) = 1; C(1,2) = 5; C(1,3) = 7; C(1,4) = 4; C(1,5) = 7;
    
    Matrix R = union_vector(B,A);
    
    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}

int m_extract_row_01() {
	
	Matrix A(2,3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	A(2,1) = 4; A(2,2) = 5;
	
	Matrix B(1,3);
	B(1,1) = 1; B(1,2) = 2; B(1,3) = 3;
	
	Matrix C = extract_row(A,1);
	
    _assert(m_equals(B, C, 1e-10));
    
    return 0;
}

int m_extract_col_01() {
	
	Matrix A(2,3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	A(2,1) = 4; A(2,2) = 5;
	
	Matrix B(2,1);
	B(1,1) = 2; B(2,1) = 5;
	
	Matrix C = extract_column(A,2);
	
    _assert(m_equals(B, C, 1e-10));
    
    return 0;
}

int m_assign_row_01() {
	
	Matrix A(2,3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	A(2,1) = 4; A(2,2) = 5;
	
	Matrix B(1,2);
	B(1,1) = 2; B(1,2) = 5;
	
	Matrix C = assign_row(A,B,2);
	
	Matrix R(2,3);
	R(1,1) = 1; R(1,2) = 2; R(1,3) = 3;
	R(2,1) = 2; R(2,2) = 5;

    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}

int m_assign_col_01() {
	
	Matrix A = zeros(2, 2);
	A(1,1) = 4; A(1,2) = 7;
	A(2,1) = 2.5; A(2,2) = 6;
	
	Matrix B= zeros(2,1);
	B(1,1) = 1.2; 
	B(2,1) = 2; 
	
	Matrix C =  assign_column(A,B,2);
	
	Matrix R= zeros(2,2);
	R(1,1) = 4; R(1,2) = 1.2;
	R(2,1) = 2.5; R(2,2) = 2;

    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}

int R_x_01() {
	Matrix A(3,3);
    A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = -0.44807361612917012694; A(2,3) = 0.89399666360055785042;
	A(3,1) = 0; A(3,2) = -0.89399666360055785042; A(3,3) = -0.44807361612917012694;
    
    Matrix B = R_x(90);
    
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int R_y_01() {
	Matrix A(3,3);
    A(1,1) = -0.44807361612917012694; A(1,2) = 0; A(1,3) = -0.89399666360055785042;
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 0.89399666360055785042; A(3,2) = 0; A(3,3) = -0.44807361612917012694;
    
    Matrix B = R_y(90);
    
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int R_z_01() {
	Matrix A(3,3);
    A(1,1) = -0.44807361612917012694; A(1,2) = 0.89399666360055785042; A(1,3) = 0;
	A(2,1) = -0.89399666360055785042; A(2,2) = -0.44807361612917012694; A(2,3) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;
    
    Matrix B = R_z(90);
    
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int AccelPointMass_01() {
	Matrix A(1,2);
	A(1,1) = 2; A(1,2) = 2;
	
	Matrix B(1,2);
	B(1,1) = 2; B(1,2) = 1;
    
	Matrix C(1,2);
	C(1,1) = -0.572433402239946; C(1,2) = -3.48621670111997;
    
	Matrix R = AccelPointMass(A,B,3.2);
	
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int Cheb3D_01() {
	Matrix A(5);
	A(1) = 1; A(2) = -0.5; A(3) = 0.25; A(4) = -0.1; A(5) = 0.05;
	
	Matrix B(5);
	B(1) = 0.5; B(2) = 0.3; B(3) = -0.15; B(4) = 0.1; B(5) = -0.05;
    
	Matrix C(5);
	C(1) = 0.2; C(2) = -0.1; C(3) = 0.1; C(4) = -0.05; C(5) = 0.02; 
    
	Matrix R(3);
	R(1) = 0.8; R(2) = 0.6; R(3) = 0.12;
	
    _assert(m_equals(Cheb3D(0.5, 5, 0, 1, A, B, C), R, 1e-10));
    
    return 0;
}

int EccAnom_01() {
	double r=  1.6946;
    
	double E = EccAnom(1.0, 0.7);
	
    ASSERT_CLOSE(E, r);
    
    return 0;
}

int Frac_01() {
	double r=  0.2;
    
	double E = Frac(1.2);
	
    ASSERT_CLOSE(E, r);
    
    return 0;
}

int MeanObliquity_01() {
	double r=  0.4091;
    
	double E = MeanObliquity(51544.5);
	
    ASSERT_CLOSE(E, r);
    
    return 0;
}

int Mjday_01() {
	double r=  6.0783e+04;
    
	int yr = 2025, mon = 4, day = 17, hr = 12, min = 30;
    double sec = 0.0;

    double E = Mjday(yr, mon, day, hr, min, sec);
	
    ASSERT_CLOSE(E, r);
    
    return 0;
}

int Mjday_TDB_01() {
	double r=  5.1544e+04;
    double E = Mjday_TDB(51544.5);
	
    ASSERT_CLOSE(E, r);
    
    return 0;
}

int Position_01() {
	double lon = 0.1745;
    double lat = 0.7854;
    double h = 100.0;

    Matrix A = Position(lon, lat, h);
	Matrix B(3);
    B(1,1) = 4.4490; B(1,2) = 0.7843; B(1,3) = 4.4874;
	Matrix R=B*1.0e+06;
	
    _assert(m_equals(A, R, 1e-4));
    
    return 0;
}

int SAT_Const_01() {
	double R_Earth= 6378137;
	
    ASSERT_CLOSE(R_Earth,SAT_Const::R_Earth);
    
    return 0;
}

int sign__01() {
	double R= -8.2;
	
    ASSERT_CLOSE(sign_(8.2,-2.7),R);
    
    return 0;
}

int timediff_01() {
	  double UT1_UTC = 0.3341;
    double TAI_UTC = 37.0;

	auto[UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC, TAI_UTC);

	double R_UT1_TAI=-36.665900;
	double R_UTC_GPS=-18;
	double R_UT1_GPS=-17.665900;
	double R_TT_UTC=69.184000;
	double R_GPS_UTC=18;
	
    ASSERT_CLOSE(UT1_TAI,R_UT1_TAI);
	ASSERT_CLOSE(UTC_GPS,R_UTC_GPS);
	ASSERT_CLOSE(UT1_GPS,R_UT1_GPS);
	ASSERT_CLOSE(TT_UTC,R_TT_UTC);
	ASSERT_CLOSE(GPS_UTC,R_GPS_UTC);
    
    return 0;
} 

int AzElPa_01() {
    Matrix s(3);
    s(1) = 1000.0;
    s(2) = 2000.0;
    s(3) = 3000.0;

    auto [Az,El,dAds,dEds] = AzElPa(s);

    double R_Az = 0.463648;
    double R_El = 0.930274;
    
    Matrix r_dAds(1, 3);
    r_dAds(1, 1) = 0.4; r_dAds(1, 2) = -0.2; r_dAds(1, 3) = 0;
    Matrix R_dAds = r_dAds * 1.0e-03;

    Matrix r_dEds(1, 3);
    r_dEds(1, 1) = -0.0958; r_dEds(1, 2) = -0.1917; r_dEds(1, 3) = 0.1597;
    Matrix R_dEds = r_dEds * 1.0e-03;

    ASSERT_CLOSE(Az, R_Az);
    ASSERT_CLOSE(El, R_El);
    _assert(m_equals(dAds, R_dAds, 1e-4));
    _assert(m_equals(dEds, R_dEds, 1e-4));

    return 0;
}


int IERS_01() {
	
	double Mjd_UTC = 37668;
	auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC]= IERS(Mjd_UTC, 'l');

	double R_x_pole=-1.0665e-07;
	double R_y_pole=1.0487e-06;
	double R_UT1_UTC=0.031144;
	double R_LOD=0.001496;
	double R_dpsi=3.0976e-07;
	double R_deps=3.1872e-08;
	double R_TAI_UTC=2;
	
	ASSERT_CLOSE(x_pole,R_x_pole);
	ASSERT_CLOSE(y_pole,R_y_pole);
	ASSERT_CLOSE(UT1_UTC,R_UT1_UTC);
	ASSERT_CLOSE(LOD,R_LOD);
	ASSERT_CLOSE(dpsi,R_dpsi);
	ASSERT_CLOSE(deps,R_deps);
	ASSERT_CLOSE(TAI_UTC,R_TAI_UTC);
    
    return 0;
} 

int Legendre_01() {
    int n = 2, m = 3;
	double fi = 1.0;
	auto [pnm, dpnm] = Legendre(n, m, fi);
	
	Matrix R_pnm(3,4);
    R_pnm(1,1) = 1; R_pnm(1,2) = 0; R_pnm(1,3) = 0; R_pnm(1,4) = 0;
	R_pnm(2,1) = 1.4574704987823 ; R_pnm(2,2) = 0.935831045210238; R_pnm(2,3) = 0; R_pnm(2,4) = 0;
	R_pnm(3,1) = 1.25691645573063; R_pnm(3,2) = 1.76084689542256; R_pnm(3,3) = 0.565313394670859; R_pnm(3,4) = 0;
	
	Matrix R_dpnm(3,4);
    R_dpnm(1,1) = 0; R_dpnm(1,2) = 0; R_dpnm(1,3) = 0; R_dpnm(1,4) = 0;
	R_dpnm(2,1) = 0.935831045210238; R_dpnm(2,2) = -1.4574704987823; R_dpnm(2,3) = 0; R_dpnm(2,4) = 0;
	R_dpnm(3,1) = 3.0498762872218; R_dpnm(3,2) = -1.61172976752398; R_dpnm(3,3) = -1.76084689542256; R_dpnm(3,4) = 0;
	
	_assert(m_equals(dpnm, R_dpnm, 1e-10));
	_assert(m_equals(pnm, R_pnm, 1e-10));
	
    return 0;
}

int NutAngles_01() {
  double R_dpsi = -5.1235e-05;
  double R_deps = 3.1716e-05;
  double Mjd_TT = 59945.0;
  auto [dpsi, deps] = NutAngles(Mjd_TT);
	ASSERT_CLOSE(dpsi,R_dpsi);
	ASSERT_CLOSE(deps,R_deps);
    return 0;
}

int TimeUpdate_01() {
   Matrix Phi(3, 3);
	Phi(1,1) = 0.9;  Phi(1,2) = 0.1;  Phi(1,3) = 0.0;
	Phi(2,1) = 0.0;  Phi(2,2) = 0.95; Phi(2,3) = 0.05;
	Phi(3,1) = 0.0;  Phi(3,2) = 0.0;  Phi(3,3) = 0.99;

	Matrix P(3, 3);
	P(1,1) = 0.1;    P(1,2) = 0.01;   P(1,3) = 0.005;
	P(2,1) = 0.01;   P(2,2) = 0.2;    P(2,3) = 0.01;
	P(3,1) = 0.005;  P(3,2) = 0.01;   P(3,3) = 0.15;
	
	Matrix Qdt(3, 3);
	Qdt(1,1) = 0.01; Qdt(1,2) = 0.0;  Qdt(1,3) = 0.0;
	Qdt(2,1) = 0.0;  Qdt(2,2) = 0.02; Qdt(2,3) = 0.0;
	Qdt(3,1) = 0.0;  Qdt(3,2) = 0.0;  Qdt(3,3) = 0.01;
	
    Matrix R(3, 3);
    R(1, 1) = 0.0948; R(1, 2) = 0.0278; R(1, 3) = 0.0054;
    R(2, 1) = 0.0278;  R(2, 2) = 0.2018; R(2, 3) = 0.0168;
    R(3, 1) = 0.0054;  R(3, 2) = 0.0168; R(3, 3) = 0.1570;
	
	Matrix P_updated = TimeUpdate(P, Phi, Qdt);
	_assert(m_equals(P_updated, R, 1e-4));
    return 0;
}

int AccelHarmonic_01() {
	Matrix E = eye(3);
	int n_max = 10;
	int m_max = 10;
	   Matrix r(3);
	r(1,1) = 7000e3;  r(1,2) = 1000e3;  r(1,3) = 1300e3;
	
	Matrix A= AccelHarmonic(r, E, n_max, m_max);
	
	Matrix R(3,1);
	R(1,1) = -7.5161;  
	R(2,1) = -1.0738;  
	R(3,1) = -1.3994;
	
	_assert(m_equals(A, R, 1e-4));
	return 0;
}

int EqnEquinox_01() {
	double Mjd_1 = 51544.5;
	double R = -6.1932e-05;
	
	ASSERT_CLOSE(R,EqnEquinox(Mjd_1));
	return 0;
}

int JPL_Eph_DE430_01() {
	auto [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun]= JPL_Eph_DE430(51544.5);
	Matrix R_Mercury(3, 1);
	R_Mercury(1, 1) =  7037307364.07061; 
    R_Mercury(2, 1) = -192685385065.559;
    R_Mercury(3, 1) = -87549491129.0525;
	
	Matrix R_Venus(3, 1);
	R_Venus(1, 1) =  -80957460292.0907; 
    R_Venus(2, 1) = -139679946184.946;
    R_Venus(3, 1) = -53870531361.9661;
	
	Matrix R_Earth(3, 1);
	R_Earth(1, 1) =  -27566632259.0857; 
    R_Earth(2, 1) = 132361428301.196;
    R_Earth(3, 1) = 57418647273.1691;
	
	Matrix R_Mars(3, 1);
	R_Mars(1, 1) =  234547174398.74; 
    R_Mars(2, 1) = -132547798063.593;
    R_Mars(3, 1) = -63085880695.2492;
	
	Matrix R_Jupiter(3, 1);
	R_Jupiter(1, 1) =  625066612307.678; 
    R_Jupiter(2, 1) = 276628952849.716;
    R_Jupiter(3, 1) = 103337598704.726;
	
	Matrix R_Saturn(3, 1);
	R_Saturn(1, 1) =  984884157624.263; 
    R_Saturn(2, 1) = 790958241347.093;
    R_Saturn(3, 1) = 282744146884.643;
	
	Matrix R_Uranus(3, 1);
	R_Uranus(1, 1) =  2185474057263.34; 
    R_Uranus(2, 1) = -2003668343398.91;
    R_Uranus(3, 1) =  -907525371048.23;
	
	Matrix R_Neptune(3, 1);
	R_Neptune(1, 1) =  2541545321752.3; 
    R_Neptune(2, 1) = -3570531546142.19;
    R_Neptune(3, 1) = -1527270114241.62;
	
	Matrix R_Pluto(3, 1);
	R_Pluto(1, 1) =  -1450832853776.53; 
    R_Pluto(2, 1) = -4318336722663.93;
    R_Pluto(3, 1) = -918296572141.713;
	
	Matrix R_Moon(3, 1);
	R_Moon(1, 1) =  -291608384.785849; 
    R_Moon(2, 1) = -266716833.450887;
    R_Moon(3, 1) = -76102486.0194983;
	
	Matrix R_Sun(3, 1);
	R_Sun(1, 1) =  26499033756.8212; 
    R_Sun(2, 1) = -132757417354.791;
    R_Sun(3, 1) = -57556718399.1903;
	
	_assert(m_equals(r_Mercury, R_Mercury, 1e-8));
	_assert(m_equals(r_Venus, R_Venus, 1e-8));
	_assert(m_equals(r_Earth, R_Earth, 1e-8));
	_assert(m_equals(r_Mars, R_Mars, 1e-8));
	_assert(m_equals(r_Jupiter, R_Jupiter, 1e-8));
	_assert(m_equals(r_Saturn, R_Saturn, 1e-8));
	_assert(m_equals(r_Uranus, R_Uranus, 1e-8));
	_assert(m_equals(r_Neptune, R_Neptune, 1e-8));
	_assert(m_equals(r_Pluto, R_Pluto, 1e-8));
	_assert(m_equals(r_Moon, R_Moon, 1e-8));
	_assert(m_equals(r_Sun, R_Sun, 1e-8));
	return 0;
}

int LTC_01() {
	double Mjd_1 = 51544.5;
	double Mjd_2 = 58000.0;

    Matrix P = LTC(Mjd_1,Mjd_2);
	Matrix R(3, 3);
    R(1, 1) = 0.3796; R(1, 2) = -0.9252; R(1, 3) = 0;
    R(2, 1) = -0.0772;  R(2, 2) = -0.0317; R(2, 3) = 0.9965;
    R(3, 1) = -0.9219;  R(3, 2) = -0.3782; R(3, 3) = -0.0835;
	_assert(m_equals(P, R, 1e-4));
	return 0;
}

int NutMatrix_01() {
	double Mjd_1 = 51544.5;

    Matrix P = NutMatrix(Mjd_1);
	Matrix R(3, 3);
    R(1, 1) =  1; R(1, 2) = 0.0001; R(1, 3) = 0;
    R(2, 1) = -0.0001;  R(2, 2) = 1; R(2, 3) = 0;
    R(3, 1) = 0;  R(3, 2) = 0; R(3, 3) = 1;
	_assert(m_equals(P, R, 1e-4));
	return 0;
}

int PoleMatrix_01() {
	double Mjd_1 = 51544.5;
    double Mjd_2 = 58000.0;

    Matrix P = PoleMatrix(Mjd_1, Mjd_2);
	Matrix R(3, 3);
    R(1, 1) =  -0.9252; R(1, 2) = 0.0317; R(1, 3) = -0.3782;
    R(2, 1) = 0;  R(2, 2) = 0.9965; R(2, 3) = 0.0835;
    R(3, 1) = 0.3796;  R(3, 2) = 0.0772; R(3, 3) = -0.9219;
	_assert(m_equals(P, R, 1e-4));
	return 0;
}

int PrecMatrix_01() {
	double Mjd_1 = 51544.5;
    double Mjd_2 = 58000.0;
	
	Matrix R(3, 3);
    R(1, 1) = 1; R(1, 2) = -0.0040; R(1, 3) = -0.0017;
    R(2, 1) = 0.0040;  R(2, 2) = 1; R(2, 3) = 0;
    R(3, 1) = 0.0017;  R(3, 2) = 0; R(3, 3) = 1;
	
	Matrix P = PrecMatrix(Mjd_1, Mjd_2);
	_assert(m_equals(P, R, 1e-4));
	return 0;
}

int gmst_01() {
	double R_gmst= 4.9065;
	double gmstime = gmst(59580.5);
	ASSERT_CLOSE(R_gmst,gmstime);
    return 0;
}

int gast_01() {
	double R_gast= 1.7998;
	double gastime = gast(37668);
	ASSERT_CLOSE(R_gast,gastime);
    return 0;
}

int MeasUpdate_01() {
	Matrix x = zeros(2, 1);
    x(1, 1) = 1.0;
    x(2, 1) = 0.5;

    Matrix z = zeros(1, 1);
    z(1, 1) = 1.2;

    Matrix g = zeros(1, 1);
    g(1, 1) = 1.0;

    Matrix s = zeros(1, 1);
    s(1, 1) = 0.1;

    Matrix G = zeros(1, 2);
    G(1, 1) = 1.0;
    G(1, 2) = 0.0;

    Matrix P = eye(2);
    P(1, 1) = 0.2;
    P(2, 2) = 0.2;

    int n = 2;
    auto [K, x_new, P_new] = MeasUpdate(x, z, g, s, G, P, n);
	
	Matrix R_K(2, 1);
    R_K(1, 1) =  0.9524;
    R_K(2, 1) = 0;
	
	Matrix R_x(2, 1);
    R_x(1, 1) =  1.1905;
    R_x(2, 1) = 0.5;
	
	Matrix R_P(2, 2);
    R_P(1, 1) =  0.0095; R_P(1, 2) = 0; 
    R_P(2, 1) = 0;  R_P(2, 2) = 0.2;
	
	_assert(m_equals(K, R_K, 1e-4));
	_assert(m_equals(x_new, R_x, 1e-4));
	_assert(m_equals(P_new, R_P, 1e-4));
    return 0;
}

int G_AccelHarmonic_01() {
	Matrix E = eye(3);
	int n_max = 10;
	int m_max = 10;
	   Matrix r(3);
	r(1,1) = 7000e3;  r(1,2) = 1000e3;  r(1,3) = 1300e3;
	
	Matrix A= G_AccelHarmonic(r, E, n_max, m_max);
	
	Matrix R=zeros(3,3);
	R(1, 1) = 1.98157940456412e-06; R(1, 2) = 4.36521842672732e-07; R(1, 3) = 5.69823020768467e-07;
    R(2, 1) = 4.36521843782955e-07;  R(2, 2) = -1.01138282149194e-06; R(2, 3) = 8.14115350689093e-08;
    R(3, 1) = 5.69823020102334e-07;  R(3, 2) = 8.14115350689093e-08; R(3, 3) = -9.70196583960359e-07;
	
	_assert(m_equals(A, R, 1e-4));
	return 0;
}

int GHAMatrix_01() {
	Matrix A= GHAMatrix(51544.5);
	
	Matrix R=zeros(3,3);
	R(1, 1) = 0.1815; R(1, 2) = -0.9834; 
    R(2, 1) = 0.9834;  R(2, 2) = 0.1815;
																	R(3, 3) = 1;
	_assert(m_equals(A, R, 1e-4));
    return 0;
}

int Accel_01() {
	AuxParam.Mjd_UTC = 58849.0;
    AuxParam.n = 10;
    AuxParam.m = 10;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;

	Matrix Y(6, 1);
    Y(1,1) = 7000000.0;
    Y(2,1) = 0.0;      
    Y(3,1) = 0.0;      
    Y(4,1) = 0.0;   
    Y(5,1) = 7500.0;
    Y(6,1) = 0.0;   

    double x = 0.0;
	
	Matrix dY = Accel(x, Y);
	Matrix R(6,1);
	R(1, 1) = 0; 
    R(2, 1) = 7500;
    R(3, 1) = 0;
	R(4, 1) = -8.14565001765738; 
    R(5, 1) = 4.97631036326491e-05;
    R(6, 1) = -5.40544911539266e-05;
	
	_assert(m_equals(dY, R, 1e-5));
	return 0;
}

int VarEqn_01() {
	AuxParam.Mjd_UTC = 58849.0;
    AuxParam.n = 10;
    AuxParam.m = 10;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;

	double x = 0.0;
    Matrix yPhi(42);
    yPhi(1) = 7000; yPhi(2) = 0; yPhi(3) = 0;
    yPhi(4) = 0; yPhi(5) = 7.5; yPhi(6) = 1.0;

    for (int j = 1; j <= 6; ++j) {
        for (int i = 1; i <= 6; ++i) {
            yPhi(6*j + i) = (i == j) ? 1.0 : 0.0;
        }
    }
	yPhi=transpose(yPhi);
    Matrix yPhip = VarEqn(x, yPhi);
	
	Matrix R(42,1);
	R(1, 1) = 0; 
    R(2, 1) = 7.5;
    R(3, 1) = 1;
	R(4, 1) =  1.57813230907696e+31; 
    R(5, 1) = -1.11848081633618e+31;
    R(6, 1) = -9.94090511338494e+29;
	R(7, 1) = 0;
	R(8, 1) = 0;
	R(9, 1) = 0;
	R(10, 1) = -2.70526998103899e+28;
	R(11, 1) = 1.91736618641012e+28;
	R(12, 1) = 1.70322436778075e+27;
	R(13, 1) = 0;
	R(14, 1) = 0;
	R(15, 1) = 0;
	R(16, 1) = 1.91736571660745e+28;
	R(17, 1) = 1.98525621889342e+28;
	R(18, 1) = 7.86077318297544e+26;
	R(19, 1) = 0;
	R(20, 1) = 0;
	R(21, 1) = 0;
	R(22, 1) = 1.70322462993725e+27;
	R(23, 1) = 7.8607776187776e+26;
	R(24, 1) = 7.20013043584687e+27;
	R(25, 1) = 1;
	R(26, 1) = 0;
	R(27, 1) = 0;
	R(28, 1) = 0;
	R(29, 1) = 0;
	R(30, 1) = 0;
	R(31, 1) = 0;
	R(32, 1) = 1;
	R(33, 1) = 0;
	R(34, 1) = 0;
	R(35, 1) = 0;
	R(36, 1) = 0;
	R(37, 1) = 0;
	R(38, 1) = 0;
	R(39, 1) = 1;
	R(40, 1) = 0;
	R(41, 1) = 0;
	R(42, 1) = 0;

	
	_assert(m_equals(yPhip, R, 1e-6));
	return 0;
}

int VarEqn_02() {
		 AuxParam.Mjd_UTC = 49746.1112847221;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;
	AuxParam.Mjd_TT =49746.1108586111;
	cout<<"-----------------------------------------"<<endl;
    int n = 42;
    double x = 0.0;
    Matrix yPhi = zeros(n,1);
    yPhi(1,1) = 554255.93722869; yPhi(2,1) = 3213514.86734919; yPhi(3,1) = 3990892.97587674;
    yPhi(4,1) = 5394.06842166295; yPhi(5,1) = -2365.21337882319; yPhi(6,1) = -7061.84554200204;
    yPhi(7,1) = 1.0; yPhi(14,1) = 1.0; yPhi(21,1) = 1.0; yPhi(28,1) = 1.0; yPhi(35,1) = 1.0; yPhi(42,1) = 1.0;

    Matrix yPhip = VarEqn(x,yPhi);
	cout<<"yPhip"<<endl<<yPhip<<endl;
	
    Matrix yPhip_sol = zeros(n,1);
    yPhip_sol(1,1) = 5394.06842166295;   yPhip_sol(2,1) = -2365.21337882319;   yPhip_sol(3,1) = -7061.84554200204;
    yPhip_sol(4,1) = -5.13438675840588;  yPhip_sol(5,1) = -2.9771762235621;    yPhip_sol(6,1) = -3.70591776714193;
    yPhip_sol(10,1) = 5.70032084907797e-07;  yPhip_sol(11,1) = 8.67651592351137e-07;  yPhip_sol(12,1) = 1.08169353918441e-06;
    yPhip_sol(16,1) = 8.6765159574781e-07;   yPhip_sol(17,1) = -4.2335910770693e-07;  yPhip_sol(18,1) = 6.2718370418322e-07;
    yPhip_sol(22,1) = 1.08169354007259e-06;  yPhip_sol(23,1) = 6.2718370490946e-07;   yPhip_sol(24,1) = -1.46672925360747e-07;
    yPhip_sol(25,1) = 1.0; yPhip_sol(32,1) = 1.0; yPhip_sol(39,1) = 1.0;
	
	cout<<"yPhip_sol"<<endl<<yPhip_sol<<endl;
		
	Matrix yPhip_R = zeros(n,1);
    yPhip_R(1,1) = 5394.06842166295;   yPhip_R(2,1) = -2365.21337882319;   yPhip_R(3,1) = -7061.84554200204;
    yPhip_R(4,1) = -1.60604427421658;  yPhip_R(5,1) = -9.31397009489115;    yPhip_R(6,1) = -11.6198734901516;
    yPhip_R(10,1) = -2.79882609643956e-06;  yPhip_R(11,1) =5.76876576730001e-07;  yPhip_R(12,1) = 7.30575216323359e-07;
    yPhip_R(16,1) = 5.76876580948849e-07;   yPhip_R(17,1) = 4.70193322499313e-07;  yPhip_R(18,1) = 4.20745706186665e-06;
    yPhip_R(22,1) = 7.3057521876585e-07;  yPhip_R(23,1) = 4.20745705653758e-06;   yPhip_R(24,1) = 2.32863277815909e-06;
    yPhip_R(25,1) = 1.0; yPhip_sol(32,1) = 1.0; yPhip_R(39,1) = 1.0;

		cout<<"yPhip_R"<<endl<<yPhip_R<<endl;
	
    _assert(m_equals(yPhip_R, yPhip, 1e-10));
    return 0;
}


int DEInteg_01() {
    Matrix Y0_apr(6);

    Y0_apr(1) = 6221397.62857869;
	Y0_apr(2) = 2867713.77965738;
	Y0_apr(3) = 3006155.98509949;
    Y0_apr(4) = 4645.04725161807;
	Y0_apr(5) = -2752.21591588205;
	Y0_apr(6) = -7507.99940987033; 
	Y0_apr=transpose(Y0_apr);
    
    Matrix Y = DEInteg(Accel,0,-134.999991953373,1e-13,1e-6,6,Y0_apr);
	
	Matrix R=zeros(6,1);
	R(1, 1) = 5542555.93722861; 
    R(2, 1) = 3213514.8673492;
    R(3, 1) = 3990892.97587686;
	R(4, 1) =  5394.06842166353; 
    R(5, 1) = -2365.21337882342;
    R(6, 1) = -7061.84554200298;
	
	_assert(m_equals(Y, R, 1e-8));
	return 0;
}

int all_tests()
{
	_verify(m_sum_01);
	_verify(m_sub_01);
	_verify(m_mult_01);
	_verify(m_div_01);
	_verify(m_assign_01);
	_verify(m_zeros_01);
	_verify(m_inv_01);
	_verify(m_eye_01);
	_verify(m_transpose_01);
	_verify(op_sum_01);
	_verify(op_sub_01);
	_verify(op_mult_01);
	_verify(op_div_01);
	_verify(m_norm_01);
	_verify(m_dot_01);
	_verify(m_cross_01);
	_verify(extract_vector_01);
	_verify(union_vector_01);
	_verify(m_extract_row_01);
	_verify(m_extract_col_01);
	_verify(m_assign_row_01);
	_verify(m_assign_col_01);
	_verify(R_x_01);
	_verify(R_y_01);
	_verify(R_z_01);
	_verify(AccelPointMass_01);
	_verify(Cheb3D_01);
	_verify(EccAnom_01);
	_verify(Frac_01);
	_verify(MeanObliquity_01);
	_verify(Mjday_01);
	_verify(Mjday_TDB_01);
	_verify(Position_01);
	_verify(SAT_Const_01);
	_verify(sign__01);
	_verify(timediff_01);
	_verify(AzElPa_01);
	_verify(IERS_01);
	_verify(Legendre_01);
	_verify(NutAngles_01);
	_verify(TimeUpdate_01);
	_verify(AccelHarmonic_01);
	_verify(EqnEquinox_01);
	_verify(JPL_Eph_DE430_01);
	_verify(LTC_01);
	_verify(NutMatrix_01);
	_verify(PoleMatrix_01);
	_verify(PrecMatrix_01);
	_verify(gmst_01);
	_verify(gast_01);
	_verify(MeasUpdate_01);
	_verify(G_AccelHarmonic_01);
	_verify(GHAMatrix_01);
	//_verify(DEInteg_01);
	_verify(Accel_01);
	_verify(VarEqn_01);
	//_verify(VarEqn_02);
	
    return 0;
}


int main()
{
	eop19620101(21413);
	GGM03S(181);
	DE430Coeff(2285,1020);
	 AuxParam.Mjd_UTC = 49746.1112847221;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;
	
    int result = all_tests();
	

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
