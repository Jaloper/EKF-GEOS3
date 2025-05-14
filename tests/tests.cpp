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
	
	Matrix A = zeros(2, 4);
	A(1,1) = 4; A(1,2) = 7;
	A(2,1) = 2; A(2,2) = 6; A(2,3) = 2; A(2,4) = 6;
	
	Matrix B= zeros(3,1);
	B(1,1) = 1; 
	B(2,1) = 2; 
	B(3,1) = 3;
	
	Matrix C =  assign_column(A,B,2);
	
	Matrix R= zeros(3,4);
	R(1,1) = 4; R(1,2) = 1;
	R(2,1) = 2; R(2,2) = 2; R(2,3) = 2; R(2,4) = 6;
				R(3,2) = 3;

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

    Matrix A = timediff(UT1_UTC, TAI_UTC);
	Matrix R(5);
	R(1,1) = -36.665900; R(1,2) = -18; R(1,3) = -17.665900;R(1,4) = 69.184000; R(1,5) = 18;
	
    _assert(m_equals(A, R, 1e-4));
    
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
	
	Matrix R_pnm(4,4);
    R_pnm(1,1) = 1; R_pnm(1,2) = 0; R_pnm(1,3) = 0; R_pnm(1,4) = 0;
	R_pnm(2,1) = 1.4575; R_pnm(2,2) = 0.9358; R_pnm(2,3) = 0; R_pnm(2,4) = 0;
	R_pnm(3,1) = 1.2569; R_pnm(3,2) = 1.7608; R_pnm(3,3) = 0.5653; R_pnm(3,4) = 0;
	
	Matrix R_dpnm(4,4);
    R_dpnm(1,1) = 0; R_dpnm(1,2) = 0; R_dpnm(1,3) = 0; R_dpnm(1,4) = 0;
	R_dpnm(2,1) = 0.9358; R_dpnm(2,2) = -1.4575; R_dpnm(2,3) = 0; R_dpnm(2,4) = 0;
	R_dpnm(3,1) = 3.0499; R_dpnm(3,2) = -1.6117; R_dpnm(3,3) = -1.7608; R_dpnm(3,4) = 0;
	
	_assert(m_equals(dpnm, R_dpnm, 1e-4));
	_assert(m_equals(pnm, R_pnm, 1e-4));
	
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
	
	Matrix R=zeros(3,3);
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
	
	_assert(m_equals(r_Mercury, R_Mercury, 1e-4));
	_assert(m_equals(r_Venus, R_Venus, 1e-4));
	_assert(m_equals(r_Earth, R_Earth, 1e-4));
	_assert(m_equals(r_Mars, R_Mars, 1e-4));
	_assert(m_equals(r_Jupiter, R_Jupiter, 1e-4));
	_assert(m_equals(r_Saturn, R_Saturn, 1e-4));
	_assert(m_equals(r_Uranus, R_Uranus, 1e-4));
	_assert(m_equals(r_Neptune, R_Neptune, 1e-4));
	_assert(m_equals(r_Pluto, R_Pluto, 1e-4));
	_assert(m_equals(r_Moon, R_Moon, 1e-4));
	_assert(m_equals(r_Sun, R_Sun, 1e-4));
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
	
    return 0;
}


int main()
{
	eop19620101(21413);
	GGM03S(181);
	DE430Coeff(2285,1020);
    int result = all_tests();
	

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
