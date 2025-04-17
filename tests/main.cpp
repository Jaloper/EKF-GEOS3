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
#include <iostream>

using namespace std;

int main() {
    Matrix M1(2, 4);
	M1(1,1) = 4; M1(1,2) = 7;
	M1(2,1) = 2; M1(2,2) = 6; M1(2,3) = 2; M1(2,4) = 6;
	
    Matrix M2(5);
	M2(1) = 1; M2(2) = -0.5; M2(3) = 0.25; M2(4) = -0.1; M2(5) = 0.05;

    Matrix M3(3,1);
	M3(1,1) = 1; 
	M3(2,1) = 2; 
	M3(3,1) = 3;
	
	Matrix M4(5);
	M4(1) = 0.5; M4(2) = 0.3; M4(3) = -0.15; M4(4) = 0.1; M4(5) = -0.05;
	
	Matrix M5(5);
	M5(1) = 0.2; M5(2) = -0.1; M5(3) = 0.1; M5(4) = -0.05; M5(5) = 0.02; 
	Matrix M6 = Cheb3D(0.5, 5, 0, 1, M2, M4, M5);
	
	int yr = 2025, mon = 4, day = 17, hr = 12, min = 30;
    double sec = 0.0;

    double E = Mjday(yr, mon, day, hr, min, sec);
	cout << "E\n" << E << "\n";
    cout << "M1\n" << M1 << "\n";
    cout << "M2\n" << M2 << "\n";
    cout << "M3\n" << M3 << "\n";
	cout << "M4\n" << M4 << "\n";
	cout << "M5\n" << M5 << "\n";
	cout << "M6\n"<< M6 << "\n";
    return 0;
}