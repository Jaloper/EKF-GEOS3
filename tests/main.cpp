#include "..\include\matrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\AccelPointMass.hpp"
#include <iostream>

using namespace std;

int main() {
    Matrix M1(2, 4);
	M1(1,1) = 4; M1(1,2) = 7;
	M1(2,1) = 2; M1(2,2) = 6; M1(2,3) = 2; M1(2,4) = 6;
	
    Matrix M2(1,2);
	M2(1,1) = 2; M2(1,2) = 2;

    Matrix M3(3,1);
	M3(1,1) = 1; 
	M3(2,1) = 2; 
	M3(3,1) = 3;
	
	Matrix M4(1,2);
	M4(1,1) = 2; M4(1,2) = 1;
	
	Matrix M5 =  union_vector(M4,M2); 
	Matrix b = AccelPointMass(M2,M4,3.2);
    cout << "M1\n" << M1 << "\n";
    cout << "M2\n" << M2 << "\n";
    cout << "M3\n" << M3 << "\n";
	cout << "M4\n" << M4 << "\n";
	cout << "M5\n" << M5 << "\n";
	cout << b << "\n";
    return 0;
}