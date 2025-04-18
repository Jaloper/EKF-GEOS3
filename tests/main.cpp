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
	double lon = 0.1745;  // Aproximadamente 10 grados en radianes
    double lat = 0.7854;  // Aproximadamente 45 grados en radianes
    double h = 100.0;     // Altura en metros

    double R_Earth = 6378137.0;                 // Radio ecuatorial de la Tierra
    double f_Earth = 1.0 / 298.257223563;       // Achatamiento de la Tierra

    // Llamamos a la función Position
    Matrix M6 = Position(lon, lat, h);
    std::cout << "Radio de la Tierra (R_Earth): " << SAT_Const::R_Earth << " m" << std::endl;
    std::cout << "1 grado en radianes (Rad): " << SAT_Const::Rad << std::endl;

    double E = Mjday_TDB(51544.5);
	cout << "E\n" << E << "\n";
    cout << "M1\n" << M1 << "\n";
    cout << "M2\n" << M2 << "\n";
    cout << "M3\n" << M3 << "\n";
	cout << "M4\n" << M4 << "\n";
	cout << "M5\n" << M5 << "\n";
	cout << "M6\n"<< M6 << "\n";
    return 0;
}