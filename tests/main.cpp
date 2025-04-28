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
#include <iostream>

using namespace std;

int main() {
    eop19620101(4); //c=21413
    cout << "eopdata\n" << eopdata << "\n";

    Matrix M1(2, 4);
    M1(1, 1) = 4; M1(1, 2) = 7;
    M1(2, 1) = 2; M1(2, 2) = 6; M1(2, 3) = 2; M1(2, 4) = 6;

    Matrix M2(5);
    M2(1) = 1; M2(2) = -0.5; M2(3) = 0.25; M2(4) = -0.1; M2(5) = 0.05;

    Matrix M3(3, 1);
    M3(1, 1) = 1;
    M3(2, 1) = 2;
    M3(3, 1) = 3;

    Matrix M4(5);
    M4(1) = 0.5; M4(2) = 0.3; M4(3) = -0.15; M4(4) = 0.1; M4(5) = -0.05;

    Matrix M5(5);
    M5(1) = 0.2; M5(2) = -0.1; M5(3) = 0.1; M5(4) = -0.05; M5(5) = 0.02;
    Matrix s(3);
    s(1) = 1000.0;
    s(2) = 2000.0;
    s(3) = 3000.0;

auto result = AzElPa(s);

// Acceder a los valores de la tupla
double Az = std::get<0>(result);
double El = std::get<1>(result);
Matrix dAds = std::get<2>(result);
Matrix dEds = std::get<3>(result);

    cout << "M1\n" << M1 << "\n";
    cout << "M2\n" << M2 << "\n";
    cout << "M3\n" << M3 << "\n";
    cout << "M4\n" << M4 << "\n";
    cout << "M5\n" << M5 << "\n";
    cout << "dAds\n" << dAds << "\n";
    cout << "dEds\n" << dEds << "\n";
    cout << "El\n" << El << "\n";
    cout << "Az\n" << Az << "\n";




    return 0;
}
