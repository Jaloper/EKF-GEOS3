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
#include <iostream>

using namespace std;

int main() {
    eop19620101(21413); //c=21413
	GGM03S(181);
	DE430Coeff(2285,1020);
    //cout << "eopdata\n" << eopdata << "\n";
 AuxParam.Mjd_UTC = 58849.0;
    AuxParam.n = 10;
    AuxParam.m = 10;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;

 double x = 0.0;

    // Crear vector yPhi de 42 elementos
    Matrix yPhi(42); // Vector columna

    // Estado inicial: posición y velocidad
    yPhi(1) = 7000; yPhi(2) = 0; yPhi(3) = 0;       // r
    yPhi(4) = 0; yPhi(5) = 7.5; yPhi(6) = 1.0;      // v

    // Matriz identidad de 6x6 en formato column-wise
    for (int j = 1; j <= 6; ++j) {
        for (int i = 1; i <= 6; ++i) {
            yPhi(6*j + i) = (i == j) ? 1.0 : 0.0;
        }
    }

    // Ejecutar la función
    Matrix yPhip = VarEqn(x, yPhi);

    // Imprimir resultado
    std::cout <<"yPhip\n"<< yPhip<<endl;


    return 0;
}
