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
#include <iostream>

using namespace std;

int main() {
    eop19620101(4); //c=21413
    cout << "eopdata\n" << eopdata << "\n";

    Matrix s(3);
    s(1) = 1000.0;
    s(2) = 2000.0;
    s(3) = 3000.0;

    auto result = AzElPa(s);
// Llamada a la función TimeUpdate
Matrix P_updated = TimeUpdate(P, Phi,Qdt);

    
    
	cout << "Phi\n" << Phi << "\n";
	cout << "P\n" << P << "\n";
	cout << "Qdt\n" << Qdt << "\n";
	cout << "P_updated\n" << P_updated << "\n";




    return 0;
}
