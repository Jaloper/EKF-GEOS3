#include "..\include\LTC.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"

Matrix LTC(double lon, double lat) {
	
    Matrix M = R_y(-1.0 * lat) * R_z(lon);

    for (int j = 1; j <= 3; j++) {
        double aux = M(1, j);
        M(1, j) = M(2, j);
        M(2, j) = M(3, j);
        M(3, j) = aux;
    }

    return M;
}
