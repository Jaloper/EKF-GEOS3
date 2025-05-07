#include "..\include\IERS.hpp"
#include "..\include\SAT_Const.hpp"
#include <stdexcept>
#include <tuple>

std::tuple<double, double, double, double, double, double, double, double, double> 
IERS(
    Matrix eop, double Mjd_UTC, char interp
) {

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;

    if (interp == 'l') {
        // linear interpolation
        int mjd = static_cast<int>(floor(Mjd_UTC));
        int i = -1;

        for (int col = 1; col <= eop.n_column; ++col) {
            if (static_cast<int>(eop(4, col)) == mjd) {
                i = col;
                break;
            }
        }

        if (i == -1 || i + 1 > eop.n_column) {
            throw std::runtime_error("MJD not found or no next data point for interpolation");
        }

        double mfme = 1440.0 * (Mjd_UTC - floor(Mjd_UTC));
        double fixf = mfme / 1440.0;

        x_pole  = eop(5,i) + (eop(5,i+1) - eop(5,i)) * fixf;
        y_pole  = eop(6,i) + (eop(6,i+1) - eop(6,i)) * fixf;
        UT1_UTC = eop(7,i) + (eop(7,i+1) - eop(7,i)) * fixf;
        LOD     = eop(8,i) + (eop(8,i+1) - eop(8,i)) * fixf;
        dpsi    = eop(9,i) + (eop(9,i+1) - eop(9,i)) * fixf;
        deps    = eop(10,i) + (eop(10,i+1) - eop(10,i)) * fixf;
        dx_pole = eop(11,i) + (eop(11,i+1) - eop(11,i)) * fixf;
        dy_pole = eop(12,i) + (eop(12,i+1) - eop(12,i)) * fixf;
        TAI_UTC = eop(13,i);

        x_pole  /= SAT_Const::Arcs;
        y_pole  /= SAT_Const::Arcs;
        dpsi    /= SAT_Const::Arcs;
        deps    /= SAT_Const::Arcs;
        dx_pole /= SAT_Const::Arcs;
        dy_pole /= SAT_Const::Arcs;
    }
    else if (interp == 'n') {
        int mjd = static_cast<int>(floor(Mjd_UTC));
        int i = -1;

        for (int col = 1; col <= eop.n_column; ++col) {
            if (static_cast<int>(eop(4, col)) == mjd) {
                i = col;
                break;
            }
        }

        if (i == -1) {
            throw std::runtime_error("MJD not found");
        }

        x_pole  = eop(5,i) / SAT_Const::Arcs;
        y_pole  = eop(6,i) / SAT_Const::Arcs;
        UT1_UTC = eop(7,i);
        LOD     = eop(8,i);
        dpsi    = eop(9,i) / SAT_Const::Arcs;
        deps    = eop(10,i) / SAT_Const::Arcs;
        dx_pole = eop(11,i) / SAT_Const::Arcs;
        dy_pole = eop(12,i) / SAT_Const::Arcs;
        TAI_UTC = eop(13,i);
    }

    return std::make_tuple(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}
