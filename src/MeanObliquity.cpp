#include "..\include\MeanObliquity.hpp"


double MeanObliquity(double Mjd_TT) {
    const double MJD_J2000 = 51544.5;
    const double Rad = M_PI / 180.0;

    double T = (Mjd_TT - MJD_J2000) / 36525.0;

    double MOblq = Rad * (84381.448 / 3600.0 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0);

    return MOblq;
}
