#include "..\include\AccelPointMass.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\global.hpp"

Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max) {
	
    const double r_ref = 6378.1363e3;         // [m] Earth's reference radius
    const double gm = 398600.4415e9;          // [m^3/s^2] gravitational constant

    // Body-fixed position 
    Matrix r_bf = transpose(E * transpose(r));
	
    // Auxiliary quantities
    double d = norm(r_bf);                     // distance
    double latgc = asin(r_bf(3)/d);
    double lon = atan2(r_bf(2),r_bf(1));
	
    auto [pnm, dpnm] = Legendre(n_max, m_max, latgc);
    double dUdr = 0.0, dUdlatgc = 0.0, dUdlon = 0.0;
    double q1 = 0, q2 = 0, q3 = 0;

    for (int n = 0; n <= n_max; ++n) {
        double b1 = (-gm / (d * d)) * pow(r_ref / d, n) * (n + 1);
        double b2 =  (gm / d) * pow(r_ref / d, n);
        double b3 =  (gm / d) * pow(r_ref / d, n);
        for (int m = 0; m <= m_max; ++m) {
            q1 = q1 + pnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
			q2 = q2 + dpnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
			q3 = q3 + m*pnm(n+1,m+1)*(Snm(n+1,m+1)*cos(m*lon)-Cnm(n+1,m+1)*sin(m*lon));
        }
         dUdr     = dUdr     + q1*b1;
		dUdlatgc = dUdlatgc + q2*b2;
		dUdlon   = dUdlon   + q3*b3;
        q1 = q2 = q3 = 0.0;
    }

    // Body-fixed acceleration
    double r2xy = r_bf(1) * r_bf(1) + r_bf(2) * r_bf(2);
    double sqrt_r2xy = std::sqrt(r2xy);
	double ax = (1.0/d*dUdr - r_bf(3)/(pow(d,2)*sqrt(r2xy)) * dUdlatgc) * r_bf(1) - (1/r2xy * dUdlon) * r_bf(2);

	double ay = (1.0/d*dUdr-r_bf(3)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1/r2xy*dUdlon)*r_bf(1);
	double az =  1.0/d*dUdr*r_bf(3)+sqrt(r2xy)/pow(d,2)*dUdlatgc;
	
    Matrix a_bf=zeros(3,3);
	a_bf(1,1)=ax; a_bf(1,2) = ay;   a_bf(1,3) = az;
	a_bf=transpose(a_bf);
    // Inertial acceleration
	Matrix aux=transpose(E) * a_bf;
    return extract_column(aux,1);
}
