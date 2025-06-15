#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\AccelHarmonic.hpp"

Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max) {

	double d = 1.0;   // Position increment [m]

	Matrix G = zeros(3,3);
	Matrix dr = zeros(1,3);

// Gradient
	for (int i = 1; i <= 3; ++i) {
		 for (int j = 1; j <= 3; ++j) dr(1, j) = 0.0;
		// Set offset in i-th component of the position vector
		 dr(1, i) = d;
		// Acceleration difference
		Matrix accel1= AccelHarmonic ( r+dr/2,U, n_max, m_max );
		Matrix accel2= AccelHarmonic ( r-dr/2,U, n_max, m_max );
		Matrix da =  accel1-accel2;
		// Derivative with respect to i-th axis
		for (int j = 1; j <= 3; ++j) {
				G(j, i) = da(j, 1) / d;
			}  		
	}
	return G;
}
