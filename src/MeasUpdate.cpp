// $Source$
//--------------------------------------------------------------------------------
// MeasUpdate
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/16
//
/** @file MeasUpdate.cpp
 *  @brief Implements the measurement update step for a Kalman filter.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\MeasUpdate.hpp"

std::tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix x, Matrix z, Matrix g, Matrix s, Matrix G, Matrix P, int n){
	int m;
	if (z.n_row>z.n_column)  m= z.n_row;
	else m=z.n_column;
	Matrix& Inv_W = zeros(m,m);

	for (int i = 1; i <= m; ++i) {
		Inv_W(i,i) = s(i)*s(i);    // Inverse weight (measurement covariance)
	}

	// Kalman gain
	Matrix& Gt = transpose(G);
    Matrix& GP = G * P;
    Matrix& S = Inv_W + GP * Gt;
    Matrix& Sinv = inv(S);
    Matrix& K = P * Gt * Sinv;

	// State update
	Matrix& x_dev = x + K*(z-g);

	// Covariance update
	Matrix& P_dev = (eye(n)-K*G)*P;
	
	
	return std::tie(K, x_dev, P_dev);
}
