// $Source$
//--------------------------------------------------------------------------------
// LTC
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file LTC.cpp
 *  @brief Implements the computation of the Local Tangent Coordinate transformation matrix.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

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
