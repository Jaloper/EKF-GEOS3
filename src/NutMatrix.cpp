// $Source$
//--------------------------------------------------------------------------------
// NutMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file NutMatrix.cpp
 *  @brief Implements the computation of the nutation transformation matrix.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\NutMatrix.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_z.hpp"

Matrix NutMatrix(double Mjd_TT) {
	
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity(Mjd_TT);

    // Nutation in longitude and obliquity
    auto [dpsi, deps] = NutAngles(Mjd_TT);

    // Transformation from mean to true equator and equinox
    Matrix NutMat = R_x(-(eps + deps)) * R_z(-dpsi) * R_x(eps);

    return NutMat;
}
