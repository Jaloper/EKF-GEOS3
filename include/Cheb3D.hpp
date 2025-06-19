// $Header$
//--------------------------------------------------------------------------------
// Cheb3D
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/11
//
/** @file Cheb3D.hpp
 *  @brief This header file declares the function to calculate the 3D Chebyshev interpolation based on the input coefficients.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _CHEB3D_
#define _CHEB3D_

#include "..\include\matrix.hpp"
#include <cmath>


//-----------------------------------------------------------------------------------------------
// Cheb3D(double t, int N, double Ta, double Tb, Matrix Cx, Matrix Cy, Matrix Cz)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the 3D Chebyshev interpolation using the provided coefficients for each axis.
 *
 *  @param [in] t  Time for which to calculate the interpolated value.
 *  @param [in] N  Number of Chebyshev terms (degree of interpolation).
 *  @param [in] Ta  Start time of the interpolation range.
 *  @param [in] Tb  End time of the interpolation range.
 *  @param [in] Cx  Matrix of coefficients for the X-axis.
 *  @param [in] Cy  Matrix of coefficients for the Y-axis.
 *  @param [in] Cz  Matrix of coefficients for the Z-axis.
 *
 *  @return Matrix  The interpolated 3D position as a matrix [x, y, z].
 */
 Matrix Cheb3D(double t, int N, double Ta,double Tb, Matrix Cx, Matrix Cy, Matrix Cz);

#endif