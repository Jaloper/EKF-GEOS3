// $Header$
//--------------------------------------------------------------------------------
// G_AccelHarmonic
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/16
//
/** @file G_AccelHarmonic.hpp
 *  @brief This header file declares the function to calculate the gradient of the harmonic acceleration.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _G_ACCELHARMONIC_
#define _G_ACCELHARMONIC_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the gradient of the harmonic acceleration with respect to position.
 *
 *  This function estimates the partial derivatives of the harmonic acceleration vector
 *  using central differences, forming a 3x3 gradient matrix. Each column of the output
 *  corresponds to the derivative of the acceleration vector with respect to one spatial axis.
 *
 *  @param [in] r       Position vector (3x1) in an Earth-centered coordinate system.
 *  @param [in] U       Gravitational potential coefficients matrix.
 *  @param [in] n_max   Maximum degree of the spherical harmonics expansion.
 *  @param [in] m_max   Maximum order of the spherical harmonics expansion.
 *  @return Matrix      A 3x3 matrix representing the gradient of the acceleration field.
 */
Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max);

#endif
