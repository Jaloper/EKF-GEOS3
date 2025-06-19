// $Header$
//--------------------------------------------------------------------------------
// AccelHarmonic
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/07
//
/** @file AccelHarmonic.hpp
 *  @brief This header declares the function to compute the acceleration due to
 *         Earth's harmonic gravity field in the inertial frame.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the acceleration of a satellite due to the Earth's gravitational
 *         potential using a spherical harmonic expansion.
 *
 *  @param [in] r      Satellite position vector in inertial frame [m].
 *  @param [in] E      Transformation matrix from inertial to body-fixed frame.
 *  @param [in] n_max  Maximum degree of the spherical harmonics.
 *  @param [in] m_max  Maximum order of the spherical harmonics.
 *
 *  @return Matrix     Acceleration vector in inertial frame [m/s^2].
 */
//-----------------------------------------------------------------------------------------------
Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max);

#endif