// $Header$
//--------------------------------------------------------------------------------
// AzElPa
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/18
//
/** @file AzElPa.hpp
 *  @brief This header file declares the function to calculate the azimuth, elevation, 
 *         and their partial derivatives with respect to the position vector.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _AZELPA_
#define _AZELPA_
#define _USE_MATH_DEFINES

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// AzElPa(Matrix s)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Calculates azimuth, elevation, and their partial derivatives with respect to position vector.
 *
 *  @param [in] s  Position vector in Cartesian coordinates (3x1).
 *
 *  @return std::tuple<double, double, Matrix&, Matrix&>  Tuple with azimuth (Az), 
 *                                                         elevation (El), and their partial 
 *                                                         derivatives (dAds, dEds).
 */
std::tuple<double, double, Matrix&, Matrix&> AzElPa(Matrix s);

#endif