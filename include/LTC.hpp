// $Header$
//--------------------------------------------------------------------------------
// LTC
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file LTC.hpp
 *  @brief This header file declares a function to compute the Local Tangent Coordinate transformation matrix.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _LTC_
#define _LTC_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// LTC(double lon, double lat)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the transformation matrix from geocentric to Local Tangent Coordinates (LTC) based on longitude and latitude.
 *
 *  @param [in] lon Longitude in radians.
 *  @param [in] lat Latitude in radians.
 *  @return Matrix  A 3x3 transformation matrix representing the rotation from geocentric to LTC system.
 */
Matrix LTC(double lon, double lat);

#endif