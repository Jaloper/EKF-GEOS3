// $Header$
//--------------------------------------------------------------------------------
// Position
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/18
//
/** @file Position.hpp
 *  @brief This header file declares a function to compute the geocentric position vector from geodetic coordinates.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _POSITION_
#define _POSITION_

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// Position(double lon, double lat, double h)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the geocentric position vector based on geodetic longitude, latitude, and height above the ellipsoid.
 *
 *  @param [in] lon Longitude in radians.
 *  @param [in] lat Geodetic latitude in radians.
 *  @param [in] h   Height above the ellipsoid in meters.
 *  @return Matrix A 3x1 position vector (x, y, z) in meters in the geocentric reference frame.
 */
Matrix Position(double lon, double lat, double h);

#endif