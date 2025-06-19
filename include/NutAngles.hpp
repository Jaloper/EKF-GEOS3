// $Header$
//--------------------------------------------------------------------------------
// NutAngles
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/29
//
/** @file NutAngles.hpp
 *  @brief This header file declares a function to compute nutation angles in longitude and obliquity.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _NUTANGLES_
#define _NUTANGLES_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// NutAngles(double Mjd_TT)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the nutation angles in longitude and obliquity for a given Modified Julian Date in Terrestrial Time.
 *
 *  @param [in] Mjd_TT Modified Julian Date in Terrestrial Time (TT).
 *  @return std::tuple<double, double> Tuple containing the nutation in longitude (dpsi) and obliquity (deps) in radians.
 */
std::tuple<double, double> NutAngles(double Mjd_TT);

#endif