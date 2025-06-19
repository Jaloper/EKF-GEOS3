// $Header$
//--------------------------------------------------------------------------------
// MeanObliquity
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/17
//
/** @file MeanObliquity.hpp
 *  @brief This header file declares a function to compute the mean obliquity of the ecliptic.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _MEANOBLIQUITY_
#define _MEANOBLIQUITY_
#define _USE_MATH_DEFINES

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// MeanObliquity(double Mjd_TT)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the mean obliquity of the ecliptic for a given Modified Julian Date in Terrestrial Time.
 *
 *  @param [in] Mjd_TT Modified Julian Date in Terrestrial Time (TT).
 *  @return double The mean obliquity of the ecliptic in radians.
 */
double MeanObliquity(double Mjd_TT);

#endif