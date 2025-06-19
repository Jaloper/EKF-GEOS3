// $Header$
//--------------------------------------------------------------------------------
// EqnEquinox
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file EqnEquinox.hpp
 *  @brief This header file declares the function to calculate the equation of the equinoxes based on the nutation in longitude and the mean obliquity.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _EQNEQUINOX_
#define _EQNEQUINOX_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// EqnEquinox(double Mjd_TT)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Calculates the equation of the equinoxes based on the nutation in longitude and the mean obliquity.
 *
 *  @param [in] Mjd_TT  Modified Julian Date in Terrestrial Time (TT).
 *
 *  @return double  The equation of the equinoxes (in radians).
 */
double EqnEquinox(double Mjd_TT);

#endif
