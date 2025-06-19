// $Header$
//--------------------------------------------------------------------------------
// gast
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/16
//
/** @file gast.hpp
 *  @brief This header file declares the function to compute the Greenwich Apparent Sidereal Time (GAST).
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _GAST_
#define _GAST_
#define _USE_MATH_DEFINES

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// gast(double Mjd_UT1)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the Greenwich Apparent Sidereal Time (GAST) for a given UT1 Julian Date.
 *
 *  This function combines the Greenwich Mean Sidereal Time (GMST) and the equation of the equinoxes
 *  to calculate the apparent sidereal time, normalized within [0, 2π).
 *
 *  @param [in] Mjd_UT1  Modified Julian Date in UT1 time scale.
 *  @return double       GAST in radians, in the range [0, 2π).
 */
double gast(double Mjd_UT1);

#endif