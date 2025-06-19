// $Header$
//--------------------------------------------------------------------------------
// gmst
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/07
//
/** @file gmst.hpp
 *  @brief This header file declares the function that computes the Greenwich Mean Sidereal Time (GMST).
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _GMST_
#define _GMST_
#define _USE_MATH_DEFINES

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// gmst(double Mjd_UT1)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the Greenwich Mean Sidereal Time (GMST) for a given UT1 Julian Date.
 *
 *  This function calculates GMST based on the IAU 2006 expression, normalizing the result
 *  into the range [0, 2π) radians.
 *
 *  @param [in] Mjd_UT1  Modified Julian Date in UT1 time scale.
 *  @return double       GMST in radians, within the range [0, 2π).
 */
double gmst(double Mjd_UT1);

#endif
