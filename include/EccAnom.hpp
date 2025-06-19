// $Header$
//--------------------------------------------------------------------------------
// EccAnom
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/17
//
/** @file EccAnom.hpp
 *  @brief This header file declares the function to calculate the eccentric anomaly from the mean anomaly and eccentricity using Newton's method.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _ECCANOM_
#define _ECCANOM_
#define _USE_MATH_DEFINES

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// EccAnom(double M, double e)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Calculates the eccentric anomaly using Newton's method based on the mean anomaly and eccentricity.
 *
 *  @param [in] M  Mean anomaly (in radians).
 *  @param [in] e  Orbital eccentricity.
 *
 *  @return double  The calculated eccentric anomaly (in radians).
 */
double EccAnom(double M, double e);

#endif
