// $Header$
//--------------------------------------------------------------------------------
// AccelPointMass
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/09
//
/** @file AccelPointMass.hpp
 *  @brief This header declares the function to compute the acceleration on a satellite
 *         due to a third-body gravitational perturbation (e.g., Sun, Moon, or planet).
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _ACCELPOINTMASS_
#define _ACCELPOINTMASS_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// AccelPointMass(Matrix r, Matrix s, double GM)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the perturbing acceleration on a satellite due to a third body
 *         using Newton's law of universal gravitation.
 *
 *  @param [in] r   Position vector of the satellite in inertial frame [m].
 *  @param [in] s   Position vector of the third body (e.g., Sun or Moon) [m].
 *  @param [in] GM  Gravitational parameter of the third body [m³/s²].
 *
 *  @return Matrix  Acceleration vector on the satellite due to the third body [m/s²].
 */
//-----------------------------------------------------------------------------------------------
Matrix AccelPointMass(Matrix r, Matrix s, double GM);

#endif