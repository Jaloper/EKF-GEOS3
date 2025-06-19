// $Header$
//--------------------------------------------------------------------------------
// Accel
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/16
//
/** @file Accel.hpp
 *  @brief This header file declares de function to compute satellite acceleration due to Earth's gravity and third-body perturbations.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _ACCEL_
#define _ACCEL_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// Accel(double x, Matrix Y)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the time derivative of the satellite state vector, including velocity and 
 *         acceleration due to Earth's harmonic gravity field and third-body gravitational 
 *         perturbations from the Sun, Moon, and planets.
 *
 *  This function accounts for Earth's precession, nutation, rotation, and polar motion using
 *  IAU/IERS models and DE430 planetary ephemerides.
 *
 *  @param [in] x  Time since the reference epoch (in seconds).
 *  @param [in] Y  State vector of the satellite [position; velocity] in the inertial frame (ICRF/EME2000).
 *
 *  @return Matrix dY  Time derivative of the state vector [velocity; acceleration] in the inertial frame.
 */
//-----------------------------------------------------------------------------------------------
Matrix Accel(double x, Matrix Y);

#endif