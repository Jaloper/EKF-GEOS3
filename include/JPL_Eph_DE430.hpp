// $Header$
//--------------------------------------------------------------------------------
// JPL_Eph_DE430
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/09
//
/** @file JPL_Eph_DE430.hpp
 *  @brief This header file declares a function to compute planetary positions using JPL DE430 ephemerides.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _JPL_EPH_DE430_
#define _JPL_EPH_DE430_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// JPL_Eph_DE430(double Mjd_TDB)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the positions of Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, Moon, and Sun relative to Earth using JPL DE430 ephemerides.
 *
 *  @param [in] Mjd_TDB  Modified Julian Date in Barycentric Dynamical Time (TDB).
 *  @return std::tuple<Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&>  Tuple containing position vectors (in meters) for Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, Moon, and Sun, each as a transposed 3x1 Matrix.
 */
std::tuple<Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&> JPL_Eph_DE430(double Mjd_TDB);

#endif