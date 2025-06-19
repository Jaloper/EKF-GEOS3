// $Header$
//--------------------------------------------------------------------------------
// NutMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file NutMatrix.hpp
 *  @brief This header file declares a function to compute the nutation transformation matrix.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _NUTMATRIX_
#define _NUTMATRIX_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// NutMatrix(double Mjd_TT)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the nutation transformation matrix from mean to true equator and equinox for a given Modified Julian Date in Terrestrial Time.
 *
 *  @param [in] Mjd_TT Modified Julian Date in Terrestrial Time (TT).
 *  @return Matrix A 3x3 transformation matrix representing the nutation from mean to true equator and equinox.
 */
Matrix NutMatrix(double Mjd_TT);

#endif