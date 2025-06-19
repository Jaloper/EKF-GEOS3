// $Header$
//--------------------------------------------------------------------------------
// PrecMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/07
//
/** @file PrecMatrix.hpp
 *  @brief This header file declares a function to compute the precession transformation matrix.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _PRECMATRIX_
#define _PRECMATRIX_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// PrecMatrix(double Mjd_1, double Mjd_2)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the precession transformation matrix between two Modified Julian Dates.
 *
 *  @param [in] Mjd_1 Initial Modified Julian Date (reference epoch).
 *  @param [in] Mjd_2 Final Modified Julian Date (target epoch).
 *  @return Matrix A 3x3 transformation matrix representing the precession from the reference epoch to the target epoch.
 */
Matrix PrecMatrix(double Mjd_1, double Mjd_2);

#endif