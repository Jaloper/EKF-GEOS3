// $Header$
//--------------------------------------------------------------------------------
// GHAMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/16
//
/** @file GHAMatrix.hpp
 *  @brief This header file declares the function that computes the Greenwich Hour Angle rotation matrix.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _GHAMATRIX_
#define _GHAMATRIX_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// GHAMatrix(double Mjd_UT1)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the rotation matrix for the Greenwich Hour Angle (GHA).
 *
 *  This function returns a rotation matrix about the Z-axis corresponding to
 *  the Greenwich Apparent Sidereal Time at the given UT1 Modified Julian Date.
 *
 *  @param [in] Mjd_UT1  Modified Julian Date in UT1 time scale.
 *  @return Matrix       3x3 rotation matrix representing Earth's rotation due to GHA.
 */
Matrix GHAMatrix(double Mjd_UT1);

#endif
