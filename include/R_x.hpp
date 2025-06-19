// $Header$
//--------------------------------------------------------------------------------
// R_x
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/09
//
/** @file R_x.hpp
 *  @brief This header file declares a function to compute the rotation matrix around the x-axis.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _R_X_
#define _R_X_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// R_x(double angle)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the 3x3 rotation matrix for a rotation around the x-axis by a given angle.
 *
 *  @param [in] angle Rotation angle in radians.
 *  @return Matrix A 3x3 rotation matrix representing the rotation around the x-axis.
 */
Matrix R_x(double angle);

#endif