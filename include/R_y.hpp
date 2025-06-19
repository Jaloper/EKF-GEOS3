// $Header$
//--------------------------------------------------------------------------------
// R_y
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/09
//
/** @file R_y.hpp
 *  @brief This header file declares a function to compute the rotation matrix around the y-axis.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _R_Y_
#define _R_Y_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// R_y(double angle)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the 3x3 rotation matrix for a rotation around the y-axis by a given angle.
 *
 *  @param [in] angle Rotation angle in radians.
 *  @return Matrix A 3x3 rotation matrix representing the rotation around the y-axis.
 */
Matrix R_y(double angle);

#endif