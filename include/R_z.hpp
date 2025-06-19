// $Header$
//--------------------------------------------------------------------------------
// R_z
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/09
//
/** @file R_z.hpp
 *  @brief This header file declares a function to compute the rotation matrix around the z-axis.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _R_Z_
#define _R_Z_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// R_z(double angle)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the 3x3 rotation matrix for a rotation around the z-axis by a given angle.
 *
 *  @param [in] angle Rotation angle in radians.
 *  @return Matrix A 3x3 rotation matrix representing the rotation around the z-axis.
 */
Matrix R_z(double angle);

#endif