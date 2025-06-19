// $Header$
//--------------------------------------------------------------------------------
// PoleMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file PoleMatrix.hpp
 *  @brief This header file declares a function to compute the polar motion transformation matrix.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _POLEMATRIX_
#define _POLEMATRIX_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// PoleMatrix(double xp, double yp)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the polar motion transformation matrix based on pole coordinates.
 *
 *  @param [in] xp Polar motion coordinate along the x-axis in radians.
 *  @param [in] yp Polar motion coordinate along the y-axis in radians.
 *  @return Matrix A 3x3 transformation matrix representing the polar motion correction.
 */
Matrix PoleMatrix(double xp, double yp);

#endif