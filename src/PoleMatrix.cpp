// $Source$
//--------------------------------------------------------------------------------
// PoleMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file PoleMatrix.cpp
 *  @brief Implements the computation of the polar motion transformation matrix.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\PoleMatrix.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_x.hpp"

Matrix PoleMatrix (double xp,double yp) {
	Matrix Ry = R_y(-xp);
	Matrix Rx = R_x(-yp);
	return Ry * Rx;
}
