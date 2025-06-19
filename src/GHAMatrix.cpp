// $Source$
//--------------------------------------------------------------------------------
// GHAMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/16
//
/** @file GHAMatrix.cpp
 *  @brief Implementation of Greenwich Hour Angle rotation matrix computation.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\GHAMatrix.hpp"
#include "..\include\R_z.hpp"
#include "..\include\gast.hpp"

 Matrix GHAMatrix(double Mjd_UT1){
	return R_z( gast(Mjd_UT1) );
 }
 