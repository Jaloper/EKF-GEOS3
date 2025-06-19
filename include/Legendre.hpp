// $Header$
//--------------------------------------------------------------------------------
// Legendre
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/29
//
/** @file Legendre.hpp
 *  @brief This header file declares a function to compute Legendre polynomials and their derivatives.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _LEGENDRE_
#define _LEGENDRE_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// Legendre(int n, int m, double fi)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the associated Legendre polynomials and their derivatives up to degree n and order m for a given colatitude.
 *
 *  @param [in] n  Maximum degree of the Legendre polynomials.
 *  @param [in] m  Maximum order of the Legendre polynomials (m <= n).
 *  @param [in] fi Colatitude angle in radians.
 *  @return std::tuple<Matrix&, Matrix&> Tuple containing two matrices: pnm (associated Legendre polynomials) and dpnm (their derivatives), both of size (n+1) x (m+1).
 */
std::tuple<Matrix&, Matrix&> Legendre(int n, int m, double fi);

#endif