// $Header$
//--------------------------------------------------------------------------------
// DEInteg
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/20
//
/** @file DEInteg.hpp
 *  @brief This header file declares the function to solve a system of differential equations using an adaptive integration method with error control.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _DEINTEG_
#define _DEINTEG_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// DEInteg(Matrix func(double t, Matrix y), double t, double tout, double relerr, double abserr, 
//         int n_eqn, Matrix& y)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Solves a system of differential equations using an adaptive integration method.
 *  
 *  @param [in] func  Function that computes the derivatives of the system at time `t` and state vector `y`.
 *  @param [in] t     Initial time of the integration.
 *  @param [in] tout  Final time at which to obtain the solution.
 *  @param [in] relerr Maximum relative error tolerance.
 *  @param [in] abserr Maximum absolute error tolerance.
 *  @param [in] n_eqn Number of equations in the system.
 *  @param [in,out] y Initial conditions of the system of equations.
 *  
 *  @return Matrix    Solution of the differential system at time `tout`.
 */
Matrix DEInteg(Matrix func(double t, Matrix y), double t, double tout, double relerr, double abserr, 
               int n_eqn, Matrix& y);

#endif
