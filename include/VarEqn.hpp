// $Header$
//--------------------------------------------------------------------------------
// VarEqn
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/18
//
/** @file VarEqn.hpp
 *  @brief This header file declares a function to compute the variational equations for orbit propagation.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _VAREQN_
#define _VAREQN_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// VarEqn(double x, Matrix yPhi)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the time derivative of the state vector and state transition matrix for variational equations.
 *
 *  @param [in] x    Time offset from the reference epoch in seconds.
 *  @param [in] yPhi Combined state vector and state transition matrix (vectorized, 42x1).
 *  @return Matrix Time derivative of the combined state vector and state transition matrix (42x1).
 */
Matrix VarEqn(double x, Matrix yPhi);

#endif