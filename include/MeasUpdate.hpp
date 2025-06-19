// $Header$
//--------------------------------------------------------------------------------
// MeasUpdate
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/16
//
/** @file MeasUpdate.hpp
 *  @brief This header file declares a function to perform a measurement update in a Kalman filter.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// MeasUpdate(Matrix x, Matrix z, Matrix g, Matrix s, Matrix G, Matrix P, int n)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Performs a measurement update step in a Kalman filter, computing the Kalman gain, updated state estimate, and updated covariance matrix.
 *
 *  @param [in] x  State estimate vector (n x 1).
 *  @param [in] z  Measurement vector (m x 1 or 1 x m).
 *  @param [in] g  Predicted measurement vector (m x 1).
 *  @param [in] s  Measurement noise standard deviation vector (m x 1).
 *  @param [in] G  Measurement matrix (m x n).
 *  @param [in] P  State covariance matrix (n x n).
 *  @param [in] n  Dimension of the state vector.
 *  @return std::tuple<Matrix&, Matrix&, Matrix&> Tuple containing the Kalman gain matrix (n x m), updated state estimate (n x 1), and updated covariance matrix (n x n).
 */
std::tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix x, Matrix z, Matrix g, Matrix s, Matrix G, Matrix P, int n);

#endif