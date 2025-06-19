// $Header$
//--------------------------------------------------------------------------------
// TimeUpdate
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/30
//
/** @file TimeUpdate.hpp
 *  @brief This header file declares functions to perform the time update step in a Kalman filter.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// TimeUpdate(Matrix P, Matrix Phi)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Performs the time update step in a Kalman filter without process noise.
 *
 *  @param [in] P   State covariance matrix.
 *  @param [in] Phi State transition matrix.
 *  @return Matrix Updated state covariance matrix.
 */
Matrix TimeUpdate(Matrix P, Matrix Phi);

//-----------------------------------------------------------------------------------------------
// TimeUpdate(Matrix P, Matrix Phi, Matrix Qdt)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Performs the time update step in a Kalman filter with process noise.
 *
 *  @param [in] P   State covariance matrix.
 *  @param [in] Phi State transition matrix.
 *  @param [in] Qdt Process noise covariance matrix.
 *  @return Matrix Updated state covariance matrix.
 */
Matrix TimeUpdate(Matrix P, Matrix Phi, Matrix Qdt);

#endif