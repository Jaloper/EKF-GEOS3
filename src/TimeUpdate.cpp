// $Source$
//--------------------------------------------------------------------------------
// TimeUpdate
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/30
//
/** @file TimeUpdate.cpp
 *  @brief Implements the time update step for a Kalman filter with and without process noise.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\TimeUpdate.hpp"


Matrix TimeUpdate(Matrix P, Matrix Phi, Matrix Qdt) {
       Matrix temp = Phi * P;
    Matrix result = temp * transpose(Phi);
    result = result + Qdt;
    P = result;
    return P;
}

Matrix TimeUpdate(Matrix P, Matrix Phi) {
        Matrix temp = Phi * P;
    Matrix result = temp * transpose(Phi);
    return result;
}
