// $Source$
//--------------------------------------------------------------------------------
// sign_
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/18
//
/** @file sign_.cpp
 *  @brief Implements the computation of the sign of a number based on another number's sign.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\sign_.hpp"

double sign_(double a, double b) {
    if (b >= 0.0)
        return fabs(a);
    else
        return -fabs(a);
}
