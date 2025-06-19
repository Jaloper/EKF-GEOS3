// $Source$
//--------------------------------------------------------------------------------
// Frac
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/17
//
/** @file Frac.cpp
 *  @brief Calculates the fractional part of a number.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\Frac.hpp"

double Frac(double x) {
    return x - floor(x);
}
