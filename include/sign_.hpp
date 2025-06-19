// $Header$
//--------------------------------------------------------------------------------
// sign_
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/18
//
/** @file sign_.hpp
 *  @brief This header file declares a function to compute the sign of a number based on another number's sign.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _SIGN__
#define _SIGN__

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// sign_(double a, double b)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Returns the absolute value of 'a' with the sign of 'b'.
 *
 *  @param [in] a The number whose absolute value is used.
 *  @param [in] b The number whose sign determines the result's sign.
 *  @return double The absolute value of 'a' with the sign of 'b' (positive if b >= 0, negative otherwise).
 */
double sign_(double a, double b);

#endif