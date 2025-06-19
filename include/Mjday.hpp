// $Header$
//--------------------------------------------------------------------------------
// Mjday
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/17
//
/** @file Mjday.hpp
 *  @brief This header file declares a function to compute the Modified Julian Date from calendar date and time.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _MJDAY_
#define _MJDAY_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// Mjday(int yr, int mon, int day, int hr, int min, double sec)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes the Modified Julian Date (MJD) from a given calendar date and time.
 *
 *  @param [in] yr   Year (e.g., 2025).
 *  @param [in] mon  Month (1 to 12).
 *  @param [in] day  Day of the month (1 to 31).
 *  @param [in] hr   Hour (0 to 23, default is 0).
 *  @param [in] min  Minute (0 to 59, default is 0).
 *  @param [in] sec  Second (0.0 to 59.999..., default is 0.0).
 *  @return double   The Modified Julian Date corresponding to the input date and time.
 */
double Mjday(int yr, int mon, int day, int hr, int min, double sec);

#endif