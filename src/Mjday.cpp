// $Source$
//--------------------------------------------------------------------------------
// Mjday
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/17
//
/** @file Mjday.cpp
 *  @brief Implements the computation of the Modified Julian Date from calendar date and time.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\Mjday.hpp"

double Mjday(int yr, int mon, int day, int hr = 0, int min = 0, double sec = 0.0) {
    double jd = 367.0 * yr
        - floor((7 * (yr + floor((mon + 9) / 12.0))) * 0.25)
        + floor(275 * mon / 9.0)
        + day + 1721013.5
        + (((sec / 60.0 + min) / 60.0 + hr) / 24.0);

    return jd - 2400000.5;
}
