// $Source$
//--------------------------------------------------------------------------------
// MeanObliquity
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/17
//
/** @file MeanObliquity.cpp
 *  @brief Implements the computation of the mean obliquity of the ecliptic.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\MeanObliquity.hpp"


double MeanObliquity(double Mjd_TT) {
    const double MJD_J2000 = 51544.5;
    const double Rad = M_PI / 180.0;

    double T = (Mjd_TT - MJD_J2000) / 36525.0;

    double MOblq = Rad * (84381.448 / 3600.0 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0);

    return MOblq;
}
