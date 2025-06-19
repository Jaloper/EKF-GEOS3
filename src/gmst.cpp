// $Source$
//--------------------------------------------------------------------------------
// gmst
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/07
//
/** @file gmst.cpp
 *  @brief Implementation of Greenwich Mean Sidereal Time computation.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\gmst.hpp"
#include "..\include\Frac.hpp"

double gmst(double Mjd_UT1) {
    const double Secs = 86400.0;                       // Seconds per day
    const double MJD_J2000 = 51544.5;

    double Mjd_0 = floor(Mjd_UT1);
    double UT1 = Secs * (Mjd_UT1 - Mjd_0);             // [s]
    double T_0 = (Mjd_0 - MJD_J2000) / 36525.0;
    double T = (Mjd_UT1 - MJD_J2000) / 36525.0;

    double gmst = 24110.54841 + 8640184.812866 * T_0 + 1.002737909350795 * UT1
                  + (0.093104 - 6.2e-6 * T) * T * T;  // [s]

    double gmstime = 2 * M_PI * Frac(gmst / Secs);       // [rad], 0..2pi
    return gmstime;
}
