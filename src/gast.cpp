// $Source$
//--------------------------------------------------------------------------------
// gast
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/16
//
/** @file gast.cpp
 *  @brief Computes the Greenwich Apparent Sidereal Time (GAST).
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\gast.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\gmst.hpp"

double gast (double Mjd_UT1){
    double gast_raw = gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1);
    double gast_mod = fmod(gast_raw, 2 * M_PI);
    if (gast_mod < 0)
        gast_mod += 2 * M_PI;
    return gast_mod;
}
