// $Source$
//--------------------------------------------------------------------------------
// EqnEquinox
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file EqnEquinox.cpp
 *  @brief Calculates the equation of the equinoxes based on the nutation in longitude and the mean obliquity.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\EqnEquinox.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\MeanObliquity.hpp"


double EqnEquinox(double Mjd_TT) {
	
    // Nutation in longitude and obliquity
    auto [dpsi, deps] = NutAngles(Mjd_TT);
	
    // Equation of the equinoxes
    return dpsi * cos(MeanObliquity(Mjd_TT));
}
