// $Header$
//--------------------------------------------------------------------------------
// IERS
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/25
//
/** @file IERS.hpp
 *  @brief This header file declares the function that retrieves Earth orientation parameters (EOP) from IERS data.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _IERS_
#define _IERS_

#include "..\include\matrix.hpp"
#include <cmath>
#include <tuple>

//-----------------------------------------------------------------------------------------------
// IERS(double Mjd_UTC, char interp='n')
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Retrieves Earth orientation parameters for a given UTC date.
 *
 *  This function extracts or interpolates Earth Orientation Parameters (EOP) from global data,
 *  provided by the IERS, corresponding to the given UTC Modified Julian Date.
 *
 *  @param [in] Mjd_UTC  Modified Julian Date in UTC time scale.
 *  @param [in] interp   Interpolation flag: 'n' for nearest neighbor, 'l' for linear interpolation.
 *  @return std::tuple<double, double, double, double, double, double, double, double, double>
 *          Tuple containing:
 *            - x_pole  [rad]  Polar motion x-component
 *            - y_pole  [rad]  Polar motion y-component
 *            - UT1_UTC [s]    Difference UT1 - UTC
 *            - LOD     [s]    Excess length of day
 *            - dpsi    [rad]  Nutation correction in longitude
 *            - deps    [rad]  Nutation correction in obliquity
 *            - dx_pole [rad]  Polar motion correction x-component
 *            - dy_pole [rad]  Polar motion correction y-component
 *            - TAI_UTC [s]    Difference TAI - UTC
 */
std::tuple<double, double, double, double, double, double, double, double, double> 
IERS(
    double Mjd_UTC, char interp = 'n'
);

#endif
