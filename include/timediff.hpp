// $Header$
//--------------------------------------------------------------------------------
// timediff
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/18
//
/** @file timediff.hpp
 *  @brief This header file declares a function to compute time differences between various time scales.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _TIMEDIFF_
#define _TIMEDIFF_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// timediff(double UT1_UTC, double TAI_UTC)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Computes time differences between various time scales (UT1, UTC, TAI, GPS, TT).
 *
 *  @param [in] UT1_UTC Difference between UT1 and UTC in seconds.
 *  @param [in] TAI_UTC Difference between TAI and UTC in seconds.
 *  @return std::tuple<double, double, double, double, double> Tuple containing time differences:
 *         - UT1_TAI: UT1 - TAI [s]
 *         - UTC_GPS: UTC - GPS [s]
 *         - UT1_GPS: UT1 - GPS [s]
 *         - TT_UTC: TT - UTC [s]
 *         - GPS_UTC: GPS - UTC [s]
 */
std::tuple<double, double, double, double, double> timediff(double UT1_UTC, double TAI_UTC);

#endif