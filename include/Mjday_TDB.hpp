// $Header$
//--------------------------------------------------------------------------------
// Mjday_TDB
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/18
//
/** @file Mjday_TDB.hpp
 *  @brief This header file declares a function to convert Modified Julian Date in Terrestrial Time to Barycentric Dynamical Time.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _MJDAY_TDB_
#define _MJDAY_TDB_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// Mjday_TDB(double Mjd_TT)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Converts a Modified Julian Date in Terrestrial Time (TT) to Modified Julian Date in Barycentric Dynamical Time (TDB).
 *
 *  @param [in] Mjd_TT Modified Julian Date in Terrestrial Time (TT).
 *  @return double The Modified Julian Date in Barycentric Dynamical Time (TDB).
 */
double Mjday_TDB(double Mjd_TT);

#endif