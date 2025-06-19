// $Header$
//--------------------------------------------------------------------------------
// global
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/25
//
/** @file global.hpp
 *  @brief This header file declares global variables and functions for handling Earth orientation parameters, gravity models, and observation data.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef _GLOBAL_
#define _GLOBAL_

#include "..\include\matrix.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <string.h>

//-----------------------------------------------------------------------------------------------
// Param structure
//-----------------------------------------------------------------------------------------------
/**
 *  @struct Param
 *  @brief Structure to hold auxiliary parameters for computations.
 *
 *  @var Mjd_UTC Modified Julian Date in Coordinated Universal Time (UTC).
 *  @var Mjd_TT  Modified Julian Date in Terrestrial Time (TT).
 *  @var n       Degree of the gravity model.
 *  @var m       Order of the gravity model.
 *  @var sun     Flag to include Sun perturbations (1 = include, 0 = exclude).
 *  @var moon    Flag to include Moon perturbations (1 = include, 0 = exclude).
 *  @var planets Flag to include planetary perturbations (1 = include, 0 = exclude).
 */
typedef struct {
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
} Param;

extern Param AuxParam;       ///< Global auxiliary parameters.
extern Matrix eopdata;       ///< Earth orientation parameters data matrix.
extern Matrix Cnm;           ///< Cosine coefficients for gravity model.
extern Matrix Snm;           ///< Sine coefficients for gravity model.
extern Matrix PC;            ///< DE430 planetary ephemeris coefficients.
extern Matrix obs;           ///< Observation data matrix.
extern int n_eqn;            ///< Number of equations.

//-----------------------------------------------------------------------------------------------
// eop19620101(int c)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Loads Earth orientation parameters from the file eop19620101.txt.
 *
 *  @param [in] c Number of columns (data points) to read.
 */
void eop19620101(int c);

//-----------------------------------------------------------------------------------------------
// GGM03S(int n)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Loads the GGM03S gravity model coefficients from the file GGM03S.txt.
 *
 *  @param [in] n Maximum degree and order of the gravity model.
 */
void GGM03S(int n);

//-----------------------------------------------------------------------------------------------
// DE430Coeff(int i, int j)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Loads DE430 planetary ephemeris coefficients from the file DE430Coeff.txt.
 *
 *  @param [in] i Number of rows to read.
 *  @param [in] j Number of columns to read.
 */
void DE430Coeff(int i, int j);

//-----------------------------------------------------------------------------------------------
// GEOS3(int n)
//-----------------------------------------------------------------------------------------------
/**
 *  @brief Loads GEOS-3 observation data from the file GEOS3.txt.
 *
 *  @param [in] n Number of observation data points to read.
 */
void GEOS3(int n);

#endif