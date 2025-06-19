// $Source$
//--------------------------------------------------------------------------------
// AccelPointMass
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/09
//
/** @file AccelPointMass.cpp
 *  @brief Computes the acceleration on a satellite due to the gravitational attraction
 *         of an external point mass (e.g., Sun, Moon, or planet).
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\AccelPointMass.hpp"

Matrix AccelPointMass(Matrix r, Matrix s, double GM){
	Matrix d = r - s;

Matrix a =  ( d/pow((norm(d)),3) + s/pow(norm(s),3))*-GM;
	return a;
}
