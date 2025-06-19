// $Source$
//--------------------------------------------------------------------------------
// AzElPa
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/18
//
/** @file AzElPa.cpp
 *  @brief Calculates azimuth, elevation, 
 *         and their partial derivatives with respect to the position vector.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include <tuple>
#include "..\include\AzElPa.hpp"

std::tuple<double, double, Matrix&, Matrix&> AzElPa(Matrix s) {
    const double pi2 = 2.0 * M_PI;

    double rho = std::sqrt(s(1) * s(1) + s(2) * s(2));

    // Angles
    double Az = std::atan2(s(1), s(2));
    if (Az < 0.0) {
        Az += pi2;
    }

    double El = std::atan(s(3) / rho);

    // Partials
    double rho2 = rho * rho;
    double s_norm2 = dot(s, s);

    Matrix &dAds = zeros(1, 3);
    Matrix &dEds = zeros(1, 3);

    dAds(1, 1) = s(2) / rho2;
    dAds(1, 2) = -s(1) / rho2;
    dAds(1, 3) = 0.0;

    dEds(1, 1) = -s(1) * s(3) / rho;
    dEds(1, 2) = -s(2) * s(3) / rho;
    dEds(1, 3) = rho;

    dEds = dEds / s_norm2;

    return std::tie(Az, El, dAds, dEds);
}
