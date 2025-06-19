// $Header$
//--------------------------------------------------------------------------------
// SAT_Const
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/18
//
/** @file SAT_Const.hpp
 *  @brief This header file defines a class containing constants for satellite orbit computations.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#ifndef M_PI
/** @brief Defines the mathematical constant pi if not already defined. */
#define M_PI 3.14159265358979323846
#endif

#ifndef SAT_CONST_HPP
#define SAT_CONST_HPP

#include <cmath>

/** @class SAT_Const
 *  @brief Class containing static constants for satellite orbit propagation and related computations.
 */
class SAT_Const {
public:
    // Mathematical constants
    /** @brief Mathematical constant pi. */
    static constexpr double pi = M_PI;

    /** @brief Two times pi (2π). */
    static constexpr double pi2 = 2.0 * M_PI;

    /** @brief Conversion factor from degrees to radians (π/180). */
    static constexpr double Rad = M_PI / 180.0;

    /** @brief Conversion factor from radians to degrees (180/π). */
    static constexpr double Deg = 180.0 / M_PI;

    /** @brief Conversion factor from radians to arcseconds (3600 * 180/π). */
    static constexpr double Arcs = 3600.0 * 180.0 / M_PI;

    // General
    /** @brief Modified Julian Date of the J2000 epoch (January 1, 2000, 12:00 TT). */
    static constexpr double MJD_J2000 = 51544.5;

    /** @brief Time difference between epoch B1950 and J2000 in Julian centuries. */
    static constexpr double T_B1950 = -0.500002108;

    /** @brief Speed of light in vacuum [m/s], as per DE430. */
    static constexpr double c_light = 299792458.0;

    /** @brief Astronomical unit [m], as per DE430. */
    static constexpr double AU = 149597870700.0;

    // Physical parameters
    /** @brief Earth's equatorial radius [m], as per DE430. */
    static constexpr double R_Earth = 6378136.3;

    /** @brief Earth's flattening factor, as per WGS-84. */
    static constexpr double f_Earth = 1.0 / 298.257223563;

    /** @brief Sun's radius [m], as per DE430. */
    static constexpr double R_Sun = 696000e3;

    /** @brief Moon's radius [m], as per DE430. */
    static constexpr double R_Moon = 1738e3;

    // Earth rotation
    /** @brief Earth's angular velocity [rad/s], as per WGS-84 (15.04106717866910 deg/hour converted to rad/s). */
    static constexpr double omega_Earth = 15.04106717866910 / 3600.0 * Rad;

    // Gravitational coefficients [m^3/s^2], as per DE430
    /** @brief Gravitational constant times Earth's mass [m^3/s^2]. */
    static constexpr double GM_Earth = 398600.435436e9;

    /** @brief Gravitational constant times Sun's mass [m^3/s^2]. */
    static constexpr double GM_Sun = 132712440041.939400e9;

    /** @brief Gravitational constant times Moon's mass [m^3/s^2], derived from Earth-Moon mass ratio. */
    static constexpr double GM_Moon = GM_Earth / 81.30056907419062;

    /** @brief Gravitational constant times Mercury's mass [m^3/s^2]. */
    static constexpr double GM_Mercury = 22031.780000e9;

    /** @brief Gravitational constant times Venus's mass [m^3/s^2]. */
    static constexpr double GM_Venus = 324858.592000e9;

    /** @brief Gravitational constant times Mars's mass [m^3/s^2]. */
    static constexpr double GM_Mars = 42828.375214e9;

    /** @brief Gravitational constant times Jupiter's mass [m^3/s^2]. */
    static constexpr double GM_Jupiter = 126712764.800000e9;

    /** @brief Gravitational constant times Saturn's mass [m^3/s^2]. */
    static constexpr double GM_Saturn = 37940585.200000e9;

    /** @brief Gravitational constant times Uranus's mass [m^3/s^2]. */
    static constexpr double GM_Uranus = 5794548.600000e9;

    /** @brief Gravitational constant times Neptune's mass [m^3/s^2]. */
    static constexpr double GM_Neptune = 6836527.100580e9;

    /** @brief Gravitational constant times Pluto's mass [m^3/s^2]. */
    static constexpr double GM_Pluto = 977.0000000000009e9;

    // Solar radiation pressure
    /** @brief Solar radiation pressure at 1 AU [N/m^2], calculated as solar constant (~1367 W/m^2) divided by speed of light, per IERS 96. */
    static constexpr double P_Sol = 1367.0 / c_light;
};

#endif