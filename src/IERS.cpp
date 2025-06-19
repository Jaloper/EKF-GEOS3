// $Source$
//--------------------------------------------------------------------------------
// IERS
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/25
//
/** @file IERS.cpp
 *  @brief Implementation of Earth orientation parameter retrieval from IERS data.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\IERS.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\global.hpp"
#include <stdexcept>
#include <tuple>

std::tuple<double, double, double, double, double, double, double, double, double> 
IERS(
    double Mjd_UTC, char interp
) {

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;

    if (interp == 'l') {
        // linear interpolation
        double mjd = (std::floor(Mjd_UTC));
        int i = 0;
        for (int j = 1; j <= eopdata.n_column; j++) {
            if (fabs(eopdata(4, j) - mjd)<0.1) {
                i = j;
                break;
            }
        }
        if (i == 0) {
            printf("IERS: Not find\n");
			exit(EXIT_FAILURE);
        }
		Matrix preeop= extract_column(eopdata,i);
		preeop=transpose(preeop);
		Matrix nexteop= extract_column(eopdata,i+1);
		nexteop=transpose(nexteop);
        double mfme = 1440.0 * (Mjd_UTC - floor(Mjd_UTC));
        double fixf = mfme / 1440.0;
		
            x_pole  = preeop(5)+(nexteop(5)-preeop(5))*fixf;
		y_pole  = preeop(6)+(nexteop(6)-preeop(6))*fixf;
		UT1_UTC = preeop(7)+(nexteop(7)-preeop(7))*fixf;
		LOD     = preeop(8)+(nexteop(8)-preeop(8))*fixf;
		dpsi    = preeop(9)+(nexteop(9)-preeop(9))*fixf;
		deps    = preeop(10)+(nexteop(10)-preeop(10))*fixf;
		dx_pole = preeop(11)+(nexteop(11)-preeop(11))*fixf;
		dy_pole = preeop(12)+(nexteop(12)-preeop(12))*fixf;
		TAI_UTC = preeop(13);
		
        x_pole  /= SAT_Const::Arcs;
        y_pole  /= SAT_Const::Arcs;
        dpsi    /= SAT_Const::Arcs;
        deps    /= SAT_Const::Arcs;
        dx_pole /= SAT_Const::Arcs;
        dy_pole /= SAT_Const::Arcs;
    }
    else if (interp == 'n') {
        int mjd = static_cast<int>(floor(Mjd_UTC));
        int i = -1;

        for (int col = 1; col <= eopdata.n_column; ++col) {
            if (static_cast<int>(eopdata(4, col)) == mjd) {
                i = col;
                break;
            }
        }

        if (i == -1) {
            throw std::runtime_error("MJD not found");
        }

        x_pole  = eopdata(5,i) / SAT_Const::Arcs;
        y_pole  = eopdata(6,i) / SAT_Const::Arcs;
        UT1_UTC = eopdata(7,i);
        LOD     = eopdata(8,i);
        dpsi    = eopdata(9,i) / SAT_Const::Arcs;
        deps    = eopdata(10,i) / SAT_Const::Arcs;
        dx_pole = eopdata(11,i) / SAT_Const::Arcs;
        dy_pole = eopdata(12,i) / SAT_Const::Arcs;
        TAI_UTC = eopdata(13,i);
    }
	
    return std::make_tuple(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}
