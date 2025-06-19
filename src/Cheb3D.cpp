// $Source$
//--------------------------------------------------------------------------------
// Cheb3D
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/04/11
//
/** @file Cheb3D.cpp
 *  @brief Calculates the 3D Chebyshev interpolation based on the input coefficients.
 *
 *  @author Jarein López Ruiz
 *  @bug No known bugs.
 */
//--------------------------------------------------------------------------------

#include "..\include\Cheb3D.hpp"

 Matrix Cheb3D(double t, int N, double Ta,double Tb, Matrix Cx, Matrix Cy, Matrix Cz){

if ( (t<Ta) || (Tb<t) )
    cout<<"ERROR: Time out of range in Cheb3D::Value\n";

double tau = (2*t-Ta-Tb)/(Tb-Ta);  

Matrix f1 = zeros(1,3);
Matrix f2 = zeros(1,3);
Matrix aux = zeros(3);
Matrix old_f1= zeros(1,3);
for (int i = N; i >= 2; i--) {
    old_f1 = f1;
	aux(1)=Cx(i);
	aux(2)=Cy(i);
	aux(3)=Cz(i);
    f1 = f1*(2*tau)-f2+aux;
    f2 = old_f1;
}
	aux(1)=Cx(1);
	aux(2)=Cy(1);
	aux(3)=Cz(1);
return f1*tau-f2+aux;
 }
 