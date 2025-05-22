#include "..\include\AccelPointMass.hpp"

Matrix AccelPointMass(Matrix r, Matrix s, double GM){
	Matrix d = r - s;

Matrix a =  ( d/pow((norm(d)),3) + s/pow(norm(s),3))*-GM;
	return a;
}