#include "..\include\EccAnom.hpp"


 double EccAnom(double M, double e){

int maxit = 15;
int i = 1;
double E,f;
double eps=1e-10;

// Starting value
M = fmod(M, 2.0*M_PI);

if (e<0.8)
    E = M; 
else
    E = M_PI;

f = E - e*sin(E) - M;
E = E - f / ( 1.0 - e*cos(E) );

// Iteration
while (fabs(f) > 100*eps){  
    f = E - e*sin(E) - M;
    E = E - f / ( 1.0 - e*cos(E) );
    i = i+1;
    if (i==maxit){
        cout<<(" convergence problems in EccAnom");
		break;
    }  
}
return E;
 }

