#include "../include/matrix.hpp"
#include <iostream>

int main(){
	Matrix M1(3,2);
	M1(1,1) = 5.0;
	
	Matrix M2(3,2);
	M2(1,1) = 3.0;
	
	Matrix M3=M1+M2;
	
	cout<<M1<<end;
	cout<<M2<<end;
	cout<<M3<<end;
	
	return 0;
}