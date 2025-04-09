#include "..\include\matrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\AccelPointMass.hpp"
#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_mult_01() {
	
    Matrix A(2, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
    
    Matrix B(3, 2);
    B(1,1) = 7;  B(1,2) = 8;
    B(2,1) = 9;  B(2,2) = 10;
    B(3,1) = 11; B(3,2) = 12;
    
    Matrix C(2, 2);
    C(1,1) = 58;  C(1,2) = 64;
    C(2,1) = 139; C(2,2) = 154;
    
    Matrix R = A * B;
    
    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_div_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;
    
    Matrix B(2, 2);
    B(1,1) = 5; B(1,2) = 6;
    B(2,1) = 7; B(2,2) = 8;
    
    Matrix C(2, 2);
    C(1,1) = 3.0;  C(1,2) = -2.0;
    C(2,1) = 2.0;  C(2,2) = -1.0;
    
    Matrix R = A / B;
    
    _assert(m_equals(C, R, 1e-8));
    return 0;
}

int m_assign_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;
    
    Matrix B=A;
    
    _assert(m_equals(A, B, 1e-8));
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_inv_01() {
    Matrix A(2, 2);
    A(1,1) = 4; A(1,2) = 7;
    A(2,1) = 2; A(2,2) = 6;
    
    Matrix B(2, 2);
    B(1,1) = 0.6;  B(1,2) = -0.7;
    B(2,1) = -0.2; B(2,2) = 0.4;
    
    Matrix R = inv(A);
    
    _assert(m_equals(B, R, 1e-10));
    return 0;
}

int m_eye_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 0;
    A(2,1) = 0; A(2,2) = 1;
    
    Matrix B=eye(2);
    
    _assert(m_equals(B, A, 1e-10));
    return 0;
}

int m_transpose_01() {
    Matrix A(2, 2);
    A(1,1) = 4; A(1,2) = 7;
    A(2,1) = 2; A(2,2) = 6;
    
    Matrix B(2, 2);
    B(1,1) = 4; B(1,2) = 2;
    B(2,1) = 7; B(2,2) = 6;
    
	Matrix R = transpose(A);
	
    _assert(m_equals(B, R, 1e-10));
    return 0;
}

int op_sum_01() {
    double d = 5.2;
	
	Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 0;
    A(2,1) = 3; A(2,2) = 4;
	
	Matrix B(2, 2);
    B(1,1) = 6.2; B(1,2) = 5.2;
    B(2,1) = 8.2; B(2,2) = 9.2;
	
	
	Matrix R = A + d;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int op_sub_01() {
    double d = -5.2;
	
	Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 0;
    A(2,1) = 3; A(2,2) = 4;
	
	Matrix B(2, 2);
    B(1,1) = 6.2; B(1,2) = 5.2;
    B(2,1) = 8.2; B(2,2) = 9.2;
	
	
	Matrix R = A - d;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int op_mult_01() {
    double d = 2;
	
	Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 0;
    A(2,1) = 3; A(2,2) = 4;
	
	Matrix B(2, 2);
    B(1,1) = 2; B(1,2) = 0;
    B(2,1) = 6; B(2,2) = 8;
	
	
	Matrix R = A * d;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int op_div_01() {
    double d = 2;
	
	Matrix A(2);
    A(1,1) = 1; A(1,2) = 0;
    A(2,1) = 3; A(2,2) = 4;
	
	Matrix B(2, 2);
    B(1,1) = 0.5; B(1,2) = 0;
    B(2,1) = 1.5; B(2,2) = 2;
	
	
	Matrix R = A / d;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_norm_01() {
	
	Matrix A(2);
    A(1,1) = 3; A(1,2) = 0;
    A(2,1) = 0; A(2,2) = 4;
	
	double result=5;
	double R = norm(A);
    _assert(m_equals(result, R, 1e-10));
    
    return 0;
}

int m_dot_01() {
	
	Matrix A(1,2);
    A(1,1) = 4; A(1,2) = 7;
	
	double result=65;
	double R = dot(A,A);
    _assert(m_equals(result, R, 1e-10));
    
    return 0;
}

int m_cross_01() {
	
	Matrix A(1,3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	
	Matrix B(1,3);
	B(1,1) = 1; B(1,2) = 5; B(1,3) = 7;
	
	Matrix C(1,3);
	C(1,1) = -1; C(1,2) = -4; C(1,3) = 3;
	
    _assert(m_equals(cross(A,B), C, 1e-10));
    
    return 0;
}

int extract_vector_01() {
	
	Matrix A(1,3);
	A(1,1) = 1; A(1,2) = 5; A(1,3) = 7;
	
	Matrix B =  extract_vector(A,2,3); 
	
	Matrix C(1,2);
	C(1,1) = 5; C(1,2) = 7;
	
    _assert(m_equals(B, C, 1e-10));
    
    return 0;
}

int union_vector_01() {
    Matrix A(1,2);
    A(1,1) = 4; A(1,2) = 7;
    
    Matrix B(1,3);
    B(1,1) = 1; B(1,2) = 5; B(1,3) = 7;
    
    Matrix C(1,5);
    C(1,1) = 1; C(1,2) = 5; C(1,3) = 7; C(1,4) = 4; C(1,5) = 7;
    
    Matrix R = union_vector(B,A);
    
    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}

int m_extract_row_01() {
	
	Matrix A(2,3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	A(2,1) = 4; A(2,2) = 5;
	
	Matrix B(1,3);
	B(1,1) = 1; B(1,2) = 2; B(1,3) = 3;
	
	Matrix C = extract_row(A,1);
	
    _assert(m_equals(B, C, 1e-10));
    
    return 0;
}

int m_extract_col_01() {
	
	Matrix A(2,3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	A(2,1) = 4; A(2,2) = 5;
	
	Matrix B(2,1);
	B(1,1) = 2; B(2,1) = 5;
	
	Matrix C = extract_column(A,2);
	
    _assert(m_equals(B, C, 1e-10));
    
    return 0;
}

int m_assign_row_01() {
	
	Matrix A(2,3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	A(2,1) = 4; A(2,2) = 5;
	
	Matrix B(1,2);
	B(1,1) = 2; B(1,2) = 5;
	
	Matrix C = assign_row(A,B,2);
	
	Matrix R(2,3);
	R(1,1) = 1; R(1,2) = 2; R(1,3) = 3;
	R(2,1) = 2; R(2,2) = 5;

    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}

int m_assign_col_01() {
	
	Matrix A(2, 4);
	A(1,1) = 4; A(1,2) = 7;
	A(2,1) = 2; A(2,2) = 6; A(2,3) = 2; A(2,4) = 6;
	
	Matrix B(3,1);
	B(1,1) = 1; 
	B(2,1) = 2; 
	B(3,1) = 3;
	
	Matrix C =  assign_column(A,B,2);
	
	Matrix R(3,4);
	R(1,1) = 4; R(1,2) = 1;
	R(2,1) = 2; R(2,2) = 2; R(2,3) = 2; R(2,4) = 6;
				R(3,2) = 3;

    _assert(m_equals(R, C, 1e-10));
    
    return 0;
}

int R_x_01() {
	Matrix A(3);
    A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = -0.44807361612917012694; A(2,3) = 0.89399666360055785042;
	A(3,1) = 0; A(3,2) = -0.89399666360055785042; A(3,3) = -0.44807361612917012694;
    
    Matrix B = R_x(90);
    
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int R_y_01() {
	Matrix A(3);
    A(1,1) = -0.44807361612917012694; A(1,2) = 0; A(1,3) = -0.89399666360055785042;
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 0.89399666360055785042; A(3,2) = 0; A(3,3) = -0.44807361612917012694;
    
    Matrix B = R_y(90);
    
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int R_z_01() {
	Matrix A(3);
    A(1,1) = -0.44807361612917012694; A(1,2) = 0.89399666360055785042; A(1,3) = 0;
	A(2,1) = -0.89399666360055785042; A(2,2) = -0.44807361612917012694; A(2,3) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;
    
    Matrix B = R_z(90);
    
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int AccelPointMass_01() {
	Matrix A(1,2);
	A(1,1) = 2; A(1,2) = 2;
	
	Matrix B(1,2);
	B(1,1) = 2; B(1,2) = 1;
    
	Matrix C(1,2);
	C(1,1) = -0.572433402239946; C(1,2) = -3.48621670111997;
    
	Matrix R = AccelPointMass(A,B,3.2);
	
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int all_tests()
{
    _verify(m_sum_01);
    _verify(m_sub_01);
	_verify(m_mult_01);
	_verify(m_div_01);
	_verify(m_assign_01);
    _verify(m_zeros_01);
	_verify(m_inv_01);
	_verify(m_eye_01);
	_verify(m_transpose_01);
	_verify(op_sum_01);
	_verify(op_sub_01);
	_verify(op_mult_01);
	_verify(op_div_01);
	_verify(m_norm_01);
	_verify(m_dot_01);
	_verify(m_cross_01);
	_verify(extract_vector_01);
	_verify(union_vector_01);
	_verify(m_extract_row_01);
	_verify(m_extract_col_01);
	_verify(m_assign_row_01);
	_verify(m_assign_col_01);
	_verify(R_x_01);
	_verify(R_y_01);
	_verify(R_z_01);
	_verify(AccelPointMass_01);
	
    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
