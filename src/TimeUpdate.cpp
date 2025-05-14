#include "..\include\TimeUpdate.hpp"


Matrix TimeUpdate(Matrix P, Matrix Phi, Matrix Qdt) {
       Matrix temp = Phi * P;
    Matrix result = temp * transpose(Phi);
    result = result + Qdt;
    P = result;
    return P;
}

Matrix TimeUpdate(Matrix P, Matrix Phi) {
        Matrix temp = Phi * P;
    Matrix result = temp * transpose(Phi);
    P = result;
    return P;
}