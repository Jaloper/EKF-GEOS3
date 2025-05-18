#include "..\include\sign_.hpp"

double sign_(double a, double b) {
    if (b >= 0.0)
        return abs(a);
    else
        return -abs(a);
}