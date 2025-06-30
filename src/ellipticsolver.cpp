// Test file to solve elliptic function numerically

#include "configHeader.h"
#include "ellipticsolver.h"

double PerrinSolverPQR(double a, double b, double c)
{
    // Function that solves the elliptical integral for Perrin's frictional
    // coefficients. This is the equation for P given a, b, c but also can be
    // used for Q given b, a, c and for R given c, a, b

    a *= a;
    b *= b;
    c *= c;
    double elliptic = boost::math::ellint_rd(b, c, a);
    double result = 2 * elliptic / 3;
    return result;
}

double PerrinSolverS(double a, double b, double c)
{
    a *= a;
    b *= b;
    c *= c;

    double elliptic = boost::math::ellint_rf(a, b, c);
    double result = 2 * elliptic;
    return result;
}
