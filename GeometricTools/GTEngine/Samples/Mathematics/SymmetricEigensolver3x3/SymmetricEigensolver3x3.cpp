// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEngine.h>
#include <iostream>
using namespace gte;

double Determinant(std::array<std::array<double, 3>, 3> const& evec)
{
    double c00 = evec[1][1] * evec[2][2] - evec[1][2] * evec[2][1];
    double c01 = evec[1][0] * evec[2][2] - evec[1][2] * evec[2][0];
    double c02 = evec[1][0] * evec[2][1] - evec[1][1] * evec[2][0];
    return evec[0][0] * c00 - evec[0][1] * c01 + evec[0][2] * c02;
}

int main(int, char const*[])
{
    std::mt19937 mte;
    std::uniform_real_distribution<double> rnd(-1000.0, 1000.0);

    double a00, a01, a02, a11, a12, a22;
    std::array<double, 3> eval;
    std::array<std::array<double, 3>, 3> evec;
    SymmetricEigensolver3x3<double> esolver;
    int iters;
    double test[3];  // should be (0,0,0)
    double det;  // should be 1

    // Aggressive (14 iterations).
    a00 = rnd(mte);
    a01 = rnd(mte);
    a02 = 0.0;
    a11 = rnd(mte);
    a12 = 0.5 * a01;  // Exercise |b12| <= |b01| case in esolver.
    a22 = rnd(mte);
    iters = esolver(a00, a01, a02, a11, a12, a22, true, 0, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eigenvalues not sorted

    iters = esolver(a00, a01, a02, a11, a12, a22, true, +1, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eval[0] <= eval[1] <= eval[2]

    iters = esolver(a00, a01, a02, a11, a12, a22, true, -1, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eval[0] >= eval[1] >= eval[2]

    a12 = 2.0 * a01;  // Exercise |b12| > |b01| case in esolver.
    iters = esolver(a00, a01, a02, a11, a12, a22, true, 0, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eigenvalues not sorted

    iters = esolver(a00, a01, a02, a11, a12, a22, true, +1, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eval[0] <= eval[1] <= eval[2]

    iters = esolver(a00, a01, a02, a11, a12, a22, true, -1, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eval[0] >= eval[1] >= eval[2]

    // Not aggressive (8 iterations).
    a12 = 0.5 * a01;  // Exercise |b12| <= |b01| case in esolver.
    iters = esolver(a00, a01, a02, a11, a12, a22, false, 0, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eigenvalues not sorted

    iters = esolver(a00, a01, a02, a11, a12, a22, false, +1, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eval[0] <= eval[1] <= eval[2]

    iters = esolver(a00, a01, a02, a11, a12, a22, false, -1, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eval[0] >= eval[1] >= eval[2]

    a12 = 2.0 * a01;  // Exercise |b12| > |b01| case in esolver.
    iters = esolver(a00, a01, a02, a11, a12, a22, false, 0, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eigenvalues not sorted

    iters = esolver(a00, a01, a02, a11, a12, a22, false, +1, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eval[0] <= eval[1] <= eval[2]

    iters = esolver(a00, a01, a02, a11, a12, a22, false, -1, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);
    // eval[0] >= eval[1] >= eval[2]

    // A is already diagonal.
    a00 = 1.0;
    a01 = 0.0;
    a02 = 0.0;
    a11 = 2.0;
    a12 = 0.0;
    a22 = 3.0;
    iters = esolver(a00, a01, a02, a11, a12, a22, false, 0, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);

    // A is block diagonal.
    a00 = 1.0;
    a01 = 2.0;
    a02 = 0.0;
    a11 = 3.0;
    a12 = 0.0;
    a22 = 4.0;
    iters = esolver(a00, a01, a02, a11, a12, a22, false, 0, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);

    // A is zero.
    a00 = 0.0;
    a01 = 0.0;
    a02 = 0.0;
    a11 = 0.0;
    a12 = 0.0;
    a22 = 0.0;
    iters = esolver(a00, a01, a02, a11, a12, a22, false, 0, eval, evec);
    test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
    test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
    test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
    test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
    test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
    test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
    test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
    test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
    test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
    det = Determinant(evec);

    // Cases where a02 is not zero.
    double maxComponentError = 0.0;
    for (int i = 0; i < 1024; ++i)
    {
        a00 = rnd(mte);
        a01 = rnd(mte);
        a02 = rnd(mte);
        a11 = rnd(mte);
        a12 = rnd(mte);
        a22 = rnd(mte);

        iters = esolver(a00, a01, a02, a11, a12, a22, false, +1, eval, evec);
        double maxError = 0.0;

        test[0] = (a00 - eval[0])*evec[0][0] + a01*evec[0][1] + a02*evec[0][2];
        if (std::abs(test[0]) > maxError)
        {
            maxError = std::abs(test[0]);
        }
        test[1] = a01*evec[0][0] + (a11 - eval[0])*evec[0][1] + a12*evec[0][2];
        if (std::abs(test[1]) > maxError)
        {
            maxError = std::abs(test[1]);
        }
        test[2] = a02*evec[0][0] + a12*evec[0][1] + (a22 - eval[0])*evec[0][2];
        if (std::abs(test[2]) > maxError)
        {
            maxError = std::abs(test[2]);
        }

        test[0] = (a00 - eval[1])*evec[1][0] + a01*evec[1][1] + a02*evec[1][2];
        if (std::abs(test[0]) > maxError)
        {
            maxError = std::abs(test[0]);
        }
        test[1] = a01*evec[1][0] + (a11 - eval[1])*evec[1][1] + a12*evec[1][2];
        if (std::abs(test[1]) > maxError)
        {
            maxError = std::abs(test[1]);
        }
        test[2] = a02*evec[1][0] + a12*evec[1][1] + (a22 - eval[1])*evec[1][2];
        if (std::abs(test[2]) > maxError)
        {
            maxError = std::abs(test[2]);
        }

        test[0] = (a00 - eval[2])*evec[2][0] + a01*evec[2][1] + a02*evec[2][2];
        if (std::abs(test[0]) > maxError)
        {
            maxError = std::abs(test[0]);
        }
        test[1] = a01*evec[2][0] + (a11 - eval[2])*evec[2][1] + a12*evec[2][2];
        if (std::abs(test[1]) > maxError)
        {
            maxError = std::abs(test[1]);
        }
        test[2] = a02*evec[2][0] + a12*evec[2][1] + (a22 - eval[2])*evec[2][2];
        if (std::abs(test[2]) > maxError)
        {
            maxError = std::abs(test[2]);
        }

        det = Determinant(evec);

        std::cout << "i = " << i << ": det = " << det << " , maxError = "
            << maxError << std::endl;

        if (maxError > maxComponentError)
        {
            maxComponentError = maxError;
        }
    }
    std::cout << "maxError = " << maxComponentError << std::endl;
    return 0;
}
