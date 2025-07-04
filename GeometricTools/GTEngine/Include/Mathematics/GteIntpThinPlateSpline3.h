// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteGMatrix.h>
#include <algorithm>
#include <cmath>
#include <limits>

// WARNING.  The implementation allows you to transform the inputs (x,y,z) to
// the unit cube and perform the interpolation in that space.  The idea is
// to keep the floating-point numbers to order 1 for numerical stability of
// the algorithm.  The classical thin-plate spline algorithm does not include
// this transformation.  The interpolation is invariant to translations and
// rotations of (x,y,z) but not to scaling.

namespace gte
{

template <typename Real>
class IntpThinPlateSpline3
{
public:
    // Construction.  Data points are (x,y,z,f(x,y,z)).  The smoothing parameter
    // must be nonnegative
    IntpThinPlateSpline3(int numPoints, Real const* X, Real const* Y,
        Real const* Z, Real const* F, Real smooth, bool transformToUnitCube);

    // Check this after the constructor call to see whether the thin plate
    // spline coefficients were successfully computed.  If so, then calls to
    // operator()(Real,Real,Real) will work properly.
    inline bool IsInitialized() const;

    // Evaluate the interpolator.  If IsInitialized()returns 'false', the
    // operator will return std::numeric_limits<Real>::max().
    Real operator()(Real x, Real y, Real z) const;

    // Compute the functional value a^T*M*a when lambda is zero or
    // lambda*w^T*(M+lambda*I)*w when lambda is positive.  See the thin plate
    // splines PDF for a description of these quantities.
    Real ComputeFunctional() const;

private:
    // Kernel(t) = -|t|
    static Real Kernel(Real t);

    // Input data.
    int mNumPoints;
    std::vector<Real> mX;
    std::vector<Real> mY;
    std::vector<Real> mZ;
    Real mSmooth;

    // Thin plate spline coefficients.  The A[] coefficients are associated
    // with the Green's functions G(x,y,z,*) and the B[] coefficients are
    // associated with the affine term B[0] + B[1]*x + B[2]*y + B[3]*z.
    std::vector<Real> mA;  // mNumPoints elements
    Real mB[4];

    // Extent of input data.
    Real mXMin, mXMax, mXInvRange;
    Real mYMin, mYMax, mYInvRange;
    Real mZMin, mZMax, mZInvRange;

    bool mInitialized;
};


template <typename Real>
IntpThinPlateSpline3<Real>::IntpThinPlateSpline3(int numPoints, Real const* X,
    Real const* Y, Real const* Z, Real const* F, Real smooth,
    bool transformToUnitCube)
    :
    mNumPoints(numPoints),
    mX(numPoints),
    mY(numPoints),
    mZ(numPoints),
    mSmooth(smooth),
    mA(numPoints),
    mInitialized(false)
{
    if (numPoints < 4 || !X || !Y || !Z || !F || smooth < (Real)0)
    {
        LogError("Invalid input.");
        return;
    }

    int i, row, col;

    if (transformToUnitCube)
    {
        // Map input (x,y,z) to unit cube.  This is not part of the classical
        // thin-plate spline algorithm, because the interpolation is not
        // invariant to scalings.
        auto extreme = std::minmax_element(X, X + mNumPoints);
        mXMin = *extreme.first;
        mXMax = *extreme.second;
        mXInvRange = ((Real)1) / (mXMax - mXMin);
        for (i = 0; i < mNumPoints; ++i)
        {
            mX[i] = (X[i] - mXMin) * mXInvRange;
        }

        extreme = std::minmax_element(Y, Y + mNumPoints);
        mYMin = *extreme.first;
        mYMax = *extreme.second;
        mYInvRange = ((Real)1) / (mYMax - mYMin);
        for (i = 0; i < mNumPoints; ++i)
        {
            mY[i] = (Y[i] - mYMin) * mYInvRange;
        }

        extreme = std::minmax_element(Z, Z + mNumPoints);
        mZMin = *extreme.first;
        mZMax = *extreme.second;
        mZInvRange = ((Real)1) / (mZMax - mZMin);
        for (i = 0; i < mNumPoints; ++i)
        {
            mZ[i] = (Z[i] - mZMin) * mZInvRange;
        }
    }
    else
    {
        // The classical thin-plate spline uses the data as is.  The values
        // mXMax, mYMax, and mZMax are not used, but they are initialized
        // anyway (to irrelevant numbers).
        mXMin = (Real)0;
        mXMax = (Real)1;
        mXInvRange = (Real)1;
        mYMin = (Real)0;
        mYMax = (Real)1;
        mYInvRange = (Real)1;
        mZMin = (Real)0;
        mZMax = (Real)1;
        mZInvRange = (Real)1;
        std::copy(X, X + mNumPoints, mX.begin());
        std::copy(Y, Y + mNumPoints, mY.begin());
        std::copy(Z, Z + mNumPoints, mZ.begin());
    }

    // Compute matrix A = M + lambda*I [NxN matrix].
    GMatrix<Real> AMat(mNumPoints, mNumPoints);
    for (row = 0; row < mNumPoints; ++row)
    {
        for (col = 0; col < mNumPoints; ++col)
        {
            if (row == col)
            {
                AMat(row, col) = mSmooth;
            }
            else
            {
                Real dx = mX[row] - mX[col];
                Real dy = mY[row] - mY[col];
                Real dz = mZ[row] - mZ[col];
                Real t = sqrt(dx * dx + dy * dy + dz * dz);
                AMat(row, col) = Kernel(t);
            }
        }
    }

    // Compute matrix B [Nx4 matrix].
    GMatrix<Real> BMat(mNumPoints, 4);
    for (row = 0; row < mNumPoints; ++row)
    {
        BMat(row, 0) = (Real)1;
        BMat(row, 1) = mX[row];
        BMat(row, 2) = mY[row];
        BMat(row, 3) = mZ[row];
    }

    // Compute A^{-1}.
    bool invertible;
    GMatrix<Real> invAMat = Inverse(AMat, &invertible);
    if (!invertible)
    {
        return;
    }

    // Compute P = B^t A^{-1} [4xN matrix].
    GMatrix<Real> PMat = MultiplyATB(BMat, invAMat);

    // Compute Q = P B = B^t A^{-1} B  [4x4 matrix].
    GMatrix<Real> QMat = PMat * BMat;

    // Compute Q^{-1}.
    GMatrix<Real> invQMat = Inverse(QMat, &invertible);
    if (!invertible)
    {
        return;
    }

    // Compute P*w.
    Real prod[4];
    for (row = 0; row < 4; ++row)
    {
        prod[row] = (Real)0;
        for (i = 0; i < mNumPoints; ++i)
        {
            prod[row] += PMat(row, i) * F[i];
        }
    }

    // Compute 'b' vector for smooth thin plate spline.
    for (row = 0; row < 4; ++row)
    {
        mB[row] = (Real)0;
        for (i = 0; i < 4; ++i)
        {
            mB[row] += invQMat(row, i) * prod[i];
        }
    }

    // Compute w-B*b.
    std::vector<Real> tmp(mNumPoints);
    for (row = 0; row < mNumPoints; ++row)
    {
        tmp[row] = F[row];
        for (i = 0; i < 4; ++i)
        {
            tmp[row] -= BMat(row, i) * mB[i];
        }
    }

    // Compute 'a' vector for smooth thin plate spline.
    for (row = 0; row < mNumPoints; ++row)
    {
        mA[row] = (Real)0;
        for (i = 0; i < mNumPoints; ++i)
        {
            mA[row] += invAMat(row, i) * tmp[i];
        }
    }

    mInitialized = true;
}

template <typename Real>
bool IntpThinPlateSpline3<Real>::IsInitialized() const
{
    return mInitialized;
}

template <typename Real>
Real IntpThinPlateSpline3<Real>::operator()(Real x, Real y, Real z) const
{
    if (mInitialized)
    {
        // Map (x,y,z) to the unit cube.
        x = (x - mXMin) * mXInvRange;
        y = (y - mYMin) * mYInvRange;
        z = (z - mZMin) * mZInvRange;

        Real result = mB[0] + mB[1] * x + mB[2] * y + mB[3] * z;
        for (int i = 0; i < mNumPoints; ++i)
        {
            Real dx = x - mX[i];
            Real dy = y - mY[i];
            Real dz = z - mZ[i];
            Real t = sqrt(dx*dx + dy*dy + dz*dz);
            result += mA[i] * Kernel(t);
        }
        return result;
    }

    return std::numeric_limits<Real>::max();
}

template <typename Real>
Real IntpThinPlateSpline3<Real>::ComputeFunctional() const
{
    Real functional = (Real)0;
    for (int row = 0; row < mNumPoints; ++row)
    {
        for (int col = 0; col < mNumPoints; ++col)
        {
            if (row == col)
            {
                functional += mSmooth * mA[row] * mA[col];
            }
            else
            {
                Real dx = mX[row] - mX[col];
                Real dy = mY[row] - mY[col];
                Real dz = mZ[row] - mZ[col];
                Real t = sqrt(dx * dx + dy * dy + dz * dz);
                functional += Kernel(t) * mA[row] * mA[col];
            }
        }
    }

    if (mSmooth > (Real)0)
    {
        functional *= mSmooth;
    }

    return functional;
}

template <typename Real>
Real IntpThinPlateSpline3<Real>::Kernel(Real t)
{
    return -std::abs(t);
}


}
