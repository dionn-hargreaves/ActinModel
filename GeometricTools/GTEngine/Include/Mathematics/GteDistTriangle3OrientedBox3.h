// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.5.0 (2016/11/23)

#pragma once

#include <Mathematics/GteDCPQuery.h>
#include <Mathematics/GteLCPSolver.h>
#include <Mathematics/GteOrientedBox.h>
#include <Mathematics/GteTriangle.h>
#include <Mathematics/GteVector3.h>
#include <Mathematics/GteFunctions.h>

namespace gte
{

template <typename Real>
class DCPQuery<Real, Triangle3<Real>, OrientedBox3<Real>>
{
public:
    struct Result
    {
        bool queryIsSuccessful;

        // These members are valid only when queryIsSuccessful is true;
        // otherwise, they are all set to zero.
        Real distance, sqrDistance;
        std::array<Real, 3> triangleParameter, boxParameter;
        Vector3<Real> closestPoint[2];

        // The number of iterations used by LCPSolver regardless of whether
        // the query is successful.
        int numLCPIterations;
    };

    // The default maximum iterations is 81 (n = 9, maxIterations = n*n).  If
    // the solver fails to converge, try increasing the maximum number of
    // iterations.
    void SetMaxLCPIterations(int maxLCPIterations);

    Result operator()(Triangle3<Real> const& triangle, OrientedBox3<Real> const& box);

private:
    LCPSolver<Real, 9> mLCP;
};


template <typename Real>
void DCPQuery<Real, Triangle3<Real>, OrientedBox3<Real>>::SetMaxLCPIterations(int maxLCPIterations)
{
    mLCP.SetMaxIterations(maxLCPIterations);
}

template <typename Real>
typename DCPQuery<Real, Triangle3<Real>, OrientedBox3<Real>>::Result
DCPQuery<Real, Triangle3<Real>, OrientedBox3<Real>>::operator()(
    Triangle3<Real> const& triangle, OrientedBox3<Real> const& box)
{
    Result result;

    // Rigidly transform the triangle and oriented box so that the oriented
    // box becomes a canonical box.
    Vector3<Real> K = box.extent * (Real)2;
    Vector3<Real> tempV = triangle.v[0] - box.center;
    Vector3<Real> tempE0 = triangle.v[1] - triangle.v[0];
    Vector3<Real> tempE1 = triangle.v[2] - triangle.v[0];
    Vector3<Real> V, E0, E1;
    for (int i = 0; i < 3; ++i)
    {
        V[i] = Dot(box.axis[i], tempV) + box.extent[i];
        E0[i] = Dot(box.axis[i], tempE0);
        E1[i] = Dot(box.axis[i], tempE1);
    }

    // Compute quantities to initialize q and M in the LCP.
    Real dotVE0 = Dot(V, E0);
    Real dotVE1 = Dot(V, E1);
    Real dotE0E0 = Dot(E0, E0);
    Real dotE0E1 = Dot(E0, E1);
    Real dotE1E1 = Dot(E1, E1);

    // The LCP has 5 variables and 4 (nontrivial) inequality constraints.
    std::array<Real, 9> q =
    {
        -V[0], -V[1], -V[2], dotVE0, dotVE1, K[0], K[1], K[2], (Real)1
    };

    std::array<std::array<Real, 9>, 9> M;
    M[0] = { (Real)1, (Real)0, (Real)0, -E0[0], -E1[0], (Real)1, (Real)0, (Real)0, (Real)0 };
    M[1] = { (Real)0, (Real)1, (Real)0, -E0[1], -E1[1], (Real)0, (Real)1, (Real)0, (Real)0 };
    M[2] = { (Real)0, (Real)0, (Real)1, -E0[2], -E1[2], (Real)0, (Real)0, (Real)1, (Real)0 };
    M[3] = { -E0[0], -E0[1], -E0[2], dotE0E0, dotE0E1, (Real)0, (Real)0, (Real)0, (Real)1 };
    M[4] = { -E1[0], -E1[1], -E1[2], dotE0E1, dotE1E1, (Real)0, (Real)0, (Real)0, (Real)1 };
    M[5] = { (Real)-1, (Real)0, (Real)0, (Real)0, (Real)0, (Real)0, (Real)0, (Real)0, (Real)0 };
    M[6] = { (Real)0, (Real)-1, (Real)0, (Real)0, (Real)0, (Real)0, (Real)0, (Real)0, (Real)0 };
    M[7] = { (Real)0, (Real)0, (Real)-1, (Real)0, (Real)0, (Real)0, (Real)0, (Real)0, (Real)0 };
    M[8] = { (Real)0, (Real)0, (Real)0, (Real)-1, (Real)-1, (Real)0, (Real)0, (Real)0, (Real)0 };

    std::array<Real, 9> w, z;
    if (mLCP.Solve(q, M, w, z))
    {
        result.queryIsSuccessful = true;

        result.triangleParameter[0] = (Real)1 - z[3] - z[4];
        result.triangleParameter[1] = z[3];
        result.triangleParameter[2] = z[4];
        result.closestPoint[0] = triangle.v[0] + z[3] * tempE0 + z[4] * tempE1;
        result.closestPoint[1] = box.center;
        for (int i = 0; i < 3; ++i)
        {
            result.boxParameter[i] = z[i] - box.extent[i];
            result.closestPoint[1] += result.boxParameter[i] * box.axis[i];
        }

        Vector3<Real> diff = result.closestPoint[1] - result.closestPoint[0];
        result.sqrDistance = Dot(diff, diff);
        result.distance = Function<Real>::Sqrt(result.sqrDistance);
    }
    else
    {
        // If you reach this case, the maximum number of iterations was not
        // specified to be large enough or there is a problem due to
        // floating-point rounding errors.  If you believe the latter is
        // true, file a bug report.
        result.queryIsSuccessful = false;

        for (int i = 0; i < 3; ++i)
        {
            result.triangleParameter[i] = (Real)0;
            result.boxParameter[i] = (Real)0;
            result.closestPoint[0][i] = (Real)0;
            result.closestPoint[1][i] = (Real)0;
        }
        result.distance = (Real)0;
        result.sqrDistance = (Real)0;
    }

    result.numLCPIterations = mLCP.GetNumIterations();
    return result;
}

}
