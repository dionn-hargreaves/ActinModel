// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.1 (2016/06/29)

#pragma once

#include <Mathematics/GteMatrix2x2.h>
#include <Mathematics/GteVector3.h>
#include <Mathematics/GteParametricSurface.h>
#include <memory>

namespace gte
{
template <typename Real>
class DarbouxFrame3
{
public:
    // Construction.  The curve must persist as long as the DarbouxFrame3
    // object does.
    DarbouxFrame3(std::shared_ptr<ParametricSurface<3, Real>> const& surface);

    // Get a coordinate frame, {T0, T1, N}.  At a nondegenerate surface
    // points, dX/du and dX/dv are linearly independent tangent vectors.
    // The frame is constructed as
    //   T0 = (dX/du)/|dX/du|
    //   N  = Cross(dX/du,dX/dv)/|Cross(dX/du,dX/dv)|
    //   T1 = Cross(N, T0)
    // so that {T0, T1, N} is a right-handed orthonormal set.
    void operator()(Real u, Real v, Vector3<Real>& position,
        Vector3<Real>& tangent0, Vector3<Real>& tangent1, Vector3<Real>& normal) const;

    // Compute the principal curvatures and principal directions.
    void GetPrincipalInformation(Real u, Real v, Real& curvature0, Real& curvature1,
        Vector3<Real>& direction0, Vector3<Real>& direction1) const;

private:
    std::shared_ptr<ParametricSurface<3, Real>> mSurface;
};


template <typename Real>
DarbouxFrame3<Real>::DarbouxFrame3(
    std::shared_ptr<ParametricSurface<3, Real>> const& surface)
    :
    mSurface(surface)
{
}

template <typename Real>
void DarbouxFrame3<Real>::operator()(Real u, Real v, Vector3<Real>& position,
    Vector3<Real>& tangent0, Vector3<Real>& tangent1, Vector3<Real>& normal) const
{
    Vector3<Real> values[6];
    mSurface->Evaluate(u, v, 1, values);
    position = values[0];
    tangent0 = values[1];
    Normalize(tangent0);
    tangent1 = values[2];
    Normalize(tangent1);
    normal = UnitCross(tangent0, tangent1);
    tangent1 = Cross(normal, tangent0);
}

template <typename Real>
void DarbouxFrame3<Real>::GetPrincipalInformation(Real u, Real v, Real& curvature0,
    Real& curvature1, Vector3<Real>& direction0, Vector3<Real>& direction1) const
{
    // Tangents:  T0 = (x_u,y_u,z_u), T1 = (x_v,y_v,z_v)
    // Normal:    N = Cross(T0,T1)/Length(Cross(T0,T1))
    // Metric Tensor:    G = +-                      -+
    //                       | Dot(T0,T0)  Dot(T0,T1) |
    //                       | Dot(T1,T0)  Dot(T1,T1) |
    //                       +-                      -+
    //
    // Curvature Tensor:  B = +-                          -+
    //                        | -Dot(N,T0_u)  -Dot(N,T0_v) |
    //                        | -Dot(N,T1_u)  -Dot(N,T1_v) |
    //                        +-                          -+
    //
    // Principal curvatures k are the generalized eigenvalues of
    //
    //     Bw = kGw
    //
    // If k is a curvature and w=(a,b) is the corresponding solution to
    // Bw = kGw, then the principal direction as a 3D vector is d = a*U+b*V.
    //
    // Let k1 and k2 be the principal curvatures.  The mean curvature
    // is (k1+k2)/2 and the Gaussian curvature is k1*k2.

    // Compute derivatives.
    Vector3<Real> values[6];
    mSurface->Evaluate(u, v, 2, values);
    Vector3<Real> derU = values[1];
    Vector3<Real> derV = values[2];
    Vector3<Real> derUU = values[3];
    Vector3<Real> derUV = values[4];
    Vector3<Real> derVV = values[5];

    // Compute the metric tensor.
    Matrix2x2<Real> metricTensor;
    metricTensor(0, 0) = Dot(values[1], values[1]);
    metricTensor(0, 1) = Dot(values[1], values[2]);
    metricTensor(1, 0) = metricTensor(0, 1);
    metricTensor(1, 1) = Dot(values[2], values[2]);

    // Compute the curvature tensor.
    Vector3<Real> normal = UnitCross(values[1], values[2]);
    Matrix2x2<Real> curvatureTensor;
    curvatureTensor(0, 0) = -Dot(normal, derUU);
    curvatureTensor(0, 1) = -Dot(normal, derUV);
    curvatureTensor(1, 0) = curvatureTensor(0, 1);
    curvatureTensor(1, 1) = -Dot(normal, derVV);

    // Characteristic polynomial is 0 = det(B-kG) = c2*k^2+c1*k+c0.
    Real c0 = Determinant(curvatureTensor);
    Real c1 = ((Real)2) * curvatureTensor(0, 1) * metricTensor(0, 1) -
        curvatureTensor(0, 0) * metricTensor(1, 1) -
        curvatureTensor(1, 1) * metricTensor(0, 0);
    Real c2 = Determinant(metricTensor);

    // Principal curvatures are roots of characteristic polynomial.
    Real temp = sqrt(std::abs(c1 * c1 - ((Real)4) * c0 * c2));
    Real mult = ((Real)0.5) / c2;
    curvature0 = -mult * (c1 + temp);
    curvature1 = -mult * (c1 - temp);

    // Principal directions are solutions to (B-kG)w = 0,
    // w1 = (b12-k1*g12,-(b11-k1*g11)) OR (b22-k1*g22,-(b12-k1*g12)).
    Real a0 = curvatureTensor(0, 1) - curvature0 * metricTensor(0, 1);
    Real a1 = curvature0 * metricTensor(0, 0) - curvatureTensor(0, 0);
    Real length = sqrt(a0 * a0 + a1 * a1);
    if (length > (Real)0)
    {
        direction0 = a0 * derU + a1 * derV;
    }
    else
    {
        a0 = curvatureTensor(1, 1) - curvature0 * metricTensor(1, 1);
        a1 = curvature0 * metricTensor(0, 1) - curvatureTensor(0, 1);
        length = sqrt(a0 * a0 + a1 * a1);
        if (length > (Real)0)
        {
            direction0 = a0*derU + a1*derV;
        }
        else
        {
            // Umbilic (surface is locally sphere, any direction principal).
            direction0 = derU;
        }
    }
    Normalize(direction0);

    // Second tangent is cross product of first tangent and normal.
    direction1 = Cross(direction0, normal);
}

}
