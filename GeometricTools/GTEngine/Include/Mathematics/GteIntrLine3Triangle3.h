// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteVector3.h>
#include <Mathematics/GteLine.h>
#include <Mathematics/GteTriangle.h>
#include <Mathematics/GteFIQuery.h>
#include <Mathematics/GteTIQuery.h>

namespace gte
{

template <typename Real>
class TIQuery<Real,Line3<Real>,Triangle3<Real>>
{
public:
    struct Result
    {
        bool intersect;
    };

    Result operator()(Line3<Real> const& line,
        Triangle3<Real> const& triangle);
};

template <typename Real>
class FIQuery<Real, Line3<Real>, Triangle3<Real>>
{
public:
    struct Result
    {
        Result();

        bool intersect;
        Real parameter;
        Real triangleBary[3];
        Vector3<Real> point;
    };

    Result operator()(Line3<Real> const& line,
        Triangle3<Real> const& triangle);
};


template <typename Real>
typename TIQuery<Real, Line3<Real>, Triangle3<Real>>::Result
TIQuery<Real, Line3<Real>, Triangle3<Real>>::operator()(
    Line3<Real> const& line, Triangle3<Real> const& triangle)
{
    Result result;

    // Compute the offset origin, edges, and normal.
    Vector3<Real> diff = line.origin - triangle.v[0];
    Vector3<Real> edge1 = triangle.v[1] - triangle.v[0];
    Vector3<Real> edge2 = triangle.v[2] - triangle.v[0];
    Vector3<Real> normal = Cross(edge1, edge2);

    // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = line direction,
    // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
    //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
    //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
    //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
    Real DdN = Dot(line.direction, normal);
    Real sign;
    if (DdN > (Real)0)
    {
        sign = (Real)1;
    }
    else if (DdN < (Real)0)
    {
        sign = (Real)-1;
        DdN = -DdN;
    }
    else
    {
        // Line and triangle are parallel, call it a "no intersection"
        // even if the line and triangle are coplanar and intersecting.
        result.intersect = false;
        return result;
    }

    Real DdQxE2 = sign*DotCross(line.direction, diff, edge2);
    if (DdQxE2 >= (Real)0)
    {
        Real DdE1xQ = sign*DotCross(line.direction, edge1, diff);
        if (DdE1xQ >= (Real)0)
        {
            if (DdQxE2 + DdE1xQ <= DdN)
            {
                // Line intersects triangle.
                result.intersect = true;
                return result;
            }
            // else: b1+b2 > 1, no intersection
        }
        // else: b2 < 0, no intersection
    }
    // else: b1 < 0, no intersection

    result.intersect = false;
    return result;
}

template <typename Real>
FIQuery<Real, Line3<Real>, Triangle3<Real>>::Result::Result()
    :
    parameter((Real)0)
{
    triangleBary[0] = (Real)0;
    triangleBary[1] = (Real)0;
    triangleBary[2] = (Real)0;
}

template <typename Real>
typename FIQuery<Real, Line3<Real>, Triangle3<Real>>::Result
FIQuery<Real, Line3<Real>, Triangle3<Real>>::operator()(
    Line3<Real> const& line, Triangle3<Real> const& triangle)
{
    Result result;

    // Compute the offset origin, edges, and normal.
    Vector3<Real> diff = line.origin - triangle.v[0];
    Vector3<Real> edge1 = triangle.v[1] - triangle.v[0];
    Vector3<Real> edge2 = triangle.v[2] - triangle.v[0];
    Vector3<Real> normal = Cross(edge1, edge2);

    // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = line direction,
    // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
    //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
    //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
    //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
    Real DdN = Dot(line.direction, normal);
    Real sign;
    if (DdN > (Real)0)
    {
        sign = (Real)1;
    }
    else if (DdN < (Real)0)
    {
        sign = (Real)-1;
        DdN = -DdN;
    }
    else
    {
        // Line and triangle are parallel, call it a "no intersection"
        // even if the line and triangle are coplanar and intersecting.
        result.intersect = false;
        return result;
    }

    Real DdQxE2 = sign*DotCross(line.direction, diff, edge2);
    if (DdQxE2 >= (Real)0)
    {
        Real DdE1xQ = sign*DotCross(line.direction, edge1, diff);
        if (DdE1xQ >= (Real)0)
        {
            if (DdQxE2 + DdE1xQ <= DdN)
            {
                // Line intersects triangle.
                Real QdN = -sign*Dot(diff, normal);
                Real inv = ((Real)1) / DdN;

                result.intersect = true;
                result.parameter = QdN*inv;
                result.triangleBary[1] = DdQxE2*inv;
                result.triangleBary[2] = DdE1xQ*inv;
                result.triangleBary[0] = (Real)1 - result.triangleBary[1]
                    - result.triangleBary[2];
                result.point = line.origin +
                    result.parameter * line.direction;
                return result;
            }
            // else: b1+b2 > 1, no intersection
        }
        // else: b2 < 0, no intersection
    }
    // else: b1 < 0, no intersection

    result.intersect = false;
    return result;
}


}
