// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteApprOrthogonalLine3.h>
#include <Mathematics/GteDistPointLine.h>
#include <Mathematics/GteCylinder3.h>

namespace gte
{

// Compute the cylinder axis segment using least-squares fit.  The radius is
// the maximum distance from points to the axis.  The height is determined by
// projection of points onto the axis and determining the containing interval.
template <typename Real>
bool GetContainer(int numPoints, Vector3<Real> const* points,
     Cylinder3<Real>& cylinder);

// Test for containment of a point by a cylinder.
template <typename Real>
bool InContainer(Vector3<Real> const& point, Cylinder3<Real> const& cylinder);


template <typename Real>
bool GetContainer(int numPoints, Vector3<Real> const* points,
    Cylinder3<Real>& cylinder)
{
    ApprOrthogonalLine3<Real> fitter;
    fitter.Fit(numPoints, points);
    Line3<Real> line = fitter.GetParameters();

    DCPQuery<Real, Vector3<Real>, Line3<Real>> plQuery;
    Real maxRadiusSqr = (Real)0;
    for (int i = 0; i < numPoints; ++i)
    {
        auto result = plQuery(points[i], line);
        if (result.sqrDistance > maxRadiusSqr)
        {
            maxRadiusSqr = result.sqrDistance;
        }
    }

    Vector3<Real> diff = points[0] - line.origin;
    Real wMin = Dot(line.direction, diff);
    Real wMax = wMin;
    for (int i = 1; i < numPoints; ++i)
    {
        diff = points[i] - line.origin;
        Real w = Dot(line.direction, diff);
        if (w < wMin)
        {
            wMin = w;
        }
        else if (w > wMax)
        {
            wMax = w;
        }
    }

    cylinder.axis.origin = line.origin +
        (((Real)0.5)*(wMax + wMin))*line.direction;
    cylinder.axis.direction = line.direction;
    cylinder.radius = sqrt(maxRadiusSqr);
    cylinder.height = wMax - wMin;
    return true;
}

template <typename Real>
bool InContainer(Vector3<Real> const& point, Cylinder3<Real> const& cylinder)
{
    Vector3<Real> diff = point - cylinder.axis.origin;
    Real zProj = Dot(diff, cylinder.axis.direction);
    if (std::abs(zProj)*(Real)2 > cylinder.height)
    {
        return false;
    }

    Vector3<Real> xyProj = diff - zProj*cylinder.axis.direction;
    return Length(xyProj) <= cylinder.radius;
}


}
