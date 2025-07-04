// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteVector3.h>
#include <Mathematics/GteHypersphere.h>

namespace gte
{

// Compute the smallest sphere whose center is the average of the input
// points.
template <typename Real>
bool GetContainer(int numPoints, Vector3<Real> const* points,
    Sphere3<Real>& sphere);

// Test for containment of a point inside a sphere.
template <typename Real>
bool InContainer(Vector3<Real> const& point, Sphere3<Real> const& sphere);

// Compute the smallest sphere that contains the input spheres.
template <typename Real>
bool MergeContainers(Sphere3<Real> const& sphere0,
    Sphere3<Real> const& sphere1, Sphere3<Real>& merge);


template <typename Real>
bool GetContainer(int numPoints, Vector3<Real> const* points,
    Sphere3<Real>& sphere)
{
    sphere.center = points[0];
    for (int i = 1; i < numPoints; ++i)
    {
        sphere.center += points[i];
    }
    sphere.center /= (Real)numPoints;

    for (int i = 0; i < numPoints; ++i)
    {
        Vector3<Real> diff = points[i] - sphere.center;
        Real radiusSqr = Dot(diff, diff);
        if (radiusSqr > sphere.radius)
        {
            sphere.radius = radiusSqr;
        }
    }

    sphere.radius = sqrt(sphere.radius);
    return true;
}

template <typename Real>
bool InContainer(Vector3<Real> const& point, Sphere3<Real> const& sphere)
{
    Vector3<Real> diff = point - sphere.center;
    return Length(diff) <= sphere.radius;
}

template <typename Real>
bool MergeContainers(Sphere3<Real> const& sphere0,
    Sphere3<Real> const& sphere1, Sphere3<Real>& merge)
{
    Vector3<Real> cenDiff = sphere1.center - sphere0.center;
    Real lenSqr = Dot(cenDiff, cenDiff);
    Real rDiff = sphere1.radius - sphere0.radius;
    Real rDiffSqr = rDiff*rDiff;

    if (rDiffSqr >= lenSqr)
    {
        merge = (rDiff >= (Real)0 ? sphere1 : sphere0);
    }
    else
    {
        Real length = sqrt(lenSqr);
        if (length > (Real)0)
        {
            Real coeff = (length + rDiff) / (((Real)2)*length);
            merge.center = sphere0.center + coeff*cenDiff;
        }
        else
        {
            merge.center = sphere0.center;
        }

        merge.radius = ((Real)0.5)*(length + sphere0.radius + sphere1.radius);
    }

    return true;
}


}
