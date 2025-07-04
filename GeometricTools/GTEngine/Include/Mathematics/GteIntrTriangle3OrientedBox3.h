// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteTIQuery.h>
#include <Mathematics/GteFIQuery.h>
#include <Mathematics/GteTriangle.h>
#include <Mathematics/GteOrientedBox.h>
#include <Mathematics/GteIntrConvexPolygonPlane.h>

namespace gte
{

template <typename Real>
class TIQuery<Real, Triangle3<Real>, OrientedBox3<Real>>
{
public:
    struct Result
    {
        bool intersect;
    };

    Result operator()(Triangle3<Real> const& triangle, OrientedBox3<Real> const& box);

private:
    void GetProjection(Vector3<Real> const& axis, Triangle3<Real> const& triangle,
        Real& imin, Real& imax);

    void GetProjection(Vector3<Real> const& axis, OrientedBox3<Real> const& box,
        Real& imin, Real& imax);
};

template <typename Real>
class FIQuery<Real, Triangle3<Real>, OrientedBox3<Real>>
{
public:
    struct Result
    {
        std::vector<Vector3<Real>> insidePolygon;
        std::vector<std::vector<Vector3<Real>>> outsidePolygons;
    };

    Result operator()(Triangle3<Real> const& triangle, OrientedBox3<Real> const& box);
};


template <typename Real>
typename TIQuery<Real, Triangle3<Real>, OrientedBox3<Real>>::Result
TIQuery<Real, Triangle3<Real>, OrientedBox3<Real>>::operator()(
    Triangle3<Real> const& triangle, OrientedBox3<Real> const& box)
{
    // The method of separating axes is the algorithm implemented here.
    Result result;

    Real min0, max0, min1, max1;
    Vector3<Real> D, edge[3];

    // Test direction of triangle normal.
    edge[0] = triangle.v[1] - triangle.v[0];
    edge[1] = triangle.v[2] - triangle.v[0];
    D = Cross(edge[0], edge[1]);
    min0 = Dot(D, triangle.v[0]);
    max0 = min0;
    GetProjection(D, box, min1, max1);
    if (max1 < min0 || max0 < min1)
    {
        result.intersect = false;
        return result;
    }

    // Test direction of box faces.
    for (int i = 0; i < 3; ++i)
    {
        D = box.axis[i];
        GetProjection(D, triangle, min0, max0);
        Real DdC = Dot(D, box.center);
        min1 = DdC - box.extent[i];
        max1 = DdC + box.extent[i];
        if (max1 < min0 || max0 < min1)
        {
            result.intersect = false;
            return result;
        }
    }

    // Test direction of triangle-box edge cross products.
    edge[2] = edge[1] - edge[0];
    for (int i0 = 0; i0 < 3; ++i0)
    {
        for (int i1 = 0; i1 < 3; ++i1)
        {
            D = Cross(edge[i0], box.axis[i1]);
            GetProjection(D, triangle, min0, max0);
            GetProjection(D, box, min1, max1);
            if (max1 < min0 || max0 < min1)
            {
                result.intersect = false;
                return result;
            }
        }
    }

    result.intersect = true;
    return result;
}

template <typename Real>
void TIQuery<Real, Triangle3<Real>, OrientedBox3<Real>>::GetProjection(
    Vector3<Real> const& axis, Triangle3<Real> const& triangle, Real& imin, Real& imax)
{
    Real dot[3] =
    {
        Dot(axis, triangle.v[0]),
        Dot(axis, triangle.v[1]),
        Dot(axis, triangle.v[2])
    };

    imin = dot[0];
    imax = imin;

    if (dot[1] < imin)
    {
        imin = dot[1];
    }
    else if (dot[1] > imax)
    {
        imax = dot[1];
    }

    if (dot[2] < imin)
    {
        imin = dot[2];
    }
    else if (dot[2] > imax)
    {
        imax = dot[2];
    }
}

template <typename Real>
void TIQuery<Real, Triangle3<Real>, OrientedBox3<Real>>::GetProjection(
    Vector3<Real> const& axis, OrientedBox3<Real> const& box, Real& imin, Real& imax)
{
    Real origin = Dot(axis, box.center);
    Real maximumExtent =
        fabs(box.extent[0] * Dot(axis, box.axis[0])) +
        fabs(box.extent[1] * Dot(axis, box.axis[1])) +
        fabs(box.extent[2] * Dot(axis, box.axis[2]));

    imin = origin - maximumExtent;
    imax = origin + maximumExtent;
}

template <typename Real>
typename FIQuery<Real, Triangle3<Real>, OrientedBox3<Real>>::Result
FIQuery<Real, Triangle3<Real>, OrientedBox3<Real>>::operator()(
    Triangle3<Real> const& triangle, OrientedBox3<Real> const& box)
{
    Result result;

    // Start with the triangle and clip it against each face of the box.
    // The largest number of vertices for the polygon of intersection is 7.
    result.insidePolygon.resize(3);
    for (int i = 0; i < 3; ++i)
    {
        result.insidePolygon[i] = triangle.v[i];
    }

    Plane3<Real> plane;
    FIQuery<Real, std::vector<Vector3<Real>>, Plane3<Real>> ppQuery;
    typename FIQuery<Real, std::vector<Vector3<Real>>, Plane3<Real>>::Result ppResult;
    for (int dir = -1; dir <= 1; dir += 2)
    {
        for (int side = 0; side < 3; ++side)
        {
            // Create a plane for the box face that points inside the box.
            plane.normal = ((Real)dir) * box.axis[side];
            plane.constant = Dot(plane.normal, box.center) - box.extent[side];

            ppResult = ppQuery(result.insidePolygon, plane);
            switch (ppResult.configuration)
            {
            case ConvexPolygonPlaneConfiguration::SPLIT_BY_PLANE:
                result.insidePolygon = ppResult.positivePolygon;
                result.outsidePolygons.push_back(ppResult.negativePolygon);
                break;
            case ConvexPolygonPlaneConfiguration::POSITIVE_SIDE_VERTEX:
            case ConvexPolygonPlaneConfiguration::POSITIVE_SIDE_EDGE:
            case ConvexPolygonPlaneConfiguration::POSITIVE_SIDE_STRICT:
                // The result.insidePolygon is already ppResult.positivePolygon.
                //result.insidePolygon = ppResult.positivePolygon;
                break;
            case ConvexPolygonPlaneConfiguration::NEGATIVE_SIDE_VERTEX:
            case ConvexPolygonPlaneConfiguration::NEGATIVE_SIDE_EDGE:
            case ConvexPolygonPlaneConfiguration::NEGATIVE_SIDE_STRICT:
                result.insidePolygon.clear();
                result.outsidePolygons.push_back(ppResult.negativePolygon);
                return result;
            case ConvexPolygonPlaneConfiguration::COINCIDENT_WITH_PLANE:
                // A triangle coplanar with a box face will be processed as if
                // it is inside the box.
                result.insidePolygon = ppResult.intersection;
                break;
            default:
                result.insidePolygon.clear();
                result.outsidePolygons.clear();
                break;
            }
        }
    }

    return result;
}

}
