// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteAlignedBox.h>
#include <Mathematics/GteFIQuery.h>
#include <Mathematics/GteTIQuery.h>

// The queries consider the box to be a solid.
//
// The aligned-aligned queries use simple min-max comparisions.  The
// interesection of aligned boxes is an aligned box, possibly degenerate,
// where min[d] == max[d] for at least one dimension d.

namespace gte
{

template <typename Real>
class TIQuery<Real, AlignedBox3<Real>, AlignedBox3<Real>>
{
public:
    struct Result
    {
        bool intersect;
    };

    Result operator()(AlignedBox3<Real> const& box0,
        AlignedBox3<Real> const& box1);
};

template <typename Real>
class FIQuery<Real, AlignedBox3<Real>, AlignedBox3<Real>>
{
public:
    struct Result
    {
        bool intersect;
        AlignedBox3<Real> box;
    };

    Result operator()(AlignedBox3<Real> const& box0,
        AlignedBox3<Real> const& box1);
};


template <typename Real>
typename TIQuery<Real, AlignedBox3<Real>, AlignedBox3<Real>>::Result
TIQuery<Real, AlignedBox3<Real>, AlignedBox3<Real>>::operator()(
    AlignedBox3<Real> const& box0, AlignedBox3<Real> const& box1)
{
    Result result;
    for (int i = 0; i < 3; i++)
    {
        if (box0.max[i] < box1.min[i] || box0.min[i] > box1.max[i])
        {
            result.intersect = false;
            return result;
        }
    }
    result.intersect = true;
    return result;
}

template <typename Real>
typename FIQuery<Real, AlignedBox3<Real>, AlignedBox3<Real>>::Result
FIQuery<Real, AlignedBox3<Real>, AlignedBox3<Real>>::operator()(
    AlignedBox3<Real> const& box0, AlignedBox3<Real> const& box1)
{
    Result result;
    for (int i = 0; i < 3; i++)
    {
        if (box0.max[i] < box1.min[i] || box0.min[i] > box1.max[i])
        {
            result.intersect = false;
            return result;
        }
    }

    for (int i = 0; i < 3; i++)
    {
        if (box0.max[i] <= box1.max[i])
        {
            result.box.max[i] = box0.max[i];
        }
        else
        {
            result.box.max[i] = box1.max[i];
        }

        if (box0.min[i] <= box1.min[i])
        {
            result.box.min[i] = box1.min[i];
        }
        else
        {
            result.box.min[i] = box0.min[i];
        }
    }
    result.intersect = true;
    return result;
}


}
