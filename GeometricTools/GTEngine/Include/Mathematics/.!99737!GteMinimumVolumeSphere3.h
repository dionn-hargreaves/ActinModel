// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <LowLevel/GteLogger.h>
#include <Mathematics/GteHypersphere.h>
#include <Mathematics/GteLinearSystem.h>
#include <algorithm>
#include <random>

// Compute the minimum volume sphere containing the input set of points.  The
// algorithm randomly permutes the input points so that the construction
// occurs in 'expected' O(N) time.  All internal minimal sphere calculations
// store the squared radius in the radius member of Sphere3.  Only at
// the end is a sqrt computed.

namespace gte
{

template <typename InputType, typename ComputeType>
class MinimumVolumeSphere3
{
public:
    bool operator()(int numPoints, Vector3<InputType> const* points,
        Sphere3<InputType>& minimal);

    // Member access.
    inline int GetNumSupport() const;
    inline std::array<int, 4> const& GetSupport() const;

private:
    // Test whether point P is inside sphere S using squared distance and
    // squared radius.
    bool Contains(int i, Sphere3<ComputeType> const& sphere) const;

    Sphere3<ComputeType> ExactSphere1(int i0) const;
    Sphere3<ComputeType> ExactSphere2(int i0, int i1) const;
    Sphere3<ComputeType> ExactSphere3(int i0, int i1, int i2) const;
    Sphere3<ComputeType> ExactSphere4(int i0, int i1, int i2, int i3) const;

    Sphere3<ComputeType> UpdateSupport1(int i);
    Sphere3<ComputeType> UpdateSupport2(int i);
    Sphere3<ComputeType> UpdateSupport3(int i);
    Sphere3<ComputeType> UpdateSupport4(int i);

    // Indices of points that support current minimum volume sphere.
    bool SupportContains(int j) const;

    int mNumSupport;
    std::array<int, 4> mSupport;

    // Random permutation of the unique input points to produce expected
    // linear time for the algorithm.
    std::vector<Vector3<ComputeType>> mComputePoints;
};


template <typename InputType, typename ComputeType>
bool MinimumVolumeSphere3<InputType, ComputeType>::operator()(int numPoints,
    Vector3<InputType> const* points, Sphere3<InputType>& minimal)
{
    if (numPoints >= 1 && points)
    {
        // Function array to avoid switch statement in the main loop.
        std::function<Sphere3<ComputeType>(int)> update[5];
        update[1] = [this](int i) { return UpdateSupport1(i); };
        update[2] = [this](int i) { return UpdateSupport2(i); };
        update[3] = [this](int i) { return UpdateSupport3(i); };
        update[4] = [this](int i) { return UpdateSupport4(i); };

        // Process only the unique points.
        std::vector<int> permuted(numPoints);
        for (int i = 0; i < numPoints; ++i)
        {
            permuted[i] = i;
        }
        std::sort(permuted.begin(), permuted.end(),
            [points](int i0, int i1) { return points[i0] < points[i1]; });
        auto end = std::unique(permuted.begin(), permuted.end(),
            [points](int i0, int i1) { return points[i0] == points[i1]; });
        permuted.erase(end, permuted.end());
        numPoints = static_cast<int>(permuted.size());

        // Create a random permutation of the points.
        std::shuffle(permuted.begin(), permuted.end(),
            std::default_random_engine());

        // Convert to the compute type, which is a simple copy when
        // ComputeType is the same as InputType.
        mComputePoints.resize(numPoints);
        for (int i = 0; i < numPoints; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                mComputePoints[i][j] = points[permuted[i]][j];
            }
        }

        // Start with the first point.
        Sphere3<ComputeType> ctMinimal = ExactSphere1(0);
        mNumSupport = 1;
        mSupport[0] = 0;

        // The loop restarts from the beginning of the point list each time
