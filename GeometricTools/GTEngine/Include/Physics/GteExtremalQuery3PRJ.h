// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Physics/GteExtremalQuery3.h>
#include <limits>

namespace gte
{

template <typename Real>
class ExtremalQuery3PRJ : public ExtremalQuery3<Real>
{
public:
    // Construction.
    ExtremalQuery3PRJ(Polyhedron3<Real> const& polytope);

    // Disallow copying and assignment.
    ExtremalQuery3PRJ(ExtremalQuery3PRJ const&) = delete;
    ExtremalQuery3PRJ& operator=(ExtremalQuery3PRJ const&) = delete;

    // Compute the extreme vertices in the specified direction and return the
    // indices of the vertices in the polyhedron vertex array.
    virtual void GetExtremeVertices(Vector3<Real> const& direction,
        int& positiveDirection, int& negativeDirection) override;

private:
    Vector3<Real> mCentroid;
};


template <typename Real>
ExtremalQuery3PRJ<Real>::ExtremalQuery3PRJ(Polyhedron3<Real> const& polytope)
    :
    ExtremalQuery3<Real>(polytope)
{
    mCentroid = this->mPolytope.ComputeVertexAverage();
}

template <typename Real>
void ExtremalQuery3PRJ<Real>::GetExtremeVertices(Vector3<Real> const& direction,
    int& positiveDirection, int& negativeDirection)
{
    Real minValue = std::numeric_limits<Real>::max(), maxValue = -minValue;
    negativeDirection = -1;
    positiveDirection = -1;

    auto vertexPool = this->mPolytope.GetVertexPool();
    for (auto i : this->mPolytope.GetUniqueIndices())
    {
        Vector3<Real> diff = vertexPool.get()->at(i) - mCentroid;
        Real dot = Dot(direction, diff);
        if (dot < minValue)
        {
            negativeDirection = i;
            minValue = dot;
        }
        if (dot > maxValue)
        {
            positiveDirection = i;
            maxValue = dot;
        }
    }
}

}
