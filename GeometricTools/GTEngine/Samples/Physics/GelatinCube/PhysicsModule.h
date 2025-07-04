// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.3.0 (2016/07/10)

#pragma once

#include <Physics/GteMassSpringVolume.h>
#include <Mathematics/GteVector3.h>
using namespace gte;

class PhysicsModule : public MassSpringVolume<3, float>
{
public:
    // Construction.
    PhysicsModule(int numSlices, int numRows, int numCols, float step, float viscosity);

    // External acceleration is due to viscous forces.
    virtual Vector3<float> ExternalAcceleration(int i, float time,
        std::vector<Vector3<float>> const& positions,
        std::vector<Vector3<float>> const& velocities);

protected:
    float mViscosity;
};
