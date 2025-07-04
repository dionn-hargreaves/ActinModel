// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Physics/GteMassSpringCurve.h>
#include <Mathematics/GteVector3.h>
#include <random>
using namespace gte;

class PhysicsModule : public MassSpringCurve<3, float>
{
public:
    // Construction.  Gravity is controlled by the input 'gravity'.
    // Mass-spring systems tend to exhibit stiffness in the sense of numerical
    // stability.  To remedy this problem, a small amount of viscous friction
    // is added to the external force, -viscosity*velocity, where 'viscosity'
    // is a small positive constant.  The initial wind force is specified by
    // the caller.  The application of wind can be toggled by 'enableWind'.
    // The member 'enableWindChange' allows the wind direction to change
    // randomly, but each new direction is nearby the old direction in order
    // to obtain some sense of continuity of direction.  The magnitude of the
    // wind force is constant, the length of the initial force.
    PhysicsModule(int numParticles, float step, Vector3<float> const& gravity,
        Vector3<float> const& wind, float windChangeAmplitude,
        float viscosity);

    bool enableWind;
    bool enableWindChange;

    // External acceleration is due to forces of gravitation, wind, and
    // viscous friction.  The wind forces are randomly generated.
    virtual Vector<3, float> ExternalAcceleration(int i, float time,
        std::vector<Vector<3, float>> const& position,
        std::vector<Vector<3, float>> const& velocity);

protected:
    Vector3<float> mGravity, mWind;
    float mWindChangeAmplitude, mViscosity;
    std::mt19937 mMte;
    std::uniform_real_distribution<float> mRnd;
};
