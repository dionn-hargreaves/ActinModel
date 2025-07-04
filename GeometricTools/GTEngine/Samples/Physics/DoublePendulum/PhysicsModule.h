// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteVector4.h>
#include <Mathematics/GteOdeRungeKutta4.h>
#include <memory>
using namespace gte;

class PhysicsModule
{
public:
    // Construction.
    PhysicsModule();

    // Initialize the differential equation solver.
    void Initialize(float time, float deltaTime, float theta1, float theta2,
        float theta1Dot, float theta2Dot);

    // Access the current state.
    void GetPositions(float& x1, float& y1, float& x2, float& y2) const;

    // Apply a single step of the solver.
    void Update();

    // The physical constants.
    float gravity;
    float mass1, mass2;
    float length1, length2;
    float jointX, jointY;

protected:
    // State and auxiliary variables.
    Vector4<float> mState;
    float mTime, mAux[4];

    // Runge-Kutta 4th-order ODE solver.
    typedef OdeRungeKutta4<float, Vector4<float>> Solver;
    std::unique_ptr<Solver> mSolver;
};
