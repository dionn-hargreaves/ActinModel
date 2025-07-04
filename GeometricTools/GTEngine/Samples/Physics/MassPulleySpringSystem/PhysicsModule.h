// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.3.0 (2016/07/13)

#pragma once

class PhysicsModule
{
public:
    // Construction.
    PhysicsModule();

    // Initialize the differential equation solver.
    void Initialize(float time, float deltaTime, float y1, float dy1, float dy3);

    // Apply a single step of the solver.
    void Update();

    // The gravitational constant.
    float gravity;

    // The left mass in Figure 3.13.
    float mass1;

    // The right mass in Figure 3.13.
    float mass2;

    // The length of rigid wire connecting mass 1 to mass 2.
    float wireLength;

    // The pulley parameters.
    float mass3, radius, inertia;

    // The spring parameters.
    float springLength, springConstant;

    // Member access.
    inline float GetCurrentY1() const;
    inline float GetCurrentY2() const;
    inline float GetCurrentY3() const;
    inline float GetAngle() const;
    inline float GetCableFraction1() const;
    inline float GetCableFraction2() const;

private:
    // time information
    float mTime, mDeltaTime;

    // derived parameters
    float mAlpha, mBeta, mGamma, mDelta, mOmega, mGDivOmegaSqr;

    // initial conditions
    float mY1, mY2, mY3, mDY1, mDY2, mDY3;

    // solution parameters
    float mLPlusGDivOmegaSqr, mK1, mK2, mTCoeff, mTSqrCoeff;
    float mDeltaDivOmegaSqr, mY1Curr, mY2Curr, mY3Curr;
};


inline float PhysicsModule::GetCurrentY1() const
{
    return mY1Curr;
}

inline float PhysicsModule::GetCurrentY2() const
{
    return mY2Curr;
}

inline float PhysicsModule::GetCurrentY3() const
{
    return mY3Curr;
}

inline float PhysicsModule::GetAngle() const
{
    return (mY1 - mY1Curr) / radius;
}

inline float PhysicsModule::GetCableFraction1() const
{
    return (mY1Curr - mY3Curr) / wireLength;
}

inline float PhysicsModule::GetCableFraction2() const
{
    return (mY2Curr - mY3Curr) / wireLength;
}
