// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <GTEngine.h>
using namespace gte;

class PlaneEstimationWindow : public Window
{
public:
    PlaneEstimationWindow(Parameters& parameters);

    virtual void OnIdle() override;

private:
    bool SetEnvironment();
    std::shared_ptr<ConstantBuffer> CreateBezierControls();

    std::string mShaderSourceEvaluateBezier;
    std::string mShaderSourcePlaneEstimation;
    std::string mShaderSourcePlaneVisualize;
    std::string mShaderSourcePositionVisualize;
    std::shared_ptr<Texture2> mPositions;
    std::shared_ptr<Texture2> mPlanes;
    std::shared_ptr<ComputeProgram> mPositionProgram;
    std::shared_ptr<ComputeProgram> mPlaneProgram;
    unsigned int mNumXGroups, mNumYGroups;
    std::shared_ptr<OverlayEffect> mOverlay[2];
};
