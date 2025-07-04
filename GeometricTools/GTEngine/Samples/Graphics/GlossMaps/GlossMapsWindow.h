// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <GTEngine.h>
#include "GlossMapEffect.h"
using namespace gte;

class GlossMapsWindow : public Window3
{
public:
    GlossMapsWindow(Parameters& parameters);

    virtual void OnIdle();

private:
    bool SetEnvironment();
    void CreateScene();
    void UpdateConstants();

    std::shared_ptr<Node> mScene;
    std::shared_ptr<Visual> mSquareNoGloss, mSquareGloss;
    std::shared_ptr<DirectionalLightEffect> mDLEffect;
    std::shared_ptr<GlossMapEffect> mGMEffect;
    Vector4<float> mLightWorldDirection;
};
