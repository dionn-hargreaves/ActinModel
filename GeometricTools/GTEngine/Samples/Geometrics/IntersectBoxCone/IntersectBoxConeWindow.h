// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <GTEngine.h>
using namespace gte;

#define USE_ORIENTED_BOX

class IntersectBoxConeWindow : public Window3
{
public:
    IntersectBoxConeWindow(Parameters& parameters);

    virtual void OnIdle() override;
    virtual bool OnCharPress(unsigned char key, int x, int y) override;

private:
    void CreateScene();
    void Translate(int direction, float delta);
    void Rotate(int direction, float delta);
    void TestIntersection();

    std::shared_ptr<RasterizerState> mNoCullState;
    std::shared_ptr<RasterizerState> mNoCullWireState;
    std::shared_ptr<BlendState> mBlendState;
    std::shared_ptr<Visual> mConeMesh, mDiskMesh, mBoxMesh;
    std::shared_ptr<ConstantColorEffect> mRedEffect, mBlueEffect;
    Cone<3, float> mCone;
#if defined(USE_ORIENTED_BOX)
    TIOrientedBox3Cone3<float> mQuery;
    OrientedBox<3, float> mBox;
#else
    TIAlignedBox3Cone3<float> mQuery;
    AlignedBox<3, float> mBox;
#endif
};
