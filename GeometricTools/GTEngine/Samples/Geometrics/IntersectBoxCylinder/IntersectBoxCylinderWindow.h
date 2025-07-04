// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.4.0 (2016/11/08)

#pragma once

#include <GTEngine.h>
using namespace gte;

class IntersectBoxCylinderWindow : public Window3
{
public:
    IntersectBoxCylinderWindow(Parameters& parameters);

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
    std::shared_ptr<Visual> mCylinderMesh, mBoxMesh;
    std::shared_ptr<ConstantColorEffect> mRedEffect, mBlueEffect;
    Cylinder3<float> mCylinder;
    OrientedBox3<float> mBox;
    TIQuery<float, OrientedBox3<float>, Cylinder3<float>> mQuery;
};
