// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.4.0 (2016/11/08)

#include "IntersectBoxCylinderWindow.h"

int main(int, char const*[])
{
#if defined(_DEBUG)
    LogReporter reporter(
        "LogReport.txt",
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ASSERTION | Logger::Listener::LISTEN_FOR_ERROR,
        Logger::Listener::LISTEN_FOR_ALL);
#endif

    Window::Parameters parameters(L"IntersectBoxCylinderWindow", 0, 0, 512, 512);
    auto window = TheWindowSystem.Create<IntersectBoxCylinderWindow>(parameters);
    TheWindowSystem.MessagePump(window, TheWindowSystem.DEFAULT_ACTION);
    TheWindowSystem.Destroy<IntersectBoxCylinderWindow>(window);
    return 0;
}

IntersectBoxCylinderWindow::IntersectBoxCylinderWindow(Parameters& parameters)
    :
    Window3(parameters)
{
    mNoCullState = std::make_shared<RasterizerState>();
    mNoCullState->cullMode = RasterizerState::CULL_NONE;
    mEngine->SetRasterizerState(mNoCullState);

    mNoCullWireState = std::make_shared<RasterizerState>();
    mNoCullWireState->cullMode = RasterizerState::CULL_NONE;
    mNoCullWireState->fillMode = RasterizerState::FILL_WIREFRAME;

    mBlendState = std::make_shared<BlendState>();
    mBlendState->target[0].enable = true;
    mBlendState->target[0].srcColor = BlendState::BM_SRC_ALPHA;
    mBlendState->target[0].dstColor = BlendState::BM_INV_SRC_ALPHA;
    mBlendState->target[0].srcAlpha = BlendState::BM_SRC_ALPHA;
    mBlendState->target[0].dstAlpha = BlendState::BM_INV_SRC_ALPHA;
    mEngine->SetBlendState(mBlendState);

    CreateScene();
    InitializeCamera(60.0f, GetAspectRatio(), 1.0f, 5000.0f, 0.1f, 0.01f,
        { 0.0f, 0.0f, -24.0f }, { 0.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 0.0f });

    TestIntersection();
    mPVWMatrices.Update();
}

void IntersectBoxCylinderWindow::OnIdle()
{
    mTimer.Measure();

    if (mCameraRig.Move())
    {
        mPVWMatrices.Update();
    }

    mEngine->ClearBuffers();
    mEngine->Draw(mCylinderMesh);
    mEngine->Draw(mBoxMesh);
    mEngine->Draw(8, mYSize - 8, { 0.0f, 0.0f, 0.0f, 1.0f }, mTimer.GetFPS());
    mEngine->DisplayColorBuffer(0);

    mTimer.UpdateFrameCount();
}

bool IntersectBoxCylinderWindow::OnCharPress(unsigned char key, int x, int y)
{
    float const delta = 0.1f;
    Quaternion<float> q, incr;

    switch (key)
    {
    case 'w':
    case 'W':
        if (mEngine->GetRasterizerState() == mNoCullState)
        {
            mEngine->SetRasterizerState(mNoCullWireState);
        }
        else
        {
            mEngine->SetRasterizerState(mNoCullState);
        }
        return true;

    case ' ':
        TestIntersection();
        return true;

    case 'x':  // decrement x-center of box
        Translate(0, -delta);
        return true;

    case 'X':  // increment x-center of box
        Translate(0, +delta);
        return true;

    case 'y':  // decrement y-center of box
        Translate(1, -delta);
        return true;

    case 'Y':  // increment y-center of box
        Translate(1, +delta);
        return true;

    case 'z':  // decrement z-center of box
        Translate(2, -delta);
        return true;

    case 'Z':  // increment z-center of box
        Translate(2, +delta);
        return true;

    case 'p':  // rotate about axis[0]
        Rotate(0, -delta);
        return true;

    case 'P':  // rotate about axis[0]
        Rotate(0, +delta);
        return true;

    case 'r':  // rotate about axis[1]
        Rotate(1, -delta);
        return true;

    case 'R':  // rotate about axis[1]
        Rotate(1, +delta);
        return true;

    case 'h':  // rotate about axis[2]
        Rotate(2, -delta);
        return true;

    case 'H':  // rotate about axis[2]
        Rotate(2, +delta);
        return true;
    }

    return Window::OnCharPress(key, x, y);
}

void IntersectBoxCylinderWindow::CreateScene()
{
    VertexFormat vformat;
    vformat.Bind(VA_POSITION, DF_R32G32B32_FLOAT, 0);
    MeshFactory mf;
    mf.SetVertexFormat(vformat);

    mCylinder.axis.origin = { 3.0f, 2.0f, 3.0f };
    mCylinder.axis.direction = { 0.0f, 0.0f, 1.0f };
    mCylinder.radius = 1.0f;
    mCylinder.height = 8.0f;

    mCylinderMesh = mf.CreateCylinderClosed(16, 16, mCylinder.radius, mCylinder.height);
    mCylinderMesh->localTransform.SetTranslation(mCylinder.axis.origin);
    auto effect = std::make_shared<ConstantColorEffect>(mProgramFactory,
            Vector4<float>({ 0.0f, 0.5f, 0.0f, 0.5f }));
    mCylinderMesh->SetEffect(effect);
    mPVWMatrices.Subscribe(mCylinderMesh->worldTransform, effect->GetPVWMatrixConstant());

    mRedEffect = std::make_shared<ConstantColorEffect>(mProgramFactory,
        Vector4<float>({ 0.5f, 0.0f, 0.0f, 0.5f }));

    mBlueEffect = std::make_shared<ConstantColorEffect>(mProgramFactory,
        Vector4<float>({ 0.0f, 0.0f, 0.5f, 0.5f }));

    mBox.center = { 0.0f, 0.0f, 0.0f };
    mBox.axis[0] = { 1.0f, 0.0f, 0.0f };
    mBox.axis[1] = { 0.0f, 1.0f, 0.0f };
    mBox.axis[2] = { 0.0f, 0.0f, 1.0f };
    mBox.extent = { 1.0f, 2.0f, 3.0f };

    mBoxMesh = mf.CreateBox(mBox.extent[0], mBox.extent[1], mBox.extent[2]);
    mBoxMesh->SetEffect(mBlueEffect);
    mPVWMatrices.Subscribe(mBoxMesh->worldTransform, mBlueEffect->GetPVWMatrixConstant());

    mTrackball.Attach(mCylinderMesh);
    mTrackball.Attach(mBoxMesh);
    mTrackball.Update();
}

void IntersectBoxCylinderWindow::Translate(int direction, float delta)
{
    mBox.center[direction] += delta;
    mBoxMesh->localTransform.SetTranslation(mBox.center);
    mBoxMesh->Update();
    TestIntersection();
    mPVWMatrices.Update();
}

void IntersectBoxCylinderWindow::Rotate(int direction, float delta)
{
    Quaternion<float> incr = Rotation<3, float>(
        AxisAngle<3, float>(mBox.axis[direction], delta));
    for (int i = 0; i < 3; ++i)
    {
        if (i != direction)
        {
            mBox.axis[i] = HProject(
                gte::Rotate(incr, HLift(mBox.axis[i], 0.0f)));
        }
    }
    Quaternion<float> q;
    mBoxMesh->localTransform.GetRotation(q);
    mBoxMesh->localTransform.SetRotation(incr * q);
    mBoxMesh->Update();
    TestIntersection();
    mPVWMatrices.Update();
}

void IntersectBoxCylinderWindow::TestIntersection()
{
    mPVWMatrices.Unsubscribe(mBoxMesh->worldTransform);

    if (mQuery(mBox, mCylinder).intersect)
    {
        mBoxMesh->SetEffect(mRedEffect);
        mPVWMatrices.Subscribe(mBoxMesh->worldTransform, mRedEffect->GetPVWMatrixConstant());
    }
    else
    {
        mBoxMesh->SetEffect(mBlueEffect);
        mPVWMatrices.Subscribe(mBoxMesh->worldTransform, mBlueEffect->GetPVWMatrixConstant());
    }
}
