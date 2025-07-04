// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include "BumpMapsWindow.h"
#include "SimpleBumpMapEffect.h"

int main(int, char const*[])
{
#if defined(_DEBUG)
    LogReporter reporter(
        "LogReport.txt",
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL);
#endif

    Window::Parameters parameters(L"BumpMapsWindow", 0, 0, 640, 480);
    auto window = TheWindowSystem.Create<BumpMapsWindow>(parameters);
    TheWindowSystem.MessagePump(window, TheWindowSystem.DEFAULT_ACTION);
    TheWindowSystem.Destroy<BumpMapsWindow>(window);
    return 0;
}

BumpMapsWindow::BumpMapsWindow(Parameters& parameters)
    :
    Window3(parameters),
    mUseBumpMap(true)
{
    if (!SetEnvironment() || !CreateScene())
    {
        parameters.created = false;
        return;
    }

    InitializeCamera(60.0f, GetAspectRatio(), 0.1f, 100.0f, 0.01f, 0.001f,
        { 0.0f, -0.25f, -2.5f }, { 0.0f, 0.0f, 1.0f }, { 0.0f, 1.0f, 0.0f });
    mPVWMatrices.Update();
    UpdateBumpMap();
}

void BumpMapsWindow::OnIdle()
{
    mTimer.Measure();

    if (mCameraRig.Move())
    {
        mPVWMatrices.Update();
    }

    mEngine->ClearBuffers();
    mEngine->Draw(mTorus);
    mEngine->Draw(8, mYSize - 8, { 0.0f, 0.0f, 0.0f, 1.0f }, mTimer.GetFPS());
    mEngine->DisplayColorBuffer(0);

    mTimer.UpdateFrameCount();
}

bool BumpMapsWindow::OnCharPress(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 'b':
    case 'B':
    {
        mUseBumpMap = !mUseBumpMap;
        mPVWMatrices.Unsubscribe(mTorus->worldTransform);
        Transform save = mTorus->localTransform;
        CreateTorus();
        mTorus->localTransform = save;
        mScene->SetChild(0, mTorus);
        mTrackball.Update();
        mPVWMatrices.Update();
        UpdateBumpMap();
        return true;
    }
    }

    return Window::OnCharPress(key, x, y);
}

bool BumpMapsWindow::OnMouseMotion(MouseButton button, int x, int y, unsigned int modifiers)
{
    if (Window3::OnMouseMotion(button, x, y, modifiers))
    {
        mPVWMatrices.Update();
        UpdateBumpMap();
    }
    return true;
}

bool BumpMapsWindow::SetEnvironment()
{
    std::string path = GetGTEPath();
    if (path == "")
    {
        return false;
    }

    mEnvironment.Insert(path + "/Samples/Graphics/BumpMaps/Shaders/");
    mEnvironment.Insert(path + "/Samples/Data/");
    std::vector<std::string> inputs =
    {
        "Bricks.png",
        "BricksNormal.png"
    };

    for (auto const& input : inputs)
    {
        if (mEnvironment.GetPath(input) == "")
        {
            LogError("Cannot find file " + input);
            return false;
        }
    }

    return true;
}

bool BumpMapsWindow::CreateScene()
{
    if (CreateTorus())
    {
        mScene = std::make_shared<Node>();
        mTorus->localTransform.SetRotation(
            AxisAngle<4, float>(Vector4<float>::Unit(0), (float)GTE_C_QUARTER_PI));
        mScene->AttachChild(mTorus);
        mTrackball.Attach(mScene);
        mTrackball.Update();
        return true;
    }
    else
    {
        return false;
    }
}

bool BumpMapsWindow::CreateTorus()
{
    struct Vertex
    {
        Vector3<float> position;
        Vector3<float> normal;
        Vector3<float> lightDirection;
        Vector2<float> baseTCoord;
        Vector2<float> normalTCoord;
    };

    VertexFormat vformat;
    vformat.Bind(VA_POSITION, DF_R32G32B32_FLOAT, 0);
    vformat.Bind(VA_NORMAL, DF_R32G32B32_FLOAT, 0);
    vformat.Bind(VA_COLOR, DF_R32G32B32_FLOAT, 0);
    vformat.Bind(VA_TEXCOORD, DF_R32G32_FLOAT, 0);
    vformat.Bind(VA_TEXCOORD, DF_R32G32_FLOAT, 1);

    MeshFactory mf;
    mf.SetVertexFormat(vformat);
    mf.SetVertexBufferUsage(Resource::DYNAMIC_UPDATE);
    mTorus = mf.CreateTorus(32, 32, 1.0f, 0.4f);
    std::shared_ptr<VertexBuffer> vbuffer = mTorus->GetVertexBuffer();
    unsigned int const numVertices = vbuffer->GetNumElements();
    Vertex* vertex = vbuffer->Get<Vertex>();
    for (unsigned int i = 0; i < numVertices; ++i)
    {
        vertex[i].baseTCoord *= 4.0f;
        if (mUseBumpMap)
        {
            vertex[i].normalTCoord *= 4.0f;
        }
    }

    if (mUseBumpMap)
    {
        bool created = false;
        std::shared_ptr<SimpleBumpMapEffect> effect =
            std::make_shared<SimpleBumpMapEffect>(mProgramFactory,
            mEnvironment, created);
        if (!created)
        {
            LogError("Failed to create the bump map effect.");
            return false;
        }

        mTorus->SetEffect(effect);
        mPVWMatrices.Subscribe(mTorus->worldTransform, effect->GetPVWMatrixConstant());

        mLightDirection = Vector4<float>::Unit(2);
        SimpleBumpMapEffect::ComputeLightVectors(mTorus, mLightDirection);
    }
    else
    {
        std::string baseName = mEnvironment.GetPath("Bricks.png");
        std::shared_ptr<Texture2> baseTexture = WICFileIO::Load(baseName, true);
        baseTexture->AutogenerateMipmaps();

        std::shared_ptr<Texture2Effect> effect =
            std::make_shared<Texture2Effect>(mProgramFactory, baseTexture,
            SamplerState::MIN_L_MAG_L_MIP_L, SamplerState::WRAP, SamplerState::WRAP);

        mTorus->SetEffect(effect);
        mPVWMatrices.Subscribe(mTorus->worldTransform, effect->GetPVWMatrixConstant());
    }

    return true;
}

void BumpMapsWindow::UpdateBumpMap()
{
    if (mUseBumpMap)
    {
        // The scene graph transformations have been updated, which means the
        // tangent-space light vectors need updating.
        SimpleBumpMapEffect::ComputeLightVectors(mTorus, mLightDirection);
        mEngine->Update(mTorus->GetVertexBuffer());
    }
}
