// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include "SharedTexturesWindow.h"

SharedTexturesWindow::~SharedTexturesWindow()
{
    mSharedTexture = nullptr;
    mOverlay = nullptr;
    mSecondEngine = nullptr;
}

SharedTexturesWindow::SharedTexturesWindow (Parameters& parameters)
    :
    Window3(parameters)
{
    // The first engine is mEngine, provided by the Window class.  Create a
    // second engine that owns a texture.  The access must be GPU_RW for
    // the texture to be shared by the first engine.
    int const tsize = 16;
    mSecondEngine = std::make_shared<DX11Engine>();
    mSharedTexture = std::make_shared<Texture2>(DF_R8G8B8A8_UNORM, tsize, tsize);
    mSharedTexture->SetUsage(Resource::SHADER_OUTPUT);

    // Set the texture values to random numbers.
    char* data = mSharedTexture->GetData();
    std::mt19937 mte;
    std::uniform_int_distribution<int> rnd(0, 255);
    for (size_t i = 0; i < mSharedTexture->GetNumBytes(); ++i)
    {
        *data++ = static_cast<char>(rnd(mte));
    }

    // Specify the texture must be shared *before* telling DX11 to create
    // a GPU resource for it.
    mSharedTexture->MakeShared();
    mSecondEngine->Bind(mSharedTexture);

    // Tell the first engine to share the texture owned by the second engine.
    std::static_pointer_cast<DX11Engine>(mEngine)->Share(mSharedTexture, mSecondEngine.get());

    // Ensure that mEngine gets the texture data that was copied to the GPU
    // via mSecondEngine->Bind(mSharedTexture).
    mSecondEngine->Flush();

    // The first engine will display the shared texture.
    mOverlay = std::make_shared<OverlayEffect>(mProgramFactory, mXSize,
        mYSize, tsize, tsize, SamplerState::MIN_L_MAG_L_MIP_P,
        SamplerState::CLAMP, SamplerState::CLAMP, true);
    mOverlay->SetTexture(mSharedTexture);
}

void SharedTexturesWindow::OnIdle ()
{
    mTimer.Measure();

    mEngine->ClearBuffers();
    mEngine->Draw(mOverlay);
    mEngine->Draw(8, mYSize - 8, { 0.0f, 0.0f, 0.0f, 1.0f }, mTimer.GetFPS());
    mEngine->DisplayColorBuffer(0);

    mTimer.UpdateFrameCount();
}
