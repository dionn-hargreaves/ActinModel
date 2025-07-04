// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <D3D11.h>
#include <D3Dcompiler.h>
#include <DXGI.h>
#include <cmath>
#include <fstream>

class Application
{
public:
    ~Application();
    Application();

    bool Create(int xOrigin, int yOrigin, int xSize, int ySize,
        D3D_FEATURE_LEVEL featureLevel, UINT flags);
    void Run();

private:
    static LRESULT CALLBACK WindowProcedure(HWND handle, UINT message,
        WPARAM wParam, LPARAM lParam);

    // Support for creation.
    bool CreateAppWindow(int xOrigin, int yOrigin, int xSize, int ySize);
    bool CreateGraphics(D3D_FEATURE_LEVEL featureLevel, UINT flags);
    bool CreateShaders();
    ID3DBlob* LoadShaderBlob(std::wstring const& filename);
    bool CreateVertexBuffers();
    bool CreateInputLayout();
    bool CreateConstantBuffers();
    bool CreateBlendState();
    bool CreateDepthStencilState();
    bool CreateRasterizerState();

    // Support for drawing.
    void ClearBuffers();
    void Draw();
    void SwapBuffers(unsigned int syncInterval);

    // Support for virtual trackball.
    void OnLeftMouseDown(int x, int y);
    void OnLeftMouseDrag(int x, int y);
    void OnLeftMouseUp(int x, int y);
    void RotateTrackball(float x0, float y0, float x1, float y1);
    void ComputeProjectionMatrix();
    void ComputeViewMatrix();
    void ComputeWorldMatrix(int whichTriangle);
    void UpdateConstantBuffers();

    // Window parameters.
    ATOM mAtom;
    HWND mHandle;
    int mXOrigin, mYOrigin, mXSize, mYSize;

    // D3D11 objects and parameters.
    ID3D11Device* mDevice;
    ID3D11DeviceContext* mImmediate;
    D3D_FEATURE_LEVEL mFeatureLevel;
    IDXGISwapChain* mSwapChain;
    ID3D11Texture2D* mColorBuffer;
    ID3D11RenderTargetView* mColorView;
    ID3D11Texture2D* mDepthStencilBuffer;
    ID3D11DepthStencilView* mDepthStencilView;
    D3D11_VIEWPORT mViewport;

    // Application-specific D3D11 objects.
    ID3D11VertexShader* mVertexShader;
    ID3DBlob* mVertexShaderBlob;
    ID3D11PixelShader* mPixelShader;
    ID3D11Buffer* mVertexBuffer[3];  // [0] for position, [1] for color0, [2] for color1
    ID3D11InputLayout* mInputLayout;
    ID3D11Buffer* mConstantBuffer[2];
    ID3D11BlendState* mBlendState;
    ID3D11DepthStencilState* mDepthStencilState;
    ID3D11RasterizerState* mRasterizerState;

    int mNumVertices;
    int mVertexOffset;

    // Camera parameters.
    float mUpFOVDegrees;
    float mAspectRatio;
    float mDMin, mDMax, mUMin, mUMax, mRMin, mRMax;
    float mEye[3], mDVector[3], mUVector[3], mRVector[3];
    float mViewMatrix[4][4];
    float mProjectionMatrix[4][4];
    float mWorldMatrix[2][4][4];  // 2 4x4 matrices, one per triangle

    // Data for a virtual trackball.
    float mTrackballMatrix[4][4];
    float mXTrack0, mYTrack0, mXTrack1, mYTrack1;
    float mSaveTrackballMatrix[4][4];
    bool mTrackBallDown;
};

extern Application* TheApplication;
