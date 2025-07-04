// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include "Delaunay2DWindow.h"

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

    Window::Parameters parameters(L"Delaunay2DWindow", 0, 0, 512, 512);
    auto window = TheWindowSystem.Create<Delaunay2DWindow>(parameters);
    TheWindowSystem.MessagePump(window, TheWindowSystem.NO_IDLE_LOOP);
    TheWindowSystem.Destroy<Delaunay2DWindow>(window);
    return 0;
}

Delaunay2DWindow::Delaunay2DWindow(Parameters& parameters)
    :
    Window2(parameters),
    mCurrentTriX(-1),
    mCurrentTriY(-1),
    mCurrentIndex(0)
{
#if 1
    // Randomly generated points.
    std::mt19937 mte;
    std::uniform_real_distribution<float> rnd(0.125f, 0.875f);
    mVertices.resize(256);
    for (auto& v : mVertices)
    {
        v[0] = mXSize*rnd(mte);
        v[1] = mYSize*rnd(mte);
    }
#endif

#if 0
    // A 3x3 square grid.
    mVertices.resize(9);
    mVertices[0] = Vector2<float>( 64.0f,  64.0f);
    mVertices[1] = Vector2<float>( 64.0f, 256.0f);
    mVertices[2] = Vector2<float>( 64.0f, 448.0f);
    mVertices[3] = Vector2<float>(256.0f,  64.0f);
    mVertices[4] = Vector2<float>(256.0f, 256.0f);
    mVertices[5] = Vector2<float>(256.0f, 448.0f);
    mVertices[6] = Vector2<float>(448.0f,  64.0f);
    mVertices[7] = Vector2<float>(448.0f, 256.0f);
    mVertices[8] = Vector2<float>(448.0f, 448.0f);
#endif

    mDelaunay(static_cast<int>(mVertices.size()), &mVertices[0], 0.001f);
    if (mDelaunay.GetDimension() == 2)
    {
        mDelaunay.GetHull(mHull);
    }
    else
    {
        LogError("Degenerate point set.");
        parameters.created = false;
        return;
    }

    mInfo.initialTriangle = -1;
    mInfo.finalTriangle = 0;
}

void Delaunay2DWindow::OnDisplay()
{
    ClearScreen(0xFFFFFFFF);

    unsigned int const white = 0xFFFFFFFF;
    unsigned int const gray = 0xFF808080;
    unsigned int const red = 0xFF0000FF;
    unsigned int const blue = 0xFFFF0000;
    unsigned int const green = 0xFF00FF00;
    unsigned int const black = 0xFF000000;

    int i, x0, y0, x1, y1, x2, y2;
    Vector2<float> v0, v1, v2;

    // Draw the triangle mesh.
    std::vector<int> const& indices = mDelaunay.GetIndices();
    int numTriangles = static_cast<int>(indices.size() / 3);
    for (i = 0; i < numTriangles; ++i)
    {
        v0 = mVertices[indices[3 * i]];
        x0 = static_cast<int>(std::lrint(v0[0]));
        y0 = static_cast<int>(std::lrint(v0[1]));

        v1 = mVertices[indices[3 * i + 1]];
        x1 = static_cast<int>(std::lrint(v1[0]));
        y1 = static_cast<int>(std::lrint(v1[1]));

        v2 = mVertices[indices[3 * i + 2]];
        x2 = static_cast<int>(std::lrint(v2[0]));
        y2 = static_cast<int>(std::lrint(v2[1]));

        DrawLine(x0, y0, x1, y1, gray);
        DrawLine(x1, y1, x2, y2, gray);
        DrawLine(x2, y2, x0, y0, gray);
    }

    // Draw the hull.
    int numEdges = static_cast<int>(mHull.size() / 2);
    for (i = 0; i < numEdges; ++i)
    {
        v0 = mVertices[mHull[2 * i]];
        x0 = static_cast<int>(std::lrint(v0[0]));
        y0 = static_cast<int>(std::lrint(v0[1]));

        v1 = mVertices[mHull[2 * i + 1]];
        x1 = static_cast<int>(std::lrint(v1[0]));
        y1 = static_cast<int>(std::lrint(v1[1]));

        DrawLine(x0, y0, x1, y1, red);
    }

    // Draw the search path.
    if (mInfo.initialTriangle != -1)
    {
        for (i = 0; i < mInfo.numPath; ++i)
        {
            int index = mInfo.path[i];

            v0 = mVertices[indices[3 * index]];
            v1 = mVertices[indices[3 * index + 1]];
            v2 = mVertices[indices[3 * index + 2]];

            Vector2<float> center = (v0 + v1 + v2) / 3.0f;
            int x = static_cast<int>(center[0] + 0.5f);
            int y = static_cast<int>(center[1] + 0.5f);
            if (i < mInfo.numPath - 1)
            {
                DrawFloodFill4(x, y, blue, white);
            }
            else
            {
                DrawFloodFill4(x, y, red, white);
            }
        }
    }

    if (mCurrentTriX >= 0)
    {
        // Draw the current triangle.
        DrawFloodFill4(mCurrentTriX, mCurrentTriY, green, red);
    }
    else if (mInfo.initialTriangle != -1)
    {
        // Draw the last edge when the selected point is outside the hull.
        v0 = mVertices[mInfo.finalV[0]];
        x0 = static_cast<int>(std::lrint(v0[0]));
        y0 = static_cast<int>(std::lrint(v0[1]));

        v1 = mVertices[mInfo.finalV[1]];
        x1 = static_cast<int>(std::lrint(v1[0]));
        y1 = static_cast<int>(std::lrint(v1[1]));

        DrawLine(x0, y0, x1, y1, black);
    }

    mScreenTextureNeedsUpdate = true;
    Window2::OnDisplay();
}

bool Delaunay2DWindow::OnMouseClick(MouseButton button, MouseState state,
    int x, int y, unsigned int)
{
    if (button != MOUSE_LEFT)
    {
        return false;
    }

    if (state == MOUSE_DOWN)
    {
        Vector2<float> pos{ (float)x, (float)y };
        mInfo.initialTriangle = mInfo.finalTriangle;
        int i = mDelaunay.GetContainingTriangle(pos, mInfo);
        if (i >= 0)
        {
            mCurrentTriX = x;
            mCurrentTriY = y;
            mInfo.finalTriangle = i;
        }
        else
        {
            mCurrentTriX = -1;
            mCurrentTriY = -1;
        }
        OnDisplay();
    }

    return true;
}
