// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <GTEngine.h>
using namespace gte;

class MedianFilteringWindow : public Window
{
public:
    MedianFilteringWindow(Parameters& parameters);

    virtual void OnIdle() override;
    virtual bool OnCharPress(unsigned char key, int x, int y) override;

private:
    bool SetEnvironment();
    bool CreatePrograms(unsigned int txWidth, unsigned int txHeight);

    std::shared_ptr<Texture2> mOriginal;
    std::shared_ptr<Texture2> mImage[2];
    std::shared_ptr<OverlayEffect> mOverlay[2];

    // 0 = median 3x3 by insertion sort
    // 1 = median 3x3 by min-max
    // 2 = median 5x5 by insertion sort
    // 3 = median 5x5 by min-max
    unsigned int mSelection;
    std::shared_ptr<ComputeProgram> mMedianProgram[4];
    std::shared_ptr<ComputeProgram> mCProgram;
    unsigned int mNumXGroups, mNumYGroups;
    static std::string msName[4];
};
