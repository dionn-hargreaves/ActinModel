// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Graphics/GteViewVolume.h>
#include <Graphics/GteLighting.h>
#include <memory>

namespace gte
{

class GTE_IMPEXP Light : public ViewVolume
{
public:
    // Construction.  The depth range for DirectX is [0,1] and for OpenGL is
    // [-1,1].  For DirectX, set isDepthRangeZeroToOne to true.  For OpenGL,
    // set isDepthRangeZeroOne to false.
    Light(bool isPerspective, bool isDepthRangeZeroOne);

    std::shared_ptr<Lighting> lighting;
};

}
