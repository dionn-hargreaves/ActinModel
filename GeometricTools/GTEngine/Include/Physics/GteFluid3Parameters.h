// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#pragma once

#include <Mathematics/GteVector4.h>

namespace gte
{

struct GTE_IMPEXP Fluid3Parameters
{
    Vector4<float> spaceDelta;    // (dx, dy, dz, 0)
    Vector4<float> halfDivDelta;  // (0.5/dx, 0.5/dy, 0.5/dz, 0)
    Vector4<float> timeDelta;     // (dt/dx, dt/dy, dt/dz, dt)
    Vector4<float> viscosityX;    // (velVX, velVX, velVX, denVX)
    Vector4<float> viscosityY;    // (velVX, velVY, velVY, denVY)
    Vector4<float> viscosityZ;    // (velVZ, velVZ, velVZ, denVZ)
    Vector4<float> epsilon;       // (epsX, epsY, epsZ, eps0)
};

}
