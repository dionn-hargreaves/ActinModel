// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

StructuredBuffer<REAL> inBuffer;  // Two subnormal numbers.
RWStructuredBuffer<REAL> outBuffer;  // The sum of inputs, supposed to be subnormal.

[numthreads(1, 1, 1)]
void CSMain(int3 t : SV_DispatchThreadID)
{
    outBuffer[0] = inBuffer[0] + inBuffer[1];
}
