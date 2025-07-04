// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

cbuffer Weights
{
    float weight[2 * RADIUS + 1];
};

Texture2D<float4> inImage;
RWTexture2D<float4> outImage;

// Maximum size of groupshared memory in D3D11 is 32768 bytes.  Each float4
// uses 16 bytes, so the maximum number of samples is 2048.  The maximum
// number of threads in a group is 1024.  The test image is 1024x768, so if
// we choose NUM_X_THREADS = 1024, then NUM_Y_THREADS = 1 is required.  See
// ConvolveSeparableHGS2.hlsl for a variation in which NUM_X_THREADS < 1024.
#define NUM_X_THREADS 1024
groupshared float4 samples[NUM_X_THREADS + 2*RADIUS];

[numthreads(NUM_X_THREADS, 1, 1)]
void CSMain(int2 dt : SV_DispatchThreadID, int2 gt : SV_GroupThreadID)
{
    // Load the texels from the input texture, store them in group-shared
    // memory, and have all threads in the group wait until all texels
    // are loaded.
    samples[gt.x + RADIUS] = inImage[dt];
    if (gt.x < RADIUS)
    {
        samples[gt.x] = inImage[dt - int2(RADIUS, 0)];
    }
    else if (gt.x >= NUM_X_THREADS - RADIUS)
    {
        samples[gt.x + 2 * RADIUS] = inImage[dt + int2(RADIUS, 0)];
    }
    GroupMemoryBarrierWithGroupSync();

    float4 result = 0.0f;
    for (int x = 0; x <= 2*RADIUS; ++x)
    {
        result += weight[x] * samples[gt.x + x];
    }
    outImage[dt] = result;
}
