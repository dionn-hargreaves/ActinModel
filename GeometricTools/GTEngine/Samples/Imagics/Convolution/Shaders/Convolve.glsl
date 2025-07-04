// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#define WSIZE (2 * RADIUS + 1)

uniform Weights
{
    // (2R+1)-by-(2R+1) mask stored in row-major order.
    float weight[WSIZE * WSIZE];
};

layout(rgba32f) uniform readonly image2D inImage;
layout(rgba32f) uniform writeonly image2D outImage;

layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;
void main()
{
    ivec2 dt = ivec2(gl_GlobalInvocationID.xy);
    vec4 result = vec4(0.0f, 0.0f, 0.0f, 0.0f);
    for (int y = -RADIUS; y <= RADIUS; ++y)
    {
        for (int x = -RADIUS; x <= RADIUS; ++x)
        {
            vec4 value = imageLoad(inImage, dt + ivec2(x, y));
            result += weight[(x + RADIUS) + WSIZE * (y + RADIUS)] * value;
        }
    }
    imageStore(outImage, dt, result);
}
