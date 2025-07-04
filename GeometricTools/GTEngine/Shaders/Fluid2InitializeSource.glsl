// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

uniform Parameters
{
    vec4 spaceDelta;    // (dx, dy, 0, 0)
    vec4 halfDivDelta;  // (0.5/dx, 0.5/dy, 0, 0)
    vec4 timeDelta;     // (dt/dx, dt/dy, 0, dt)
    vec4 viscosityX;    // (velVX, velVX, 0, denVX)
    vec4 viscosityY;    // (velVX, velVY, 0, denVY)
    vec4 epsilon;       // (epsilonX, epsilonY, 0, epsilon0)
};

uniform External
{
    vec4 densityProducer;  // (x, y, variance, amplitude)
    vec4 densityConsumer;  // (x, y, variance, amplitude)
    vec4 gravity;          // (x, y, *, *)
    vec4 wind;             // (x, y, variance, amplitude)
};

layout(rg32f) uniform readonly image2D vortexVelocity;
layout(rgba32f) uniform writeonly image2D source;

layout (local_size_x = NUM_X_THREADS, local_size_y = NUM_Y_THREADS, local_size_z = 1) in;
void main()
{
    ivec2 c = ivec2(gl_GlobalInvocationID.xy);

    // Compute the location of the pixel (x,y) in normalized [0,1]^2.
    vec2 location = spaceDelta.xy*(c + 0.5f);

    // Compute an input to the fluid simulation consisting of a producer of
    // density and a consumer of density.
    vec2 diff = location - densityProducer.xy;
    float arg = -dot(diff, diff) / densityProducer.z;
    float density = densityProducer.w*exp(arg);
    diff = location - densityConsumer.xy;
    arg = -dot(diff, diff) / densityConsumer.z;
    density -= densityConsumer.w*exp(arg);

    // Compute an input to the fluid simulation consisting of gravity,
    // a single wind source, and vortex impulses.
    float windDiff = location.y - wind.y;
    float windArg = -windDiff*windDiff / wind.z;
    vec2 windVelocity = vec2(wind.w*exp(windArg), 0.0f);
    vec2 velocity = gravity.xy + windVelocity + imageLoad(vortexVelocity, c).xy;

    imageStore(source, c, vec4(velocity, 0.0f, density));
};
