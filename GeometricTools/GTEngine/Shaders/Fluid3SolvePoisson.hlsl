// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

cbuffer Parameters
{
    float4 spaceDelta;    // (dx, dy, dz, 0)
    float4 halfDivDelta;  // (0.5/dx, 0.5/dy, 0.5/dz, 0)
    float4 timeDelta;     // (dt/dx, dt/dy, dt/dz, dt)
    float4 viscosityX;    // (velVX, velVX, velVX, denVX)
    float4 viscosityY;    // (velVX, velVY, velVY, denVY)
    float4 viscosityZ;    // (velVZ, velVZ, velVZ, denVZ)
    float4 epsilon;       // (epsilonX, epsilonY, epsilonZ, epsilon0)
};

Texture3D<float> divergence;
Texture3D<float> poisson;
RWTexture3D<float> outPoisson;

[numthreads(NUM_X_THREADS, NUM_Y_THREADS, NUM_Z_THREADS)]
void CSMain(uint3 c : SV_DispatchThreadID)
{
    uint3 dim;
    divergence.GetDimensions(dim.x, dim.y, dim.z);

    int x = int(c.x);
    int y = int(c.y);
    int z = int(c.z);
    int xm = max(x - 1, 0);
    int xp = min(x + 1, dim.x - 1);
    int ym = max(y - 1, 0);
    int yp = min(y + 1, dim.y - 1);
    int zm = max(z - 1, 0);
    int zp = min(z + 1, dim.z - 1);

    // Sample the divergence at (x,y,z).
    float div = divergence[int3(x, y, z)];

    // Sample Poisson values at (x,y) and immediate neighbors.
    float poisPZZ = poisson[int3(xp, y, z)];
    float poisMZZ = poisson[int3(xm, y, z)];
    float poisZPZ = poisson[int3(x, yp, z)];
    float poisZMZ = poisson[int3(x, ym, z)];
    float poisZZP = poisson[int3(x, y, zp)];
    float poisZZM = poisson[int3(x, y, zm)];

    float4 temp = float4(
        poisPZZ + poisMZZ, poisZPZ + poisZMZ, poisZZP + poisZZM, div);
    outPoisson[c.xyz] = dot(epsilon, temp);
}
