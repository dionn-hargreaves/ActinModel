// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

cbuffer PVWMatrix
{
    float4x4 pvwMatrix;
};

cbuffer VWMatrix
{
    float4x4 vwMatrix;
};

struct VS_INPUT
{
    float3 modelPosition : POSITION;
    float3 modelNormal : NORMAL;
};

struct VS_OUTPUT
{
    float2 vertexTCoord : TEXCOORD0;
    float4 clipPosition : SV_POSITION;
};

VS_OUTPUT VSMain(VS_INPUT input)
{
    VS_OUTPUT output;

    float4 modelPosition = float4(input.modelPosition, 1.0f);
    float4 modelNormal = float4(input.modelNormal, 0.0f);

#if GTE_USE_MAT_VEC
    float4 cameraSpacePosition = mul(vwMatrix, modelPosition);
    float3 cameraSpaceNormal = normalize(mul(vwMatrix, modelNormal).xyz);
    output.clipPosition = mul(pvwMatrix, modelPosition);
#else
    float4 cameraSpacePosition = mul(modelPosition, vwMatrix);
    float3 cameraSpaceNormal = normalize(mul(modelNormal, vwMatrix).xyz);
    output.clipPosition = mul(modelPosition, pvwMatrix);
#endif

    float3 eyeDirection = normalize(cameraSpacePosition.xyz);
    float3 r = reflect(eyeDirection, cameraSpaceNormal);

    float oneMRZ = 1.0f - r.z;
    float invLength = 1.0f / sqrt(r.x * r.x + r.y * r.y + oneMRZ * oneMRZ);
    output.vertexTCoord = 0.5f * (r.xy * invLength + 1.0f);

    return output;
}

struct PS_INPUT
{
    float2 vertexTCoord : TEXCOORD0;
};

struct PS_OUTPUT
{
    float4 pixelColor : SV_TARGET0;
};

Texture2D<float4> baseTexture;
SamplerState baseSampler;

PS_OUTPUT PSMain(PS_INPUT input)
{
    PS_OUTPUT output;
    output.pixelColor = baseTexture.Sample(baseSampler, input.vertexTCoord);
    return output;
}
