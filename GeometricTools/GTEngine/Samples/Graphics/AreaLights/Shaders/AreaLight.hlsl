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

struct VS_INPUT
{
    float3 modelPosition : POSITION;
    float2 modelTCoord : TEXCOORD0;
};

struct VS_OUTPUT
{
    float4 vertexPosition : TEXCOORD0;
    float2 vertexTCoord : TEXCOORD1;
    float4 clipPosition : SV_POSITION;
};

VS_OUTPUT VSMain(VS_INPUT input)
{
    VS_OUTPUT output;
    output.vertexPosition = float4(input.modelPosition, 1.0f);
    output.vertexTCoord = input.modelTCoord;
#if GTE_USE_MAT_VEC
    output.clipPosition = mul(pvwMatrix, output.vertexPosition);
#else
    output.clipPosition = mul(output.vertexPosition, pvwMatrix);
#endif
    return output;
}

cbuffer Material
{
    // Properties of the surface.
    float4 materialEmissive;
    float4 materialAmbient;
    float4 materialDiffuse;
    float4 materialSpecular;
};

cbuffer Camera
{
    float4 cameraModelPosition;
};

cbuffer AreaLight
{
    float4 lightingAmbient;
    float4 lightingDiffuse;
    float4 lightingSpecular;
    float4 lightingAttenuation;

    // The light is rectangular in shape, represented by an oriented bounding
    // rectangle in 3D.  The application is responsible for setting the
    // rectangle members in the model space of the surface to be illuminated.
    float4 rectPosition;  // (x,y,z,1)
    float4 rectNormal, rectAxis0, rectAxis1;  // (x,y,z,0)
    float4 rectExtent;  // (extent0, extent1, *, *)
};

Texture2D<float4> baseTexture;
Texture2D<float4> normalTexture;
SamplerState commonSampler;

struct PS_INPUT
{
    float4 vertexPosition : TEXCOORD0;
    float2 vertexTCoord : TEXCOORD1;
};

struct PS_OUTPUT
{
    float4 pixelColor : SV_TARGET0;
};

float4 ComputeLightModelPosition(in float4 vertexPosition)
{
    float4 diff = vertexPosition - rectPosition;
    float2 crd = float2(dot(rectAxis0, diff), dot(rectAxis1, diff));
    crd = clamp(crd, -rectExtent.xy, rectExtent.xy);
    float4 closest = rectPosition + crd.x * rectAxis0 + crd.y * rectAxis1;
    return closest;
}

PS_OUTPUT PSMain(PS_INPUT input)
{
    PS_OUTPUT output;

    float4 lightModelPosition = ComputeLightModelPosition(input.vertexPosition);
    float3 modelLightDiff = input.vertexPosition.xyz - lightModelPosition.xyz;
    float distance = length(modelLightDiff);
    float attenuation = lightingAttenuation.w / (lightingAttenuation.x + distance *
        (lightingAttenuation.y + distance * lightingAttenuation.z));

    // This code is specific to the brick texture used by the application.  The
    // brick geometry is in a plane whose normal is consistent with the normal
    // map generated for the brick texture.  Thus, the normal may be looked up
    // directly without computing tangent-space quantities.
    float3 vertexNormal = normalTexture.Sample(commonSampler, input.vertexTCoord).rgb;
    vertexNormal = normalize(2.0f * vertexNormal - 1.0f);

    float3 vertexDirection = normalize(modelLightDiff);
    float NDotL = -dot(vertexNormal, vertexDirection);
    float3 viewVector = normalize(cameraModelPosition.xyz - input.vertexPosition.xyz);
    float3 halfVector = normalize(viewVector - vertexDirection);
    float NDotH = dot(vertexNormal, halfVector);
    float4 lighting = lit(NDotL, NDotH, materialSpecular.a);

    float3 emissive = materialEmissive.rgb;
    float3 ambient = materialAmbient.rgb * lightingAmbient.rgb;
    float4 textureDiffuse = baseTexture.Sample(commonSampler, input.vertexTCoord);
    float3 diffuse = materialDiffuse.rgb * textureDiffuse.rgb * lightingDiffuse.rgb;
    float3 specular = materialSpecular.rgb * lightingSpecular.rgb;

    float3 colorRGB = emissive +
        attenuation * (ambient + lighting.y * diffuse + lighting.z * specular);
    float colorA = materialDiffuse.a * textureDiffuse.a;
    output.pixelColor = float4(colorRGB, colorA);

    return output;
}
