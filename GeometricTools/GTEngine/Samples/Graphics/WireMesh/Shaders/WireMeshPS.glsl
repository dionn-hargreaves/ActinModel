// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

uniform WireParameters
{
    vec4 meshColor;
    vec4 edgeColor;
    vec2 windowSize;
};

layout(location = 0) in vec4 pixelColor;
layout(location = 1) noperspective in vec3 edgeDistance;
layout(location = 0) out vec4 finalColor;

void main()
{
    float dmin = min(min(edgeDistance[0], edgeDistance[1]), edgeDistance[2]);
    float blend = smoothstep(0.0f, 1.0f, dmin);
    finalColor = mix(edgeColor, pixelColor, blend);
    finalColor.a = dmin;
}
