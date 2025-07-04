// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

uniform PVWMatrix
{
    mat4 pvwMatrix;
};

layout(location = 0) in vec3 modelPosition;
layout(location = 1) in vec2 modelTCoord;
layout(location = 0) out vec4 vertexPosition;
layout(location = 1) out vec2 vertexTCoord;

void main()
{
    vertexPosition = vec4(modelPosition, 1.0f);
    vertexTCoord = modelTCoord;
#if GTE_USE_MAT_VEC
    gl_Position = pvwMatrix * vertexPosition;
#else
    gl_Position = vertexPosition * pvwMatrix;
#endif
}
