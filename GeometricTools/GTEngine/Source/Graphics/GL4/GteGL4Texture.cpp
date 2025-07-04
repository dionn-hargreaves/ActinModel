// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <Graphics/GL4/GteGL4Texture.h>
using namespace gte;

GL4Texture::GL4Texture(Texture const* texture, GLenum target, GLenum targetBinding)
    :
    GL4Resource(texture),
    mTarget(target),
    mTargetBinding(targetBinding),
    mNumLevels(texture->GetNumLevels()),
    mInternalFormat(msGLTextureInternalFormat[texture->GetFormat()]),
    mExternalFormat(msGLTextureExternalFormat[texture->GetFormat()]),
    mExternalType(msGLTextureExternalType[DataFormat::GetChannelType(texture->GetFormat())])
{
}

GLuint const GL4Texture::msGLTextureInternalFormat[DF_NUM_FORMATS]
{
    0,  // DF_UNKNOWN
    0,  // DF_R32G32B32A32_TYPELESS
    GL_RGBA32F,  // DF_R32G32B32A32_FLOAT
    GL_RGBA32UI,  // DF_R32G32B32A32_UINT
    GL_RGBA32I,  // DF_R32G32B32A32_SINT
    0,  // DF_R32G32B32_TYPELESS
    GL_RGB32F,  // DF_R32G32B32_FLOAT
    GL_RGB32UI,  // DF_R32G32B32_UINT
    GL_RGB32I,  // DF_R32G32B32_SINT
    0,  // DF_R16G16B16A16_TYPELESS
    GL_RGBA16F,  // DF_R16G16B16A16_FLOAT
    GL_RGBA16,  // DF_R16G16B16A16_UNORM
    GL_RGBA16UI,  // DF_R16G16B16A16_UINT
    GL_RGBA16_SNORM,  // DF_R16G16B16A16_SNORM
    GL_RGBA16I,  // DF_R16G16B16A16_SINT
    0,  // DF_R32G32_TYPELESS
    GL_RG32F,  // DF_R32G32_FLOAT
    GL_RG32UI,  // DF_R32G32_UINT
    GL_RG32I,  // DF_R32G32_SINT
    0,  // DF_R32G8X24_TYPELESS
    0,  // DF_D32_FLOAT_S8X24_UINT
    0,  // DF_R32_FLOAT_X8X24_TYPELESS
    0,  // DF_X32_TYPELESS_G8X24_UINT
    0,  // DF_R10G10B10A2_TYPELESS
    0,  // DF_R10G10B10A2_UNORM
    0,  // DF_R10G10B10A2_UINT
    GL_R11F_G11F_B10F,  // DF_R11G11B10_FLOAT
    0,  // DF_R8G8B8A8_TYPELESS
    GL_RGBA8,  // DF_R8G8B8A8_UNORM
    GL_RGBA8,  // DF_R8G8B8A8_UNORM_SRGB
    GL_RGBA8UI,  // DF_R8G8B8A8_UINT
    GL_RGBA8_SNORM,  // DF_R8G8B8A8_SNORM
    GL_RGBA8I,  // DF_R8G8B8A8_SINT
    0,  // DF_R16G16_TYPELESS
    GL_RG16F,  // DF_R16G16_FLOAT
    GL_RG16,  // DF_R16G16_UNORM
    GL_RG16UI,  // DF_R16G16_UINT
    GL_R16_SNORM,  // DF_R16G16_SNORM
    GL_R16I,  // DF_R16G16_SINT
    0,  // DF_R32_TYPELESS
    GL_DEPTH_COMPONENT32F,  // DF_D32_FLOAT
    GL_R32F,  // DF_R32_FLOAT
    GL_R32UI,  // DF_R32_UINT
    GL_R32I,  // DF_R32_SINT
    0,  // DF_R24G8_TYPELESS
    GL_DEPTH24_STENCIL8,  // DF_D24_UNORM_S8_UINT
    0,  // DF_R24_UNORM_X8_TYPELESS
    0,  // DF_X24_TYPELESS_G8_UINT
    0,  // DF_R8G8_TYPELESS
    GL_RG8,  // DF_R8G8_UNORM
    GL_RG8UI,  // DF_R8G8_UINT
    GL_RG8_SNORM,  // DF_R8G8_SNORM
    GL_RG8I,  // DF_R8G8_SINT
    0,  // DF_R16_TYPELESS
    GL_R16F,  // DF_R16_FLOAT
    GL_DEPTH_COMPONENT16,  // DF_D16_UNORM
    GL_R16,  // DF_R16_UNORM
    GL_R16UI,  // DF_R16_UINT
    GL_R16_SNORM,  // DF_R16_SNORM
    GL_R16I,  // DF_R16_SINT
    0,  // DF_R8_TYPELESS
    GL_R8,  // DF_R8_UNORM
    GL_R8UI,  // DF_R8_UINT
    GL_R8_SNORM,  // DF_R8_SNORM
    GL_R8I,  // DF_R8_SINT
    0,  // DF_A8_UNORM
    0,  // DF_R1_UNORM
    GL_RGB9_E5,  // DF_R9G9B9E5_SHAREDEXP
    0,  // DF_R8G8_B8G8_UNORM
    0,  // DF_G8R8_G8B8_UNORM
    0,  // DF_BC1_TYPELESS
    0,  // DF_BC1_UNORM
    0,  // DF_BC1_UNORM_SRGB
    0,  // DF_BC2_TYPELESS
    0,  // DF_BC2_UNORM
    0,  // DF_BC2_UNORM_SRGB
    0,  // DF_BC3_TYPELESS
    0,  // DF_BC3_UNORM
    0,  // DF_BC3_UNORM_SRGB
    0,  // DF_BC4_TYPELESS
    0,  // DF_BC4_UNORM
    0,  // DF_BC4_SNORM
    0,  // DF_BC5_TYPELESS
    0,  // DF_BC5_UNORM
    0,  // DF_BC5_SNORM
    GL_RGB565,  // DF_B5G6R5_UNORM
    GL_RGB5_A1,  // DF_B5G5R5A1_UNORM
    GL_RGBA8,  // DF_B8G8R8A8_UNORM
    GL_RGBA8,  // DF_B8G8R8X8_UNORM
    GL_RGB10_A2,  // DF_R10G10B10_XR_BIAS_A2_UNORM
    0,  // DF_B8G8R8A8_TYPELESS
    GL_RGBA8,  // DF_B8G8R8A8_UNORM_SRGB
    0,  // DF_B8G8R8X8_TYPELESS
    GL_RGBA8,  // DF_B8G8R8X8_UNORM_SRGB
    0,  // DF_BC6H_TYPELESS
    0,  // DF_BC6H_UF16
    0,  // DF_BC6H_SF16
    0,  // DF_BC7_TYPELESS
    0,  // DF_BC7_UNORM
    0,  // DF_BC7_UNORM_SRGB
    // DX11.1 formats (TODO: Determine number of channels)
    0,  // DF_AYUV
    0,  // DF_Y410
    0,  // DF_Y416
    0,  // DF_NV12
    0,  // DF_P010
    0,  // DF_P016
    0,  // DF_420_OPAQUE
    0,  // DF_YUY2
    0,  // DF_Y210
    0,  // DF_Y216
    0,  // DF_NV11
    0,  // DF_AI44
    0,  // DF_IA44
    0,  // DF_P8
    0,  // DF_A8P8
    0   // DF_B4G4R4A4_UNORM
};

GLuint const GL4Texture::msGLTextureExternalFormat[DF_NUM_FORMATS]
{
    0,  // DF_UNKNOWN
    GL_RGBA,  // DF_R32G32B32A32_TYPELESS
    GL_RGBA,  // DF_R32G32B32A32_FLOAT
    GL_RGBA_INTEGER,  // DF_R32G32B32A32_UINT
    GL_RGBA_INTEGER,  // DF_R32G32B32A32_SINT
    GL_RGB,  // DF_R32G32B32_TYPELESS
    GL_RGB,  // DF_R32G32B32_FLOAT
    GL_RGB_INTEGER,  // DF_R32G32B32_UINT
    GL_RGB_INTEGER,  // DF_R32G32B32_SINT
    GL_RGBA,  // DF_R16G16B16A16_TYPELESS
    GL_RGBA,  // DF_R16G16B16A16_FLOAT
    GL_RGBA,  // DF_R16G16B16A16_UNORM
    GL_RGBA_INTEGER,  // DF_R16G16B16A16_UINT
    GL_RGBA,  // DF_R16G16B16A16_SNORM
    GL_RGBA_INTEGER,  // DF_R16G16B16A16_SINT
    GL_RG,  // DF_R32G32_TYPELESS
    GL_RG,  // DF_R32G32_FLOAT
    GL_RG_INTEGER,  // DF_R32G32_UINT
    GL_RG_INTEGER,  // DF_R32G32_SINT
    0,  // DF_R32G8X24_TYPELESS
    0,  // DF_D32_FLOAT_S8X24_UINT
    0,  // DF_R32_FLOAT_X8X24_TYPELESS
    0,  // DF_X32_TYPELESS_G8X24_UINT
    GL_RGBA,  // DF_R10G10B10A2_TYPELESS
    GL_RGBA,  // DF_R10G10B10A2_UNORM
    GL_RGBA_INTEGER,  // DF_R10G10B10A2_UINT
    GL_RGB,  // DF_R11G11B10_FLOAT
    GL_RGBA,  // DF_R8G8B8A8_TYPELESS
    GL_RGBA,  // DF_R8G8B8A8_UNORM
    GL_RGBA,  // DF_R8G8B8A8_UNORM_SRGB
    GL_RGBA_INTEGER,  // DF_R8G8B8A8_UINT
    GL_RGBA,  // DF_R8G8B8A8_SNORM
    GL_RGBA_INTEGER,  // DF_R8G8B8A8_SINT
    GL_RG,  // DF_R16G16_TYPELESS
    GL_RG,  // DF_R16G16_FLOAT
    GL_RG,  // DF_R16G16_UNORM
    GL_RG_INTEGER,  // DF_R16G16_UINT
    GL_RG,  // DF_R16G16_SNORM
    GL_RG_INTEGER,  // DF_R16G16_SINT
    GL_RED,  // DF_R32_TYPELESS
    GL_DEPTH_COMPONENT,  // DF_D32_FLOAT
    GL_RED,  // DF_R32_FLOAT
    GL_RED_INTEGER,  // DF_R32_UINT
    GL_RED_INTEGER,  // DF_R32_SINT
    GL_RG,  // DF_R24G8_TYPELESS
    GL_DEPTH_COMPONENT,  // DF_D24_UNORM_S8_UINT
    0,  // DF_R24_UNORM_X8_TYPELESS
    0,  // DF_X24_TYPELESS_G8_UINT
    GL_RG,  // DF_R8G8_TYPELESS
    GL_RG,  // DF_R8G8_UNORM
    GL_RG_INTEGER,  // DF_R8G8_UINT
    GL_RG,  // DF_R8G8_SNORM
    GL_RG_INTEGER,  // DF_R8G8_SINT
    GL_RED,  // DF_R16_TYPELESS
    GL_RED,  // DF_R16_FLOAT
    GL_DEPTH_COMPONENT,  // DF_D16_UNORM
    GL_RED,  // DF_R16_UNORM
    GL_RED_INTEGER,  // DF_R16_UINT
    GL_RED,  // DF_R16_SNORM
    GL_RED_INTEGER,  // DF_R16_SINT
    GL_RED,  // DF_R8_TYPELESS
    GL_RED,  // DF_R8_UNORM
    GL_RED_INTEGER,  // DF_R8_UINT
    GL_RED,  // DF_R8_SNORM
    GL_RED_INTEGER,  // DF_R8_SINT
    0,  // DF_A8_UNORM
    0,  // DF_R1_UNORM
    0,  // DF_R9G9B9E5_SHAREDEXP
    0,  // DF_R8G8_B8G8_UNORM
    0,  // DF_G8R8_G8B8_UNORM
    0,  // DF_BC1_TYPELESS
    0,  // DF_BC1_UNORM
    0,  // DF_BC1_UNORM_SRGB
    0,  // DF_BC2_TYPELESS
    0,  // DF_BC2_UNORM
    0,  // DF_BC2_UNORM_SRGB
    0,  // DF_BC3_TYPELESS
    0,  // DF_BC3_UNORM
    0,  // DF_BC3_UNORM_SRGB
    0,  // DF_BC4_TYPELESS
    0,  // DF_BC4_UNORM
    0,  // DF_BC4_SNORM
    0,  // DF_BC5_TYPELESS
    0,  // DF_BC5_UNORM
    0,  // DF_BC5_SNORM
    GL_BGR,  // DF_B5G6R5_UNORM
    GL_BGRA,  // DF_B5G5R5A1_UNORM
    GL_BGRA,  // DF_B8G8R8A8_UNORM
    GL_BGRA,  // DF_B8G8R8X8_UNORM
    0,  // DF_R10G10B10_XR_BIAS_A2_UNORM
    GL_BGRA,  // DF_B8G8R8A8_TYPELESS
    GL_BGRA,  // DF_B8G8R8A8_UNORM_SRGB
    GL_BGRA,  // DF_B8G8R8X8_TYPELESS
    GL_BGRA,  // DF_B8G8R8X8_UNORM_SRGB
    0,  // DF_BC6H_TYPELESS
    0,  // DF_BC6H_UF16
    0,  // DF_BC6H_SF16
    0,  // DF_BC7_TYPELESS
    0,  // DF_BC7_UNORM
    0,  // DF_BC7_UNORM_SRGB
    // DX11.1 formats (TODO: Determine number of channels)
    0,  // DF_AYUV
    0,  // DF_Y410
    0,  // DF_Y416
    0,  // DF_NV12
    0,  // DF_P010
    0,  // DF_P016
    0,  // DF_420_OPAQUE
    0,  // DF_YUY2
    0,  // DF_Y210
    0,  // DF_Y216
    0,  // DF_NV11
    0,  // DF_AI44
    0,  // DF_IA44
    0,  // DF_P8
    0,  // DF_A8P8
    0   // DF_B4G4R4A4_UNORM
};

GLuint const GL4Texture::msGLTextureExternalType[DF_NUM_CHANNEL_TYPES] =
{
    GL_ZERO,                        // DF_UNSUPPORTED
    GL_BYTE,                        // DF_BYTE
    GL_UNSIGNED_BYTE,               // DF_UBYTE
    GL_SHORT,                       // DF_SHORT
    GL_UNSIGNED_SHORT,              // DF_USHORT
    GL_INT,                         // DF_INT
    GL_UNSIGNED_INT,                // DF_UINT
    GL_HALF_FLOAT,                  // DF_HALF_FLOAT
    GL_FLOAT,                       // DF_FLOAT
    GL_DOUBLE,                      // DF_DOUBLE
    GL_INT_2_10_10_10_REV,          // DF_INT_10_10_2
    GL_UNSIGNED_INT_2_10_10_10_REV, // DF_UINT_10_10_2
    GL_UNSIGNED_INT_10F_11F_11F_REV // DF_FLOAT_11_11_10
};
