// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEngine.h>
#if defined(__LINUX__)
#include <Graphics/GL4/GteGLSLProgramFactory.h>
#include <Graphics/GL4/GLX/GteGLXEngine.h>
#endif
using namespace gte;

union Float
{
    float number;
    unsigned int encoding;
};

union Double
{
    double number;
    uint64_t encoding;
};

template <typename Real, typename Binary>
class TestSubnormals
{
public:
    TestSubnormals(std::string const& filename, std::string const& realname, Binary& result)
    {
#if defined(GTE_DEV_OPENGL)
#if defined(__MSWINDOWS__)
        WGLEngine engine(false);
#else
        GLXEngine engine(false);
#endif
        GLSLProgramFactory factory;
#else
        DX11Engine engine;
        HLSLProgramFactory factory;
#endif

        std::shared_ptr<StructuredBuffer> inputBuffer = std::make_shared<StructuredBuffer>(2, sizeof(Real));
        Real* input = inputBuffer->Get<Real>();
        Binary v0, v1;
        v0.encoding = 1;
        v1.encoding = 1;
        input[0] = v0.number;  // Smallest positive subnormal.
        input[1] = v1.number;  // Same as v0.

        // Compute v0+v1 and store in this buffer.
        std::shared_ptr<StructuredBuffer> outputBuffer = std::make_shared<StructuredBuffer>(1, sizeof(Real));
        outputBuffer->SetUsage(Resource::SHADER_OUTPUT);
        outputBuffer->SetCopyType(Resource::COPY_STAGING_TO_CPU);
        Real* output = outputBuffer->Get<Real>();
        output[0] = (Real)0;

        factory.defines.Set("REAL", realname);
        std::shared_ptr<ComputeProgram> cprogram = factory.CreateFromFile(filename);
        if (!cprogram)
        {
            LogError("Cannot load or compile cshader.");
            inputBuffer = nullptr;
            outputBuffer = nullptr;
            return;
        }
        std::shared_ptr<ComputeShader> cshader = cprogram->GetCShader();
        cshader->Set("inBuffer", inputBuffer);
        cshader->Set("outBuffer", outputBuffer);

        engine.Execute(cprogram, 1, 1, 1);
        engine.CopyGpuToCpu(outputBuffer);

        result.number = output[0];

        inputBuffer = nullptr;
        outputBuffer = nullptr;
        cshader = nullptr;
    }
};

int main(int, char const*[])
{
#if defined(_DEBUG)
    LogReporter reporter(
        "LogReport.txt",
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL,
        Logger::Listener::LISTEN_FOR_ALL);
#endif

    Environment env;
    std::string gtpath = env.GetVariable("GTE_PATH");
    if (gtpath == "")
    {
        LogError("You must create the environment variable GTE_PATH.");
        return -1;
    }
    gtpath += "/Samples/Basics/IEEEFloatingPoint/Shaders/";
#if defined(GTE_DEV_OPENGL)
    gtpath += "TestSubnormals.glsl";
#else
    gtpath += "TestSubnormals.hlsl";
#endif

    // With IEEE 754-2008 behavior that preserves subnormals, the output
    // fresult should have encoding 2 (number is 2^{-148}).  Instead
    // fresult.encoding = 0, which means that the GPU has flushed the
    // subnormal result to zero.
    Float fresult;
    TestSubnormals<float,Float> ftest(gtpath, "float", fresult);

    // With IEEE 754-2008 behavior that preserves subnormals, the output
    // dresult should have encoding 2 (number is 2^{-1073}).  Indeed,
    // dresult.encoding = 2.
    Double dresult;
    TestSubnormals<double,Double> dtest(gtpath, "double", dresult);

    return 0;
}
