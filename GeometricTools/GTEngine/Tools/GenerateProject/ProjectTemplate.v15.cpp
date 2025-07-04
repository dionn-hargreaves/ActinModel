// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.8.0 (2017/03/18)

#include "ProjectTemplate.v15.h"
#include <cstdio>
#include <cstring>
#include <fstream>
#include <regex>
#include <Rpc.h>

TemplateV15::TemplateV15(std::string const& name, std::string const& gtPath, bool& success)
    :
    mName(name),
    mGTPath(gtPath),
    mProjectGUID(GetGuidString()),
    mRequiredGUID(GetGuidString())
{
    std::string solnName = mName + ".v15.sln";
    std::string projName = mName + ".v15.vcxproj";
    std::string filtName = mName + ".v15.vcxproj.filters";
    std::string winhName = mName + "Window.h";
    std::string wincName = mName + "Window.cpp";

    success =
        Create(solnName, msSolutionLines, true) &&
        Create(projName, msProjectLines, false) &&
        Create(filtName, msFilterLines, false) &&
        Create(winhName, msWinHLines, false) &&
        Create(wincName, msWinCLines, false);
}

bool TemplateV15::Create(std::string const& name, std::vector<std::string> const& lines, bool useUT8)
{
    std::ofstream outFile(name.c_str());
    if (outFile)
    {
        if (useUT8)
        {
            unsigned char ut8tag[3] = { 0xEF, 0xBB, 0xBF };
            for (int j = 0; j < 3; ++j)
            {
                outFile << ut8tag[j];
            }
            outFile << std::endl;
        }
        for (auto const& s : lines)
        {
            std::string line = s;
            line = std::regex_replace(line, mGPPattern, mGTPath);
            line = std::regex_replace(line, mPNPattern, mName);
            line = std::regex_replace(line, mPGPattern, mProjectGUID);
            line = std::regex_replace(line, mRQPattern, mRequiredGUID);
            line = std::regex_replace(line, mGTPattern, msGTGUID);
            outFile << line << std::endl;
        }

        outFile.close();
        return true;
    }
    return false;
}

std::string TemplateV15::GetGuidString()
{
    UUID guid;
    RPC_STATUS status = UuidCreate(&guid);
    if (RPC_S_OK != status)
    {
        return "";
    }

    RPC_CSTR stringUuid = nullptr;
    status = UuidToString(&guid, &stringUuid);
    if (RPC_S_OK != status || nullptr == stringUuid)
    {
        return "";
    }

    std::string stringGuid(reinterpret_cast<char*>(stringUuid));
    RpcStringFree(&stringUuid);
    return stringGuid;
}


std::regex const TemplateV15::mGPPattern("GTPATH");
std::regex const TemplateV15::mPNPattern("PROJECTNAME");
std::regex const TemplateV15::mPGPattern("PROJECTGUID");
std::regex const TemplateV15::mRQPattern("REQUIREDGUID");
std::regex const TemplateV15::mGTPattern("GTGUID");
std::string const TemplateV15::msGTGUID("43A54DE9-1F9B-4BC6-A7BC-C3FD13B2C829");

#if _MSC_VER == 1900
// MSVS 2015 Update 1 was shipped with a bug whereby a warning for static
// class members involving std::vector> is generated incorrectly.  The bug
// fix will appear in Update 2.  For now, turn off the warning.
#pragma warning(push)
#pragma warning(disable : 4592)
#endif
std::vector<std::string> const TemplateV15::msSolutionLines =
{
"Microsoft Visual Studio Solution File, Format Version 12.00",
"# Visual Studio 15",
"VisualStudioVersion = 15.0.26228.9",
"MinimumVisualStudioVersion = 10.0.40219.1",
"Project(\"{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}\") = \"PROJECTNAME.v15\", \"PROJECTNAME.v15.vcxproj\", \"{PROJECTGUID}\"",
"EndProject",
"Project(\"{2150E333-8FDC-42A3-9474-1A3956D46DE8}\") = \"Required\", \"Required\", \"{REQUIREDGUID}\"",
"EndProject",
"Project(\"{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}\") = \"GTEngine.v15\", \"GTPATHGTEngine.v15.vcxproj\", \"{GTGUID}\"",
"EndProject",
"Global",
"	GlobalSection(SolutionConfigurationPlatforms) = preSolution",
"		Debug|x86 = Debug|x86",
"		Debug|x64 = Debug|x64",
"		DebugGL4|x86 = DebugGL4|x86",
"		DebugGL4|x64 = DebugGL4|x64",
"		Release|x86 = Release|x86",
"		Release|x64 = Release|x64",
"		ReleaseGL4|x86 = ReleaseGL4|x86",
"		ReleaseGL4|x64 = ReleaseGL4|x64",
"	EndGlobalSection",
"	GlobalSection(ProjectConfigurationPlatforms) = postSolution",
"		{PROJECTGUID}.Debug|x64.ActiveCfg = Debug|x64",
"		{PROJECTGUID}.Debug|x64.Build.0 = Debug|x64",
"		{PROJECTGUID}.Debug|x86.ActiveCfg = Debug|Win32",
"		{PROJECTGUID}.Debug|x86.Build.0 = Debug|Win32",
"		{PROJECTGUID}.DebugGL4|x64.ActiveCfg = DebugGL4|x64",
"		{PROJECTGUID}.DebugGL4|x64.Build.0 = DebugGL4|x64",
"		{PROJECTGUID}.DebugGL4|x86.ActiveCfg = DebugGL4|Win32",
"		{PROJECTGUID}.DebugGL4|x86.Build.0 = DebugGL4|Win32",
"		{PROJECTGUID}.Release|x64.ActiveCfg = Release|x64",
"		{PROJECTGUID}.Release|x64.Build.0 = Release|x64",
"		{PROJECTGUID}.Release|x86.ActiveCfg = Release|Win32",
"		{PROJECTGUID}.Release|x86.Build.0 = Release|Win32",
"		{PROJECTGUID}.ReleaseGL4|x64.ActiveCfg = ReleaseGL4|x64",
"		{PROJECTGUID}.ReleaseGL4|x64.Build.0 = ReleaseGL4|x64",
"		{PROJECTGUID}.ReleaseGL4|x86.ActiveCfg = ReleaseGL4|Win32",
"		{PROJECTGUID}.ReleaseGL4|x86.Build.0 = ReleaseGL4|Win32",
"		{GTGUID}.Debug|x64.ActiveCfg = Debug|x64",
"		{GTGUID}.Debug|x64.Build.0 = Debug|x64",
"		{GTGUID}.Debug|x86.ActiveCfg = Debug|Win32",
"		{GTGUID}.Debug|x86.Build.0 = Debug|Win32",
"		{GTGUID}.DebugGL4|x64.ActiveCfg = DebugGL4|x64",
"		{GTGUID}.DebugGL4|x64.Build.0 = DebugGL4|x64",
"		{GTGUID}.DebugGL4|x86.ActiveCfg = DebugGL4|Win32",
"		{GTGUID}.DebugGL4|x86.Build.0 = DebugGL4|Win32",
"		{GTGUID}.Release|x64.ActiveCfg = Release|x64",
"		{GTGUID}.Release|x64.Build.0 = Release|x64",
"		{GTGUID}.Release|x86.ActiveCfg = Release|Win32",
"		{GTGUID}.Release|x86.Build.0 = Release|Win32",
"		{GTGUID}.ReleaseGL4|x64.ActiveCfg = ReleaseGL4|x64",
"		{GTGUID}.ReleaseGL4|x64.Build.0 = ReleaseGL4|x64",
"		{GTGUID}.ReleaseGL4|x86.ActiveCfg = ReleaseGL4|Win32",
"		{GTGUID}.ReleaseGL4|x86.Build.0 = ReleaseGL4|Win32",
"	EndGlobalSection",
"	GlobalSection(SolutionProperties) = preSolution",
"		HideSolutionNode = FALSE",
"	EndGlobalSection",
"	GlobalSection(NestedProjects) = preSolution",
"		{GTGUID} = {REQUIREDGUID}",
"	EndGlobalSection",
"EndGlobal"
};

std::vector<std::string> const TemplateV15::msProjectLines =
{
"﻿<?xml version=\"1.0\" encoding=\"utf-8\"?>",
"<Project DefaultTargets=\"Build\" ToolsVersion=\"15.0\" xmlns=\"http://schemas.microsoft.com/developer/msbuild/2003\">",
"  <ItemGroup Label=\"ProjectConfigurations\">",
"    <ProjectConfiguration Include=\"Debug|Win32\">",
"      <Configuration>Debug</Configuration>",
"      <Platform>Win32</Platform>",
"    </ProjectConfiguration>",
"    <ProjectConfiguration Include=\"Debug|x64\">",
"      <Configuration>Debug</Configuration>",
"      <Platform>x64</Platform>",
"    </ProjectConfiguration>",
"    <ProjectConfiguration Include=\"DebugGL4|Win32\">",
"      <Configuration>DebugGL4</Configuration>",
"      <Platform>Win32</Platform>",
"    </ProjectConfiguration>",
"    <ProjectConfiguration Include=\"DebugGL4|x64\">",
"      <Configuration>DebugGL4</Configuration>",
"      <Platform>x64</Platform>",
"    </ProjectConfiguration>",
"    <ProjectConfiguration Include=\"ReleaseGL4|Win32\">",
"      <Configuration>ReleaseGL4</Configuration>",
"      <Platform>Win32</Platform>",
"    </ProjectConfiguration>",
"    <ProjectConfiguration Include=\"ReleaseGL4|x64\">",
"      <Configuration>ReleaseGL4</Configuration>",
"      <Platform>x64</Platform>",
"    </ProjectConfiguration>",
"    <ProjectConfiguration Include=\"Release|Win32\">",
"      <Configuration>Release</Configuration>",
"      <Platform>Win32</Platform>",
"    </ProjectConfiguration>",
"    <ProjectConfiguration Include=\"Release|x64\">",
"      <Configuration>Release</Configuration>",
"      <Platform>x64</Platform>",
"    </ProjectConfiguration>",
"  </ItemGroup>",
"  <PropertyGroup Label=\"Globals\">",
"    <VCProjectVersion>15.0</VCProjectVersion>",
"    <ProjectGuid>{PROJECTGUID}</ProjectGuid>",
"    <Keyword>Win32Proj</Keyword>",
"    <RootNamespace>PROJECTNAME</RootNamespace>",
"    <WindowsTargetPlatformVersion>10.0.14393.0</WindowsTargetPlatformVersion>",
"  </PropertyGroup>",
"  <Import Project=\"$(VCTargetsPath)\\Microsoft.Cpp.Default.props\" />",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='Debug|Win32'\" Label=\"Configuration\">",
"    <ConfigurationType>Application</ConfigurationType>",
"    <UseDebugLibraries>true</UseDebugLibraries>",
"    <PlatformToolset>v141</PlatformToolset>",
"    <CharacterSet>Unicode</CharacterSet>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='DebugGL4|Win32'\" Label=\"Configuration\">",
"    <ConfigurationType>Application</ConfigurationType>",
"    <UseDebugLibraries>true</UseDebugLibraries>",
"    <PlatformToolset>v141</PlatformToolset>",
"    <CharacterSet>Unicode</CharacterSet>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='Debug|x64'\" Label=\"Configuration\">",
"    <ConfigurationType>Application</ConfigurationType>",
"    <UseDebugLibraries>true</UseDebugLibraries>",
"    <PlatformToolset>v141</PlatformToolset>",
"    <CharacterSet>Unicode</CharacterSet>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='DebugGL4|x64'\" Label=\"Configuration\">",
"    <ConfigurationType>Application</ConfigurationType>",
"    <UseDebugLibraries>true</UseDebugLibraries>",
"    <PlatformToolset>v141</PlatformToolset>",
"    <CharacterSet>Unicode</CharacterSet>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='Release|Win32'\" Label=\"Configuration\">",
"    <ConfigurationType>Application</ConfigurationType>",
"    <UseDebugLibraries>false</UseDebugLibraries>",
"    <PlatformToolset>v141</PlatformToolset>",
"    <WholeProgramOptimization>true</WholeProgramOptimization>",
"    <CharacterSet>Unicode</CharacterSet>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='ReleaseGL4|Win32'\" Label=\"Configuration\">",
"    <ConfigurationType>Application</ConfigurationType>",
"    <UseDebugLibraries>false</UseDebugLibraries>",
"    <PlatformToolset>v141</PlatformToolset>",
"    <WholeProgramOptimization>true</WholeProgramOptimization>",
"    <CharacterSet>Unicode</CharacterSet>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='Release|x64'\" Label=\"Configuration\">",
"    <ConfigurationType>Application</ConfigurationType>",
"    <UseDebugLibraries>false</UseDebugLibraries>",
"    <PlatformToolset>v141</PlatformToolset>",
"    <WholeProgramOptimization>true</WholeProgramOptimization>",
"    <CharacterSet>Unicode</CharacterSet>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='ReleaseGL4|x64'\" Label=\"Configuration\">",
"    <ConfigurationType>Application</ConfigurationType>",
"    <UseDebugLibraries>false</UseDebugLibraries>",
"    <PlatformToolset>v141</PlatformToolset>",
"    <WholeProgramOptimization>true</WholeProgramOptimization>",
"    <CharacterSet>Unicode</CharacterSet>",
"  </PropertyGroup>",
"  <Import Project=\"$(VCTargetsPath)\\Microsoft.Cpp.props\" />",
"  <ImportGroup Label=\"ExtensionSettings\">",
"  </ImportGroup>",
"  <ImportGroup Label = \"Shared\">",
"  </ImportGroup>",
"  <ImportGroup Label = \"PropertySheets\" Condition = \"'$(Configuration)|$(Platform)'=='Debug|Win32'\">",
"    <Import Project = \"$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props\" Condition = \"exists('$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props')\" Label = \"LocalAppDataPlatform\" />",
"  </ImportGroup>",
"  <ImportGroup Condition=\"'$(Configuration)|$(Platform)'=='DebugGL4|Win32'\" Label=\"PropertySheets\">",
"    <Import Project=\"$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props\" Condition=\"exists('$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props')\" Label=\"LocalAppDataPlatform\" />",
"  </ImportGroup>",
"  <ImportGroup Condition=\"'$(Configuration)|$(Platform)'=='Debug|x64'\" Label=\"PropertySheets\">",
"    <Import Project=\"$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props\" Condition=\"exists('$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props')\" Label=\"LocalAppDataPlatform\" />",
"  </ImportGroup>",
"  <ImportGroup Condition=\"'$(Configuration)|$(Platform)'=='DebugGL4|x64'\" Label=\"PropertySheets\">",
"    <Import Project=\"$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props\" Condition=\"exists('$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props')\" Label=\"LocalAppDataPlatform\" />",
"  </ImportGroup>",
"  <ImportGroup Label = \"PropertySheets\" Condition = \"'$(Configuration)|$(Platform)'=='Release|Win32'\">",
"    <Import Project = \"$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props\" Condition = \"exists('$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props')\" Label = \"LocalAppDataPlatform\" />",
"  </ImportGroup>",
"  <ImportGroup Condition=\"'$(Configuration)|$(Platform)'=='ReleaseGL4|Win32'\" Label=\"PropertySheets\">",
"    <Import Project=\"$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props\" Condition=\"exists('$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props')\" Label=\"LocalAppDataPlatform\" />",
"  </ImportGroup>",
"  <ImportGroup Condition=\"'$(Configuration)|$(Platform)'=='Release|x64'\" Label=\"PropertySheets\">",
"    <Import Project=\"$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props\" Condition=\"exists('$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props')\" Label=\"LocalAppDataPlatform\" />",
"  </ImportGroup>",
"  <ImportGroup Condition=\"'$(Configuration)|$(Platform)'=='ReleaseGL4|x64'\" Label=\"PropertySheets\">",
"    <Import Project=\"$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props\" Condition=\"exists('$(UserRootDir)\\Microsoft.Cpp.$(Platform).user.props')\" Label=\"LocalAppDataPlatform\" />",
"  </ImportGroup>",
"  <PropertyGroup Label=\"UserMacros\" />",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='Debug|Win32'\">",
"    <LinkIncremental>true</LinkIncremental>",
"    <OutDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</OutDir>",
"    <IntDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</IntDir>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='DebugGL4|Win32'\">",
"    <LinkIncremental>true</LinkIncremental>",
"    <OutDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</OutDir>",
"    <IntDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</IntDir>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='Debug|x64'\">",
"    <LinkIncremental>true</LinkIncremental>",
"    <OutDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</OutDir>",
"    <IntDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</IntDir>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='DebugGL4|x64'\">",
"    <LinkIncremental>true</LinkIncremental>",
"    <OutDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</OutDir>",
"    <IntDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</IntDir>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='Release|Win32'\">",
"    <LinkIncremental>false</LinkIncremental>",
"    <OutDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</OutDir>",
"    <IntDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</IntDir>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='ReleaseGL4|Win32'\">",
"    <LinkIncremental>false</LinkIncremental>",
"    <OutDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</OutDir>",
"    <IntDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</IntDir>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='Release|x64'\">",
"    <LinkIncremental>false</LinkIncremental>",
"    <OutDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</OutDir>",
"    <IntDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</IntDir>",
"  </PropertyGroup>",
"  <PropertyGroup Condition=\"'$(Configuration)|$(Platform)'=='ReleaseGL4|x64'\">",
"    <LinkIncremental>false</LinkIncremental>",
"    <OutDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</OutDir>",
"    <IntDir>_Output\\$(PlatformToolset)\\$(Platform)\\$(Configuration)\\</IntDir>",
"  </PropertyGroup>",
"  <ItemDefinitionGroup Condition=\"'$(Configuration)|$(Platform)'=='Debug|Win32'\">",
"    <ClCompile>",
"      <PrecompiledHeader>",
"      </PrecompiledHeader>",
"      <WarningLevel>Level4</WarningLevel>",
"      <Optimization>Disabled</Optimization>",
"      <PreprocessorDefinitions>WIN32;_DEBUG;_CONSOLE;%(PreprocessorDefinitions)</PreprocessorDefinitions>",
"      <AdditionalIncludeDirectories>GTPATHInclude</AdditionalIncludeDirectories>",
"      <TreatWarningAsError>true</TreatWarningAsError>",
"    </ClCompile>",
"    <Link>",
"      <SubSystem>Console</SubSystem>",
"      <GenerateDebugInformation>true</GenerateDebugInformation>",
"      <AdditionalDependencies>d3d11.lib;d3dcompiler.lib;dxgi.lib;dxguid.lib;Windowscodecs.lib;%(AdditionalDependencies)</AdditionalDependencies>",
"    </Link>",
"  </ItemDefinitionGroup>",
"  <ItemDefinitionGroup Condition=\"'$(Configuration)|$(Platform)'=='DebugGL4|Win32'\">",
"    <ClCompile>",
"      <PrecompiledHeader>",
"      </PrecompiledHeader>",
"      <WarningLevel>Level4</WarningLevel>",
"      <Optimization>Disabled</Optimization>",
"      <PreprocessorDefinitions>WIN32;_DEBUG;_CONSOLE;GTE_DEV_OPENGL;%(PreprocessorDefinitions)</PreprocessorDefinitions>",
"      <AdditionalIncludeDirectories>GTPATHInclude</AdditionalIncludeDirectories>",
"      <TreatWarningAsError>true</TreatWarningAsError>",
"    </ClCompile>",
"    <Link>",
"      <SubSystem>Console</SubSystem>",
"      <GenerateDebugInformation>true</GenerateDebugInformation>",
"      <AdditionalDependencies>d3d11.lib;d3dcompiler.lib;dxgi.lib;dxguid.lib;Windowscodecs.lib;opengl32.lib;%(AdditionalDependencies)</AdditionalDependencies>",
"    </Link>",
"  </ItemDefinitionGroup>",
"  <ItemDefinitionGroup Condition=\"'$(Configuration)|$(Platform)'=='Debug|x64'\">",
"    <ClCompile>",
"      <PrecompiledHeader>",
"      </PrecompiledHeader>",
"      <WarningLevel>Level4</WarningLevel>",
"      <Optimization>Disabled</Optimization>",
"      <PreprocessorDefinitions>WIN32;_DEBUG;_CONSOLE;%(PreprocessorDefinitions)</PreprocessorDefinitions>",
"      <AdditionalIncludeDirectories>GTPATHInclude</AdditionalIncludeDirectories>",
"      <TreatWarningAsError>true</TreatWarningAsError>",
"    </ClCompile>",
"    <Link>",
"      <SubSystem>Console</SubSystem>",
"      <GenerateDebugInformation>true</GenerateDebugInformation>",
"      <AdditionalDependencies>d3d11.lib;d3dcompiler.lib;dxgi.lib;dxguid.lib;Windowscodecs.lib;%(AdditionalDependencies)</AdditionalDependencies>",
"    </Link>",
"  </ItemDefinitionGroup>",
"  <ItemDefinitionGroup Condition=\"'$(Configuration)|$(Platform)'=='DebugGL4|x64'\">",
"    <ClCompile>",
"      <PrecompiledHeader>",
"      </PrecompiledHeader>",
"      <WarningLevel>Level4</WarningLevel>",
"      <Optimization>Disabled</Optimization>",
"      <PreprocessorDefinitions>WIN32;_DEBUG;_CONSOLE;GTE_DEV_OPENGL;%(PreprocessorDefinitions)</PreprocessorDefinitions>",
"      <AdditionalIncludeDirectories>GTPATHInclude</AdditionalIncludeDirectories>",
"      <TreatWarningAsError>true</TreatWarningAsError>",
"    </ClCompile>",
"    <Link>",
"      <SubSystem>Console</SubSystem>",
"      <GenerateDebugInformation>true</GenerateDebugInformation>",
"      <AdditionalDependencies>d3d11.lib;d3dcompiler.lib;dxgi.lib;dxguid.lib;Windowscodecs.lib;opengl32.lib;%(AdditionalDependencies)</AdditionalDependencies>",
"    </Link>",
"  </ItemDefinitionGroup>",
"  <ItemDefinitionGroup Condition=\"'$(Configuration)|$(Platform)'=='Release|Win32'\">",
"    <ClCompile>",
"      <WarningLevel>Level4</WarningLevel>",
"      <PrecompiledHeader>",
"      </PrecompiledHeader>",
"      <Optimization>MaxSpeed</Optimization>",
"      <FunctionLevelLinking>true</FunctionLevelLinking>",
"      <IntrinsicFunctions>true</IntrinsicFunctions>",
"      <PreprocessorDefinitions>WIN32;NDEBUG;_CONSOLE;%(PreprocessorDefinitions)</PreprocessorDefinitions>",
"      <AdditionalIncludeDirectories>GTPATHInclude</AdditionalIncludeDirectories>",
"      <TreatWarningAsError>true</TreatWarningAsError>",
"    </ClCompile>",
"    <Link>",
"      <SubSystem>Console</SubSystem>",
"      <GenerateDebugInformation>true</GenerateDebugInformation>",
"      <EnableCOMDATFolding>true</EnableCOMDATFolding>",
"      <OptimizeReferences>true</OptimizeReferences>",
"      <AdditionalDependencies>d3d11.lib;d3dcompiler.lib;dxgi.lib;dxguid.lib;Windowscodecs.lib;%(AdditionalDependencies)</AdditionalDependencies>",
"    </Link>",
"  </ItemDefinitionGroup>",
"  <ItemDefinitionGroup Condition=\"'$(Configuration)|$(Platform)'=='ReleaseGL4|Win32'\">",
"    <ClCompile>",
"      <WarningLevel>Level4</WarningLevel>",
"      <PrecompiledHeader>",
"      </PrecompiledHeader>",
"      <Optimization>MaxSpeed</Optimization>",
"      <FunctionLevelLinking>true</FunctionLevelLinking>",
"      <IntrinsicFunctions>true</IntrinsicFunctions>",
"      <PreprocessorDefinitions>WIN32;NDEBUG;_CONSOLE;GTE_DEV_OPENGL;%(PreprocessorDefinitions)</PreprocessorDefinitions>",
"      <AdditionalIncludeDirectories>GTPATHInclude</AdditionalIncludeDirectories>",
"      <TreatWarningAsError>true</TreatWarningAsError>",
"    </ClCompile>",
"    <Link>",
"      <SubSystem>Console</SubSystem>",
"      <GenerateDebugInformation>true</GenerateDebugInformation>",
"      <EnableCOMDATFolding>true</EnableCOMDATFolding>",
"      <OptimizeReferences>true</OptimizeReferences>",
"      <AdditionalDependencies>d3d11.lib;d3dcompiler.lib;dxgi.lib;dxguid.lib;Windowscodecs.lib;opengl32.lib;%(AdditionalDependencies)</AdditionalDependencies>",
"    </Link>",
"  </ItemDefinitionGroup>",
"  <ItemDefinitionGroup Condition=\"'$(Configuration)|$(Platform)'=='Release|x64'\">",
"    <ClCompile>",
"      <WarningLevel>Level4</WarningLevel>",
"      <PrecompiledHeader>",
"      </PrecompiledHeader>",
"      <Optimization>MaxSpeed</Optimization>",
"      <FunctionLevelLinking>true</FunctionLevelLinking>",
"      <IntrinsicFunctions>true</IntrinsicFunctions>",
"      <PreprocessorDefinitions>WIN32;NDEBUG;_CONSOLE;%(PreprocessorDefinitions)</PreprocessorDefinitions>",
"      <AdditionalIncludeDirectories>GTPATHInclude</AdditionalIncludeDirectories>",
"      <TreatWarningAsError>true</TreatWarningAsError>",
"    </ClCompile>",
"    <Link>",
"      <SubSystem>Console</SubSystem>",
"      <GenerateDebugInformation>true</GenerateDebugInformation>",
"      <EnableCOMDATFolding>true</EnableCOMDATFolding>",
"      <OptimizeReferences>true</OptimizeReferences>",
"      <AdditionalDependencies>d3d11.lib;d3dcompiler.lib;dxgi.lib;dxguid.lib;Windowscodecs.lib;%(AdditionalDependencies)</AdditionalDependencies>",
"    </Link>",
"  </ItemDefinitionGroup>",
"  <ItemDefinitionGroup Condition=\"'$(Configuration)|$(Platform)'=='ReleaseGL4|x64'\">",
"    <ClCompile>",
"      <WarningLevel>Level4</WarningLevel>",
"      <PrecompiledHeader>",
"      </PrecompiledHeader>",
"      <Optimization>MaxSpeed</Optimization>",
"      <FunctionLevelLinking>true</FunctionLevelLinking>",
"      <IntrinsicFunctions>true</IntrinsicFunctions>",
"      <PreprocessorDefinitions>WIN32;NDEBUG;_CONSOLE;GTE_DEV_OPENGL;%(PreprocessorDefinitions)</PreprocessorDefinitions>",
"      <AdditionalIncludeDirectories>GTPATHInclude</AdditionalIncludeDirectories>",
"      <TreatWarningAsError>true</TreatWarningAsError>",
"    </ClCompile>",
"    <Link>",
"      <SubSystem>Console</SubSystem>",
"      <GenerateDebugInformation>true</GenerateDebugInformation>",
"      <EnableCOMDATFolding>true</EnableCOMDATFolding>",
"      <OptimizeReferences>true</OptimizeReferences>",
"      <AdditionalDependencies>d3d11.lib;d3dcompiler.lib;dxgi.lib;dxguid.lib;Windowscodecs.lib;opengl32.lib;%(AdditionalDependencies)</AdditionalDependencies>",
"    </Link>",
"  </ItemDefinitionGroup>",
"  <ItemGroup>",
"    <ClCompile Include=\"PROJECTNAMEWindow.cpp\" />",
"  </ItemGroup>",
"  <ItemGroup>",
"    <ClInclude Include=\"PROJECTNAMEWindow.h\" />",
"  </ItemGroup>",
"  <ItemGroup>",
"    <ProjectReference Include=\"GTPATHGTEngine.v15.vcxproj\">",
"      <Project>{GTGUID}</Project>",
"    </ProjectReference>",
"  </ItemGroup>",
"  <Import Project=\"$(VCTargetsPath)\\Microsoft.Cpp.targets\" />",
"  <ImportGroup Label=\"ExtensionTargets\">",
"  </ImportGroup>",
"</Project>"
};

std::vector<std::string> const TemplateV15::msFilterLines =
{
"﻿<?xml version=\"1.0\" encoding=\"utf-8\"?>",
"<Project ToolsVersion=\"4.0\" xmlns=\"http://schemas.microsoft.com/developer/msbuild/2003\">",
"  <ItemGroup>",
"    <Filter Include=\"Source Files\">",
"      <UniqueIdentifier>{4FC737F1-C7A5-4376-A066-2A32D752A2FF}</UniqueIdentifier>",
"      <Extensions>cpp;c;cc;cxx;def;odl;idl;hpj;bat;asm;asmx</Extensions>",
"    </Filter>",
"    <Filter Include=\"Header Files\">",
"      <UniqueIdentifier>{93995380-89BD-4b04-88EB-625FBE52EBFB}</UniqueIdentifier>",
"      <Extensions>h;hh;hpp;hxx;hm;inl;inc;xsd</Extensions>",
"    </Filter>",
"  </ItemGroup>",
"  <ItemGroup>",
"    <ClCompile Include=\"PROJECTNAMEWindow.cpp\">",
"      <Filter>Source Files</Filter>",
"    </ClCompile>",
"  </ItemGroup>",
"  <ItemGroup>",
"    <ClInclude Include=\"PROJECTNAMEWindow.h\">",
"      <Filter>Header Files</Filter>",
"    </ClInclude>",
"  </ItemGroup>",
"</Project>"
};

std::vector<std::string> const TemplateV15::msWinHLines =
{
"#pragma once",
"",
"#include <GTEngine.h>",
"using namespace gte;",
"",
"class PROJECTNAMEWindow : public Window3",
"{",
"public:",
"    PROJECTNAMEWindow(Parameters& parameters);",
"",
"    virtual void OnIdle() override;",
"",
"private:",
"};"
};

std::vector<std::string> const TemplateV15::msWinCLines =
{
"#include \"PROJECTNAMEWindow.h\"",
"",
"int main(int, char const*[])",
"{",
"#if defined(_DEBUG)",
"    LogReporter reporter(",
"        \"LogReport.txt\",",
"        Logger::Listener::LISTEN_FOR_ALL,",
"        Logger::Listener::LISTEN_FOR_ALL,",
"        Logger::Listener::LISTEN_FOR_ALL,",
"        Logger::Listener::LISTEN_FOR_ALL);",
"#endif",
"",
"    Window::Parameters parameters(L\"PROJECTNAMEWindow\", 0, 0, 512, 512);",
"    auto window = TheWindowSystem.Create<PROJECTNAMEWindow>(parameters);",
"    TheWindowSystem.MessagePump(window, TheWindowSystem.DEFAULT_ACTION);",
"    TheWindowSystem.Destroy<PROJECTNAMEWindow>(window);",
"    return 0;",
"}",
"",
"PROJECTNAMEWindow::PROJECTNAMEWindow(Parameters& parameters)",
"    :",
"    Window3(parameters)",
"{",
"}",
"",
"void PROJECTNAMEWindow::OnIdle()",
"{",
"    mTimer.Measure();",
"",
"    mEngine->ClearBuffers();",
"    mEngine->Draw(8, mYSize - 8, { 0.0f, 0.0f, 0.0f, 1.0f }, mTimer.GetFPS());",
"    mEngine->DisplayColorBuffer(0);",
"",
"    mTimer.UpdateFrameCount();",
"}"
};

#if _MSC_VER == 1900
#pragma warning(pop)
#endif
