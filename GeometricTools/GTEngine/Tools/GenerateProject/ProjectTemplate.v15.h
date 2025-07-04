// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.8.0 (2017/03/18)

#pragma once

#include <regex>
#include <string>

class TemplateV15
{
public:
    TemplateV15(std::string const& name, std::string const& gtPath, bool& success);

private:
    bool Create(std::string const& name, std::vector<std::string> const& lines, bool useUT8);
    static std::string GetGuidString();

    // Matching patterns.
    static std::regex const mGPPattern;
    static std::regex const mPNPattern;
    static std::regex const mPGPattern;
    static std::regex const mRQPattern;
    static std::regex const mGTPattern;

    // Replacements for patterns.
    std::string mName;
    std::string mGTPath;
    std::string mProjectGUID;
    std::string mRequiredGUID;
    static std::string const msGTGUID;

    // Templates for solution, project, project filter, and source code.
    static std::vector<std::string> const msSolutionLines;
    static std::vector<std::string> const msProjectLines;
    static std::vector<std::string> const msFilterLines;
    static std::vector<std::string> const msWinHLines;
    static std::vector<std::string> const msWinCLines;
};
