// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <LowLevel/GteLogToStringArray.h>
using namespace gte;


LogToStringArray::LogToStringArray(std::string const& name, int flags)
    :
    Logger::Listener(flags),
    mName(name)
{
}

std::string const& LogToStringArray::GetName() const
{
    return mName;
}

std::vector<std::string> const& LogToStringArray::GetMessages() const
{
    return mMessages;
}

std::vector<std::string>& LogToStringArray::GetMessages()
{
    return mMessages;
}

void LogToStringArray::Report(std::string const& message)
{
    mMessages.push_back(message);
}

