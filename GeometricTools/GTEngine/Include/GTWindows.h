// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.1 (2016/09/12)

#pragma once

#include <Applications/GteWindowBase.h>

#if defined(__MSWINDOWS__)

#include <Applications/MSW/GteWICFileIO.h>
#include <Applications/MSW/GteMSWWindow.h>
#include <Applications/MSW/GteMSWWindowSystem.h>

#if defined(GTE_DEV_OPENGL)
#include <Applications/MSW/WGL/GteWindow.h>
#include <Applications/MSW/WGL/GteWindowSystem.h>
#elif defined(GTE_USE_DX12)
#include <Applications/MSW/DX12/GteWindow.h>
#include <Applications/MSW/DX12/GteWindowSystem.h>
#else
#include <Applications/MSW/DX11/GteWindow.h>
#include <Applications/MSW/DX11/GteWindowSystem.h>
#endif

#endif

#if defined(__LINUX__)
#include <Applications/GLX/GteWICFileIO.h>
#include <Applications/GLX/GteWindow.h>
#include <Applications/GLX/GteWindowSystem.h>
#endif
