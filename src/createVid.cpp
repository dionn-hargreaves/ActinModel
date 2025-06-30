/*
 *  createVid.cpp
 *
 *
 *  Simple subroutine that calls a bash script to create a video if requested
 *
 *  James Bradford
 *  University of Sheffield
 *  Oct 2018
 */

#include "configHeader.h"
#include "createVid.h"

void callToCreateVid(std::string outdir, double FPS, bool test)
{
        char* USER = getenv("USER"); // get linux username

        std::string scriptPath = "plottingScripts/stackToVid.sh";
        std::string outPath = "/home/"+boost::lexical_cast<std::string>(USER)+"/ActinModelling";
        std::string filename;
        std::string call;
        if (test)
        {
                outdir.erase(0,14);
                filename = outPath + "/testPlotting/testVideos/" + outdir;
                outdir = outPath + "/testPlotting/testFrames/" + outdir;
        }
        else
        {
                outdir.erase(0,10);
                filename = outPath + "/plotting/videos/" + outdir;
                outdir = outPath + "/plotting/frames/" + outdir;
        }

        call = scriptPath + " " + outdir + " " + filename + " " + boost::lexical_cast<std::string>(FPS) + " " + boost::lexical_cast<std::string>(test);
        std::cout << call << std::endl;
        int systemRet = system(call.c_str());
        if (systemRet == -1)
        {
          // The system method failed
          std::cout << "System call to create vid failed" << std::endl;
        }
        else
        {
            std::cout << "Video created" << std::endl;
        }

}
