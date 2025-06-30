/*
 *  fileoutput.cpp
 *
 *  Quick function file I wrote to create a new directory if one doesn't exist
 *  using "yyyy-mm-dd-runnumber" system. This will only work on Linux machines.
 *  The forward declarations of these functions are found in the associated
 *  header file "fileoutput.h"
 *
 *  James Bradford
 *  University of Sheffield
 *  March 2020
 */

#include "configHeader.h"



std::string newDirectory(bool test)
{
    /* Function that checks whether the directory exists and loops until it
     * finds a runnumber that doesn't exist. It then creates that directory
     * so that files can be outputted.
     */

    auto t = std::time(nullptr);
    char datetmp[100];
    std::strftime(datetmp, 100, "%Y-%m-%d", std::localtime(&t));
    std::string date (datetmp);

    std::string fulldir("./");
    if (test)
    {
        fulldir += std::string("testOutput/");
    }
    else
    {
        fulldir += std::string("output/");
    }
    fulldir += std::string(date);
    fulldir += std::string("-run");


    struct stat sb;
    int runnum;
    for (runnum = 1; ; ++runnum)
    {
        if (stat(((fulldir+std::to_string(runnum)).c_str()), &sb) == 0 && S_ISDIR(sb.st_mode))
        {
            //std::cout << "The directory exists." << std::endl;
        }
        else
        {
            //std::cout << "The directory does not exist." << std::endl;
            int systemRet = system(("mkdir "+((fulldir+std::to_string(runnum)))).c_str());
            if (systemRet == -1)
            {
              // The system method failed
              std::cout << "System call to create directory failed" << std::endl;
              exit(1);
            }
            else
            {
              std::cout << "Saving any output files to: " << fulldir+std::to_string(runnum);
              std::cout << std::endl;
            }
            break;
        }
    }
    fulldir += std::to_string(runnum);

    return fulldir;
}

std::ofstream openResults(std::string directory)
{
    // Function that creates a results file for printing to

    std::string results_dir = directory + std::string("/results.dat");
    std::ofstream file_results(results_dir);
    if (!file_results)
    {
        std::cout << "Could not open results.dat for writing" << std::endl;
        exit(1);
    }
    return file_results;
}

std::ofstream openLengths(std::string directory)
{
    // Function that creates a results file for printing to

    std::string lengths_dir = directory + std::string("/totLengths.dat");
    std::ofstream file_lengths(lengths_dir);
    if (!file_lengths)
    {
        std::cout << "Could not open totLengths.dat for writing" << std::endl;
        exit(1);
    }
    return file_lengths;
}

std::ofstream openLog(std::string directory)
{
    // Function that creates a parameters file for printing to
    std::string log_dir = directory + std::string("/log.dat");
    std::ofstream file_log(log_dir);
    if (!file_log)
    {
        std::cout << "Could not open log.dat for writing" << std::endl;
        exit(1);
    }
    return file_log;

}

std::ofstream openMono(std::string directory)
{
    // Function that creates a results file for printing to

    std::string mono_dir = directory + std::string("/monoTimes.dat");
    std::ofstream file_mono(mono_dir);
    if (!file_mono)
    {
        std::cout << "Could not open monoTimes.dat for writing" << std::endl;
        exit(1);
    }
    return file_mono;
}


std::ofstream openRegions(std::string directory)
{
    // Function that creates a regions file for printing to
    std::string regions_dir = directory + std::string("/regions.dat");
    std::ofstream file_regions(regions_dir);
    if (!file_regions)
    {
        std::cout << "Could not open regions.dat for writing" << std::endl;
        exit(1);
    }
    return file_regions;
}

std::ofstream openConcMap(std::string directory)
{
    // Function that creates a measurement regions file for printing to
    std::string conc_dir = directory + std::string("/concMap.dat");
    std::ofstream file_conc(conc_dir);
    if (!file_conc)
    {
        std::cout << "Could not open concMap.dat for writing" << std::endl;
        exit(1);
    }
    return file_conc;
}

std::ofstream openLatAConcMap(std::string directory)
{
    // Function that creates a measurement regions file for printing to
    std::string conc_dir = directory + std::string("/latAConcMap.dat");
    std::ofstream file_conc(conc_dir);
    if (!file_conc)
    {
        std::cout << "Could not open concMap.dat for writing" << std::endl;
        exit(1);
    }
    return file_conc;
}

std::ofstream openLatActinConcMap(std::string directory)
{
    // Function that creates a measurement regions file for printing to
    std::string conc_dir = directory + std::string("/latActinConcMap.dat");
    std::ofstream file_conc(conc_dir);
    if (!file_conc)
    {
        std::cout << "Could not open concMap.dat for writing" << std::endl;
        exit(1);
    }
    return file_conc;
}

std::ofstream openArpConcMap(std::string directory)
{
    // Function that creates a measurement regions file for printing to
    std::string conc_dir = directory + std::string("/arpConcMap.dat");
    std::ofstream file_conc(conc_dir);
    if (!file_conc)
    {
        std::cout << "Could not open arpConcMap.dat for writing" << std::endl;
        exit(1);
    }
    return file_conc;
}

std::ofstream openMemLengths(std::string directory)
{
    // Function that creates a measurement regions file for printing to
    std::string dir = directory + std::string("/memLengths.dat");
    std::ofstream file_mem(dir);
    if (!file_mem)
    {
        std::cout << "Could not open memLengths.dat for writing" << std::endl;
        exit(1);
    }
    return file_mem;
}

std::ofstream openAngles(std::string directory)
{
    // Function that creates a measurement regions file for printing to
    std::string angles_dir = directory + std::string("/angles.dat");
    std::ofstream file_angles(angles_dir);
    if (!file_angles)
    {
        std::cout << "Could not open angles.dat for writing" << std::endl;
        exit(1);
    }
    return file_angles;
}

std::ofstream openEngulfment(std::string directory)
{
    // Function that creates a measurement regions file for printing to
    std::string engulf_dir = directory + std::string("/engulfment.dat");
    std::ofstream file_engulf(engulf_dir);
    if (!file_engulf)
    {
        std::cout << "Could not open engulfment.dat for writing" << std::endl;
        exit(1);
    }
    return file_engulf;
}
std::ofstream openNumContactTarget(std::string directory)
{
    // Function that creates a measurement regions file for printing to
    std::string contact_dir = directory + std::string("/contactTar.dat");
    std::ofstream file_contact(contact_dir);
    if (!file_contact)
    {
        std::cout << "Could not open contactTar.dat for writing" << std::endl;
        exit(1);
    }
    return file_contact;
}
std::ofstream openNumContactMem(std::string directory)
{
    // Function that creates a measurement regions file for printing to
    std::string contact_dir = directory + std::string("/contactMem.dat");
    std::ofstream file_contact(contact_dir);
    if (!file_contact)
    {
        std::cout << "Could not open contactMem.dat for writing" << std::endl;
        exit(1);
    }
    return file_contact;
}
std::ofstream openFusionTime(std::string directory)
{
    // Function that creates a measurement regions file for printing to
    std::string fusion_dir = directory + std::string("/fusionTimes.dat");
    std::ofstream file_fusion(fusion_dir);
    if (!file_fusion)
    {
        std::cout << "Could not open fusionTimes.dat for writing" << std::endl;
        exit(1);
    }
    return file_fusion;
}

void copyConfig(std::string directory, std::ifstream &input_file)
{
    // Function that creates a new config file to copy the new one to
    std::string config_dir = directory + std::string("/configInput.cfg");
    std::ofstream file_config(config_dir);
    if (!file_config)
    {
        std::cout << "Could not open config.cfg for writing" << std::endl;
        exit(1);
    }

    input_file.clear();
    input_file.seekg (0, std::ios::beg);
    // While there's still stuff left to read
    while (input_file)
    {
        // read stuff from the file into output file
        std::string strInput;
        std::getline(input_file, strInput);
        file_config << strInput << std::endl;
    }
}
