/*
 *  fileoutput.h
 *
 *  Header file containing the declaration of the fileoutput functions
 *  The full definition of the functions are found in the associated C++ file
 *  "fileoutput.cpp"
 *
 *  James Bradford
 *  University of Sheffield
 *  Jan 2017
 */


#ifndef FILEOUTPUT_H
#define FILEOUTPUT_H

std::string newDirectory(bool test);
std::ofstream openResults(std::string directory);
std::ofstream openLengths(std::string directory);
std::ofstream openLog(std::string directory);
std::ofstream openMono(std::string directory);
std::ofstream openRegions(std::string directory);
std::ofstream openConcMap(std::string directory);
std::ofstream openLatAConcMap(std::string directory);
std::ofstream openLatActinConcMap(std::string directory);
std::ofstream openArpConcMap(std::string directory);
std::ofstream openMemLengths(std::string directory);
std::ofstream openAngles(std::string directory);
std::ofstream openEngulfment(std::string directory);
std::ofstream openNumContactTarget(std::string directory);
std::ofstream openNumContactMem(std::string directory);
std::ofstream openFusionTime(std::string directory);
void copyConfig(std::string directory, std::ifstream &input_file);



#endif
