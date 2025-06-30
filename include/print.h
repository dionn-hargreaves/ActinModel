/*
 *  print.h
 *
 *  Header file containing the declaration of the print functions
 *  The full definition of the functions are found in the associated C++ file
 *  "print.cpp"
 *
 *  James Bradford
 *  University of Sheffield
 *  Jan 2017
 */

#ifndef PRINT_H
#define PRINT_H
#include "ExcZone.h"
#include "ProteinRegion.h"
#include "MembraneWall.h"
#include "Membrane.h"
#include "GactinGrid.h"
#include "ArpGrid.h"

void printLog(std::ofstream &outf, double dt, std::uint32_t origSeed,
              double fracErr);
void printRunTime(std::ofstream &outf, double runTime);
void printActinframes(std::ofstream &outf, int nActin,
                      const std::vector<Actin> &actinVec, double currTime,
                      double dt, const std::vector<ExcZone> &excZones,
                      const std::vector<ProteinRegion> &nucRegions,
                      const std::vector<ProteinRegion> &branchRegions,
                      const std::vector<ProteinRegion> &capRegions,
                      const std::vector<ProteinRegion> &antiCapRegions,
                      const std::vector<ProteinRegion> &sevRegions,
                      const std::vector<MembraneWall> &memWalls,
                      const std::vector<Membrane> &membranes,
                      bool tether, bool crossLinking, const Cortex &cortex);

void printActinHeaders(std::ofstream &outf, double dt_bw_f);

void printGActinGridHeader(std::ofstream &outf, double dt_bw_f, const GactinGrid &gActinGrid);

void printGActinGrid(std::ofstream &outf, const GactinGrid &gActinGrid,
                     double currTime);

void printLatAGrid(std::ofstream &outf, const GactinGrid &gActinGrid,
                     double currTime);

void printLatActinGrid(std::ofstream &outf, const GactinGrid &gActinGrid,
                  double currTime);

void printArpGridHeader(std::ofstream &outf, double dt_bw_f, const ArpGrid &arpGrid);

void printArpGrid(std::ofstream &outf, const ArpGrid &arpGrid,
                  double currTime);

void printMemLenHeader(std::ofstream &outf, double dt_bw_f);

void printMemLen(std::ofstream &outf, const std::vector<Membrane> &membranes, double currTime);

void printAnglesHeader(std::ofstream &outf, double dt_bw_f);

void printAngles(std::ofstream &outf, std::vector<Actin> &actinvec);
void printEngulfment(std::ofstream &outf, double currTime, double engulfment);
void printContact(std::ofstream &outf, double currTime, int numContact);
void printLength(std::ofstream &outf, std::vector<Actin> &actinvec, double currTime);
void printMonoTimeFrames(std::ofstream &outf, int nActin,
                         const std::vector<Actin> &actinVec, double currTime,
                         double dt);
#endif
