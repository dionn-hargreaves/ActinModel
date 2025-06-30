/*
 *  branching.h
 *
 *
 *  James Bradford
 *  University of Sheffield
 *  May 2017
 */

#ifndef BRANCHING_H
#define BRANCHING_H

#include <vector>
#include "Actin.h"
#include "ProteinRegion.h"


void deBranch(std::vector<Actin> &actinVec, const double p_detach,
              int &nActin, bool arpPool, ArpGrid &arpGrid,
              const bool steric, StericGrid &stericGrid);

void deBranch2(Actin &actin, std::vector<Actin> &actinVec, int &nActin, bool arpPool,
               ArpGrid &arpGrid, const bool steric, StericGrid &stericGrid);

void branchingRoutine(std::vector<Actin> &actinVec, const double k_branch,
                      const double tStep,
                      const double arpConc, const double gActinConc,
                      int &nActin,
                      const int nActinLimit, const double currTime,
                      const std::vector<ExcZone> &excZones,
                      const std::vector<Membrane> &membranes,
                      const std::vector<MembraneWall> &memWalls,
                      Cortex &cortex,
                      const bool steric, StericGrid &stericGrid,
                      const bool tether,
                      std::normal_distribution<double> tetherDistri);

void branchRegionsRoutine(std::vector<Actin> &actinVec, const double k_branch,
                          const double tStep,
                          const double gActinConc,
                          int &nActin,
                          const int nActinLimit, const double currTime,
                          std::vector<ProteinRegion> &branchingRegions,
                          const std::vector<ExcZone> &excZones,
                          const std::vector<Membrane> &membranes,
                          const std::vector<MembraneWall> &memWalls,
                          Cortex &cortex,
                          const bool steric, StericGrid &stericGrid);
/*
void branchRegGridRoutine(std::vector<Actin> &actinVec, const double k_branch,
                      const double arpConc, const double gActinConc,
                      bool &branching, int &nActin,
                      const int nActinLimit, const double currTime,
                      std::vector<ProteinRegion> &branchingRegions,
                      const std::vector<ExcZone> &excZones,
                      const bool steric, GactinGrid &gActinGrid);
*/

void branchArpGridRoutine(std::vector<Actin> &actinVec, const double k_branch,
                          const double tStep,
                          const double gActinConc,
                          int &nActin,
                          const int nActinLimit, const double currTime,
                          ArpGrid &arpGrid,
                          const std::vector<ExcZone> &excZones,
                          const std::vector<Membrane> &membranes,
                          const std::vector<MembraneWall> &memWalls,
                          Cortex &cortex,
                          const bool steric, StericGrid &stericGrid);

void branchArpAndActinGridRoutine(std::vector<Actin> &actinVec, const double k_branch,
                                  const double tStep,
                                  int &nActin,
                                  const int nActinLimit, const double currTime,
                                  ArpGrid &arpGrid,GactinGrid &gActinGrid,
                                  const std::vector<ExcZone> &excZones,
                                  const std::vector<Membrane> &membranes,
                                  const std::vector<MembraneWall> &memWalls,
                                  Cortex &cortex,
                                  const bool steric, StericGrid &stericGrid);
#endif
