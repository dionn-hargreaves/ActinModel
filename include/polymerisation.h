/*
 *  polymerisation.h
 *
 *
 *  James Bradford
 *  University of Sheffield
 *  May 2017
 */

#ifndef POLYMERISATION_H
#define POLYMERISATION_H
#include "ExcZone.h"
#include "GactinGrid.h"

void polymerisationBarbed(std::vector<Actin> &actinVec, const double gActinConc,
                          const double k_on, const double dt, const double currTime,
                          const std::vector<ExcZone> &excZones,
                          const bool steric,
                          std::vector<MembraneWall> &memWalls,
                          std::vector<Membrane> &membranes,
                          const double temp, const bool tether,
                          std::normal_distribution<double> tetherDistri,
                          std::vector<ProteinRegion> &branchRegions,
                          std::vector<ProteinRegion> &nucRegions,
                          std::vector<ProteinRegion> &capRegions,
                          std::vector<ProteinRegion> &sevRegions,
                          StericGrid &stericGrid, const bool bDynamics,
                          Cortex &cortex, int &nActin, bool arpPool, ArpGrid &arpGrid);


void polymerisationBarbedGrid(std::vector<Actin> &actinVec, GactinGrid &gActinGrid,
                              const double k_on, const double dt, const double currTime,
                              const std::vector<ExcZone> &excZones, const bool steric,
                              std::vector<MembraneWall> &memWalls,
                              std::vector<Membrane> &membranes,
                              const double temp,
                              const bool tether,
                              std::normal_distribution<double> tetherDistri,
                              std::vector<ProteinRegion> &branchRegions,
                              std::vector<ProteinRegion> &nucRegions,
                              std::vector<ProteinRegion> &capRegions,
                              std::vector<ProteinRegion> &sevRegions,
                              StericGrid &stericGrid,
                              const bool bDynamics, Cortex &cortex,
                              int &nActin, bool arpPool, ArpGrid &arpGrid);

void polymerisationPointed(std::vector<Actin> &actinVec, const double gActinConc,
                           const double k_on, const double dt, const double currTime,
                           const std::vector<ExcZone> &excZones, const bool steric,
                           std::vector<MembraneWall> &memWalls,
                           std::vector<Membrane> &membranes,
                           const double temp, const bool tether,
                           std::normal_distribution<double> tetherDistri,
                           std::vector<ProteinRegion> &branchRegions,
                           std::vector<ProteinRegion> &nucRegions,
                           std::vector<ProteinRegion> &capRegions,
                           std::vector<ProteinRegion> &sevRegions,
                           StericGrid &stericGrid, const bool bDynamics,
                           Cortex &cortex,
                           int &nActin, bool arpPool, ArpGrid &arpGrid);

void polymerisationPointedGrid(std::vector<Actin> &actinVec, GactinGrid &gActinGrid,
                               const double k_on, const double dt, const double currTime,
                               const std::vector<ExcZone> &excZones, const bool steric,
                               std::vector<MembraneWall> &memWalls,
                               std::vector<Membrane> &membranes,
                               const double temp,
                               const bool tether,
                               std::normal_distribution<double> tetherDistri,
                               std::vector<ProteinRegion> &branchRegions,
                               std::vector<ProteinRegion> &nucRegions,
                               std::vector<ProteinRegion> &capRegions,
                               std::vector<ProteinRegion> &sevRegions,
                               StericGrid &stericGrid, const bool bDynamics,
                               Cortex &cortex,
                               int &nActin, bool arpPool, ArpGrid &arpGrid);


void depolymerisationBarbed(std::vector<Actin> &actinVec, const double k_off,
                            const double dt, const bool tether, int &nActin,
                            const bool dissociate,
                            bool arpPool, ArpGrid &arpGrid,
                            const bool steric, StericGrid &stericGrid,
                            const std::vector<ExcZone> &excZones,
                            std::vector<MembraneWall> &memWalls,
                            std::vector<Membrane> &membranes,
                            Cortex &cortex,
                            const bool bDynamics,
                            std::vector<ProteinRegion> &branchRegions,
                            std::vector<ProteinRegion> &nucRegions,
                            std::vector<ProteinRegion> &capRegions,
                            std::vector<ProteinRegion> &sevRegions);

void depolymerisationBarbedForce(int i, std::vector<Actin> &actinVec,
                                 const bool tether, int &nActin,
                                 bool arpPool, ArpGrid &arpGrid, const bool steric,
                                 StericGrid &stericGrid,
                                 std::vector<ExcZone> &excZones,
                                 std::vector<MembraneWall> &memWalls,
                                 std::vector<Membrane> &membranes,
                                 Cortex &cortex,
                                 const bool bDynamics,
                                 std::vector<ProteinRegion> &branchRegions,
                                 std::vector<ProteinRegion> &nucRegions,
                                 std::vector<ProteinRegion> &capRegions,
                                 std::vector<ProteinRegion> &sevRegions);

void depolymerisationBarbedGrid(std::vector<Actin> &actinVec, const double k_off,
                                const double dt, const bool tether, int &nActin,
                                const bool dissociate, GactinGrid &gActinGrid,
                                bool arpPool, ArpGrid &arpGrid,
                                const bool steric, StericGrid &stericGrid,
                                std::vector<ExcZone> &excZones,
                                std::vector<MembraneWall> &memWalls,
                                std::vector<Membrane> &membranes,
                                Cortex &cortex,
                                const bool bDynamics,
                                std::vector<ProteinRegion> &branchRegions,
                                std::vector<ProteinRegion> &nucRegions,
                                std::vector<ProteinRegion> &capRegions,
                                std::vector<ProteinRegion> &sevRegions);

void depolymerisationPointed(std::vector<Actin> &actinVec, const double k_off,
                             const double dt, const bool tether, int &nActin,
                             const bool dissociate,
                             bool arpPool, ArpGrid &arpGrid,
                             const bool steric, StericGrid &stericGrid,
                             std::vector<ExcZone> &excZones,
                             std::vector<MembraneWall> &memWalls,
                             std::vector<Membrane> &membranes,
                             Cortex &cortex,
                             const bool bDynamics,
                             std::vector<ProteinRegion> &branchRegions,
                             std::vector<ProteinRegion> &nucRegions,
                             std::vector<ProteinRegion> &capRegions,
                             std::vector<ProteinRegion> &sevRegions);

void depolymerisationPointedForce(int i, std::vector<Actin> &actinVec,
                                  const bool tether, int &nActin,
                                  bool arpPool, ArpGrid &arpGrid,
                                  const bool steric,
                                  StericGrid &stericGrid,
                                  std::vector<ExcZone> &excZones,
                                  std::vector<MembraneWall> &memWalls,
                                  std::vector<Membrane> &membranes,
                                  Cortex &cortex,
                                  const bool bDynamics,
                                  std::vector<ProteinRegion> &branchRegions,
                                  std::vector<ProteinRegion> &nucRegions,
                                  std::vector<ProteinRegion> &capRegions,
                                  std::vector<ProteinRegion> &sevRegions);

void depolymerisationPointedGrid(std::vector<Actin> &actinVec, const double k_off,
                                 const double dt, const bool tether, int &nActin,
                                 const bool dissociate, GactinGrid &gActinGrid,
                                 bool arpPool, ArpGrid &arpGrid,
                                 const bool steric, StericGrid &stericGrid,
                                 std::vector<ExcZone> &excZones,
                                 std::vector<MembraneWall> &memWalls,
                                 std::vector<Membrane> &membranes,
                                 Cortex &cortex,
                                 const bool bDynamics,
                                 std::vector<ProteinRegion> &branchRegions,
                                 std::vector<ProteinRegion> &nucRegions,
                                 std::vector<ProteinRegion> &capRegions,
                                 std::vector<ProteinRegion> &sevRegions);

void autoDebranchBarb(Actin &parent, int &nActin, std::vector<Actin> &actinVec,
                      bool arpPool, ArpGrid &arpGrid,
                      const bool steric, StericGrid &stericGrid);

void autoDebranchPoint(Actin &parent, int &nActin, std::vector<Actin> &actinVec,
                       bool arpPool, ArpGrid &arpGrid,
                       const bool steric, StericGrid &stericGrid);

void dissociationRout(std::vector<int> dissociateIDs,
                      std::vector<Actin> &actinVec, int &nActin,
                      const bool steric, StericGrid &stericGrid);

#endif
