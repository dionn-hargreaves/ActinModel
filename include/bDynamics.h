/*
 *  bDynamics.h
 *
 *
 *  James Bradford
 *  University of Sheffield
 *  September 2017
 */

#ifndef BDYNAMICS_H
#define BDYNAMICS_H

#include <vector>
#include "Actin.h"

void BrownianDynamics(std::vector<Actin> &actinVec,
                      const std::vector<ExcZone> &excZones,
                      const std::vector<MembraneWall> &memWalls,
                      const std::vector<Membrane> &membranes,
                      Cortex &cortex, const bool steric,
                      StericGrid &stericGrid,
                      const double dt, double currTime, double temp,
                      double viscosity, bool tethering, bool crossLinking);

void brownianMotionRandWalk(Actin &actin, std::vector<Actin> &actinVec, const double viscosity,
                            const double temperature, const double dt,
                            const std::vector<ExcZone> &excZones,
                            const std::vector<MembraneWall> &memWalls,
                            const std::vector<Membrane> &membranes,
                            Cortex &cortex,
                            const bool steric, StericGrid &stericGrid, bool tethering,
                            bool crossLinking);

void initConfiguration(std::vector<Actin> &actinVec, const std::vector<ExcZone> &excZones,
                       const std::vector<MembraneWall> &memWalls,
                       const std::vector<Membrane> &membranes, Cortex &cortex,
                       const bool steric, StericGrid &stericGrid);


#endif
