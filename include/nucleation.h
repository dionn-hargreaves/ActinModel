/*
 *  nucleation.h
 *
 *
 *  James Bradford
 *  University of Sheffield
 *  May 2017
 */

#ifndef NUCLEATION_H
#define NUCLEATION_H

#include "Actin.h"
#include "ProteinRegion.h"
#include "ExcZone.h"

void initNuc(std::vector<Actin> &actinvec, const int nActin,
             const std::vector<ProteinRegion> &nucRegions,
             const std::vector<ExcZone> &excZones,
             const std::vector<Membrane> &membranes,
             const std::vector<MembraneWall> &memWalls,
             Cortex &cortex,
             const bool steric, StericGrid &stericGrid,
             const bool tether,
             std::normal_distribution<double> tetherDistri,
             bool nucDir=0, double mean=0, double stdev=0);



void nuc(std::vector<Actin> &actinvec, int &nActin, const double k_NucDens,
         const double dt,
         const double gActinConc, const double Time,
         const int nActinLimit,
         const std::vector<ProteinRegion> &nucRegions,
         const std::vector<ExcZone> &excZones,
         const std::vector<Membrane> &membranes,
         const std::vector<MembraneWall> &memWalls,
         Cortex &cortex,
         const bool steric, StericGrid &stericGrid,
         const bool tether, std::normal_distribution<double> tetherDistri,
         bool nucDir=0, double mean=0, double stdev=0);

void nucGrid(std::vector<Actin> &actinvec, int &nActin,
             const double k_NucDens, const double dt, const double Time,
             const int nActinLimit,
             const std::vector<ProteinRegion> &nucRegions,
             const std::vector<ExcZone> &excZones,
             const std::vector<Membrane> &membranes,
             const std::vector<MembraneWall> &memWalls,
             Cortex &cortex,
             const bool steric, StericGrid &stericGrid,
             const bool tether, std::normal_distribution<double> tetherDistri,
             GactinGrid &gActinGrid, bool nucDir=0, double mean=0, double stdev=0);

void setSteric(Actin &filament, const double quasi3dProb);

#endif
