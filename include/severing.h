/*
 *  severing.h
 *
 *
 *  James Bradford
 *  University of Sheffield
 *  June 2017
 */

#ifndef SEVERING_H
#define SEVERING_H

#include "configHeader.h"

void severfilament(std::vector<Actin> &actinvec, const double k_sever,
                   const double t_step, int &nActin, const int nActinLimit,
                   bool &severing, const double currTime, const bool steric,
                   StericGrid &stericGrid,
                   GactinGrid &gActinGrid, const bool arpPool, ArpGrid &arpGrid);

void severfilamentRegions(std::vector<Actin> &actinvec, const double k_sever,
                          const double t_step,
                          int &nActin, const int nActinLimit, bool &severing,
                          const double currTime,
                          std::vector<ProteinRegion> &sevRegions);

void severfilamentRegions(std::vector<Actin> &actinvec, const double k_sever,
                          const double t_step, int &nActin, const int nActinLimit,
                          bool &severing, const double currTime, const bool steric,
                          StericGrid &stericGrid,
                          GactinGrid &gActinGrid, const bool arpPool, ArpGrid &arpGrid,
                          std::vector<ProteinRegion> &sevRegions);
#endif
