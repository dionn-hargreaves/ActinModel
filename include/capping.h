/*
 *  capping.h
 *
 *
 *  James Bradford
 *  University of Sheffield
 *  May 2017
 */

#ifndef CAPPING_H
#define CAPPING_H

#include <vector>
#include "Actin.h"

void capBarbed(std::vector<Actin> &actinvec, const double pCap);
void capBarbedRegion(std::vector<Actin> &actinvec, double pCap,
                     const std::vector<ProteinRegion> &capRegions);

void capBarbedAntiRegion(std::vector<Actin> &actinvec, double pCap,
                     const std::vector<ProteinRegion> &antiCapRegions);
void capPointed(std::vector<Actin> &actinvec, const double pCap);
void capPointedRegion(std::vector<Actin> &actinvec, double pCap,
                      const std::vector<ProteinRegion> &capRegions);

void capPointedAntiRegion(std::vector<Actin> &actinvec, double pCap,
                     const std::vector<ProteinRegion> &antiCapRegions);
void uncapBarbed(std::vector<Actin> &actinvec, const double pUncap);
void uncapPointed(std::vector<Actin> &actinvec, const double pUncap);


#endif
