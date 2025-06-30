/*
 *  crosslinking.h
 *
 *
 *  James Bradford
 *  University of Sheffield
 *  Feb 2020
 */

#ifndef CROSSLINKING_H
#define CROSSLINKING_H

void crosslinkingRoutine(std::vector<Actin> &actinVec, const double k_cLink,
                      const double tStep, const double currTime,
                      const bool steric, StericGrid &stericGrid,
                      bool &crossLinking);

void unLinkBase(int i, int j, std::vector<Actin> &actinVec);
void unLink(std::vector<Actin> &actinVec, const double k_unLink,
            const double tStep, const double currTime,
            const bool steric, StericGrid &stericGrid,
            bool &crossLinking);

void autoUnlinkBarb(Actin &actin, std::vector<Actin> &actinVec);
void autoUnlinkPoint(Actin &actin, std::vector<Actin> &actinVec);

#endif
