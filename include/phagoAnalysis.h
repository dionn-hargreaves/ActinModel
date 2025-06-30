





#ifndef PHAGOANALYSIS_H
#define PHAGOANALYSIS_H

double calcEngulfment(std::vector<Membrane> &membranes,
                      std::vector<ExcZone> &targets);

int calcNumContactTips(std::vector<Membrane> &membranes,
                       std::vector<ExcZone> &targets, std::vector<Actin> &actins);

int calcNumContactTipsInRegion(std::vector<Membrane> &membranes,
                               std::vector<Actin> &actins, double regXmin,
                               double regXmax, double regYmin, double regYmax);
#endif
