/*
 *  Actin.h
 *
 *  Header file containing the declaration of the Actin class including
 *  the full definitions of the trivial member functions ("getters").
 *  The full definition of the class is found in the associated C++ file
 *  "Actin.cpp"
 *
 *  James Bradford
 *  University of Sheffield
 *  Jan 2017
 */
#ifndef ACTIN_H
#define ACTIN_H

#include "configHeader.h"
#include "GactinGrid.h"
#include "ProteinRegion.h"
#include "ExcZone.h"
#include "RNG.h"
#include "MembraneWall.h"
#include "ArpGrid.h"
#include "StericGrid.h"

class Actin
{
private:
    int m_id;
    double m_radius;
    bool m_steric;
    double m_stericRadius;
    int m_daughter_num;
    bool m_capped_b; // Barbed end capping
    bool m_capped_p; // Pointed end capping
    double m_length;
    int m_num_monomers;
    int m_num_subunits;
    double m_des_subunit_len;
    int m_branchdir;
    double m_birthtime;
    int m_parent_id;
    int m_structure_id;
    int m_subStructure_id;
    int m_crossStructure_id; // filaments on same side of crosslink
    double m_branchSubLen;
    int m_motherMonomer;
    int m_regionID;
    bool m_changed; // true if structure has changed since last dt
    bool m_flexible; // true if long enough to bend


    std::vector< std::array<double,3> > m_points;
    std::vector< std::array<double,3> > m_prevPoints; // points from a previous timestep - used for resetting
    // Vector of daughter IDs for each subunit therefore its a vector of vectors
    std::vector< std::vector<int> > m_daughterIDs;
    int m_branchSubUnit;
    std::vector<double> m_actualSubLengths;
    std::vector<double> m_prescribedSubLengths; // used for bending calculations

    // Bead-Rod model extras
    std::vector< std::array<double,2> > m_subunit_unitVecs; // unit vectors of the subunits

    std::array<double,3> m_Darray;
    std::array<double,3> m_EllipCentrePoint;
    Eigen::Matrix3d m_rotationMatrix;

    // Tethering
    std::vector< std::array<double,4> > m_tetherPoints;
    std::vector<int> m_tetherSubunits;
    //std::vector<double> m_tetherPos;
    //double m_distToTetBarb;
    int m_indexLeadingBarbTet;
    //double m_distToTetPoint;
    int m_indexLeadingPointTet;
    // Pre-determined distances, these are chosen randomly according to distribution
    // estimated from Tirf experiment
    int m_preDetDistToTetBarb;
    int m_preDetDistToTetPoint;

    int m_distToTetBarb;
    int m_distToTetPoint;
    // Discrete tether
    std::vector<int> m_tetherMono;

    // For steric grid
    // vector for each subunit, vector for each cell that sub is in
    std::vector< std::vector<int> > m_stericCells;

    // For crosslinking
    //std::vector<bool> m_cLinkFree; // vector of bools for each monomer false = occupied, true = free
    //int m_distToLeadCL; // distance of barbed end to leading crosslink/pointed end in number of monomers
    //std::vector<int> m_availMonoCLink; // vector of idxs of monomers available to link

    std::vector< std::array<int,6> > m_cLinkActinAndSites; // vector of arrays
    // which contain 0: monomer id, 1: subunit id, 2:other filament id, 3: other filaments crosslink idx, 4: other monomerid, 5: other filament subid

    std::vector< std::array<double,4> > m_cLinkPoints; // vector of arrays containing points in 2d space which are the location of crosslink ends
    // 0, 1 are crosslink location on the filament, 2,3 are crosslink location on the other filament
    bool m_cLinkMaster; // bool that is true for just one filament in the crosslink structure - used to avoid moving crosslink structures more than once per dt

    std::vector<int> m_availMonos; // vector of idxs of monomers available to branch or link
    std::vector<bool> m_monosCompat; // vector of bools for each monomer false = occupied, true = free
    int m_distToLeadBR;
    int m_distToLeadCL;

    std::vector<double> m_monoBirthTime; // vector of doubles for each monomer, the time the monomer was "created"



public:
    // Default
    Actin();

    // Regional filaments
    Actin(int ID, std::uniform_real_distribution<double> dist1,
          std::uniform_real_distribution<double> dist2, int regionID,
          const std::vector<ProteinRegion> &nucRegions,
          double birthtime = 0.0);

    Actin(int ID, std::uniform_real_distribution<double> xpos_dist,
          std::uniform_real_distribution<double> ypos_dist, int regionID,
          const std::vector<ProteinRegion> &nucRegions,
          std::normal_distribution<double> angle_dist, double birthtime = 0.0);

    // Branched filaments
    Actin(int ID, double birthtime, Actin &parent, int motherMonomer);

    // severing
    Actin(int ID, const std::array<double,2> &severPoint, int severMonomer,
          int severSub, double lenAlongSub, double severTime, const Actin &original);

    // Destructor (default)
     ~Actin();

    // Static constants belong to the class not the object, these are shared
    // and accessible publicly
    static const double s_stericRadius;
    static const double s_trueRadius;
    static const double s_monomerLength;
    static const int s_seedSize;
    static const int s_branchSeedSize;
    static int s_branchSpacing;
    static const double s_persistenceLength;
    static const double s_branchAngle;
    static double s_segmentationLength;
    static int s_structure_idGenerator;
    static int s_subStructure_idGenerator;
    static int s_crossStructure_idGenerator;
    static int s_total_subunits;
    static double s_tetherStiff;
    static double s_branchStiff;
    static double s_torsionCoeff; // angular stiffness
    static double s_cLStiff;
    static int s_cLSpacing;
    static double s_cLDist; // crosslink distance threshold and restLength

    static int s_maxSpacing; // this will be equal to the maximum of s_branchSpacing and s_cLSpacing


    // Forward declaration of non-trivial member functions
    int branchDir();
    void polymerise_barbed(std::vector<Actin> &actinVec, const bool steric,
                           StericGrid &stericGrid, const bool bDynamics,
                           const bool tether, const double currTime,
                           std::normal_distribution<double> tetherDistri,
                           std::vector<MembraneWall> &memWalls,
                           bool moveWall, double d_i,
                           std::vector<ProteinRegion> &branchRegions,
                           std::vector<ProteinRegion> &nucRegions,
                           std::vector<ProteinRegion> &capRegions,
                           std::vector<ProteinRegion> &sevRegions);

    void polymerise_pointed(std::vector<Actin> &actinVec, const bool steric,
                            StericGrid &stericGrid, const bool bDynamics,
                            const bool tether, const double currTime,
                            std::normal_distribution<double> tetherDistri,
                            std::vector<MembraneWall> &memWalls,
                            bool moveWall, double d_i,
                            std::vector<ProteinRegion> &branchRegions,
                            std::vector<ProteinRegion> &nucRegions,
                            std::vector<ProteinRegion> &capRegions,
                            std::vector<ProteinRegion> &sevRegions);

    void depolymerise_barbed(std::vector<Actin> &actinVec, int &nActin,
                             bool arpPool, ArpGrid &arpGrid, const bool steric,
                             StericGrid &stericGrid,
                             const std::vector<ExcZone> &excZones,
                             std::vector<MembraneWall> &memWalls,
                             std::vector<Membrane> &membranes,
                             Cortex &cortex,
                             const bool bDynamics,
                             const bool tether,
                             bool moveWall, double d_i,
                             std::vector<ProteinRegion> &branchRegions,
                             std::vector<ProteinRegion> &nucRegions,
                             std::vector<ProteinRegion> &capRegions,
                             std::vector<ProteinRegion> &sevRegions);

    void cancelBarbDepol(std::vector<Actin> &actinVec, const bool steric,
                         StericGrid &stericGrid);

    void cancelPointDepol(std::vector<Actin> &actinVec, const bool steric,
                          StericGrid &stericGrid);

    void depolymerise_pointed(std::vector<Actin> &actinVec, int &nActin,
                              bool arpPool, ArpGrid &arpGrid, const bool steric,
                              StericGrid &stericGrid,
                              const std::vector<ExcZone> &excZones,
                              std::vector<MembraneWall> &memWalls,
                              std::vector<Membrane> &membranes,
                              Cortex &cortex,
                              const bool bDynamics,
                              const bool tether,
                              bool moveWall, double d_i,
                              std::vector<ProteinRegion> &branchRegions,
                              std::vector<ProteinRegion> &nucRegions,
                              std::vector<ProteinRegion> &capRegions,
                              std::vector<ProteinRegion> &sevRegions);

    void addPointsBarbed(std::vector<Actin> &actinVec, const bool steric,
                         StericGrid &stericGrid);
    void removePointsBarbed(std::vector<Actin> &actinVec, const bool steric,
                            StericGrid &stericGrid, bool severing=false);
    void addPointsPointed(std::vector<Actin> &actinVec, const bool steric,
                          StericGrid &stericGrid);
    void removePointsPointed(std::vector<Actin> &actinVec, const bool steric,
                             StericGrid &stericGrid, bool severing=false);
    bool straightenFila(std::vector<Actin> &actinVec,
                        const std::vector<ExcZone> &excZones,
                        const std::vector<MembraneWall> &memWalls,
                        const std::vector<Membrane> &membranes,
                        Cortex &cortex,
                        const bool steric, StericGrid &stericGrid,
                        bool barbedEnd);

    bool straightenFila();
    void straightenEnd(bool barb);


    void updateSubLengths();
    void updateSubLength(int p);
    void updatePointsUp(int startpoint);
    void updatePointsDown(int startpoint);
    double getSumAngle(int subUnit, int centreSub) const;
    void regenerate_pos(const ProteinRegion &region);

    bool check_barbed_polymerisation_Grid(const std::vector<Actin> &actinVec,
                                          const bool steric,
                                          StericGrid &stericGrid);

    bool check_pointed_polymerisation_excVol_Grid(const std::vector<Actin> &actinVec,
                                                  const bool steric,
                                                  StericGrid &stericGrid);

    bool check_nucleation_excVol_Grid(const std::vector<Actin> &actinVec,
                                      const std::vector<ExcZone> &excZones,
                                      const std::vector<Membrane> &membranes,
                                      const std::vector<MembraneWall> &memWalls,
                                      Cortex &cortex,
                                      const bool steric, StericGrid &stericGrid);

    bool checkDiffusionStructure(std::vector<Actin> &actinVec,
                                 const std::vector<ExcZone> &excZones,
                                 const std::vector<MembraneWall> &memWalls,
                                 const std::vector<Membrane> &membranes,
                                 Cortex &cortex, const bool steric,
                                 StericGrid &stericGrid,
                                 const bool diffBool);

    bool checkDiffusionClStructure(std::vector<Actin> &actinVec,
                                 const std::vector<ExcZone> &excZones,
                                 const std::vector<MembraneWall> &memWalls,
                                 const std::vector<Membrane> &membranes,
                                 Cortex &cortex, const bool steric,
                                 StericGrid &stericGrid,
                                 const bool diffBool);

    bool check_diffusion_excVolGrid(const std::vector<Actin> &actinVec,
                                    const std::vector<ExcZone> &excZones,
                                    const std::vector<MembraneWall> &memWalls,
                                    const std::vector<Membrane> &membranes,
                                    Cortex &cortex, const bool steric,
                                    StericGrid &stericGrid, const bool diffBool);

    bool check_diffusion_excVolGrid2(const std::vector<Actin> &actinVec,
                                    const std::vector<ExcZone> &excZones,
                                    const std::vector<MembraneWall> &memWalls,
                                    const std::vector<Membrane> &membranes,
                                    Cortex &cortex, const bool steric,
                                    StericGrid &stericGrid, const bool diffBool,
                                    std::vector<int> &checkedFilas);

    bool check_bend(std::vector<Actin> &actinVec,
                    const std::vector<ExcZone> &excZones,
                    const std::vector<MembraneWall> &memWalls,
                    const std::vector<Membrane> &membranes,
                    Cortex &cortex,
                    const bool steric, StericGrid &stericGrid);

    bool check_bend(std::vector<Actin> &actinVec,
                    const std::vector<ExcZone> &excZones,
                    const std::vector<MembraneWall> &memWalls,
                    const std::vector<Membrane> &membranes,
                    Cortex &cortex,
                    const bool steric, StericGrid &stericGrid,
                    std::vector<int> &checkedFilas,
                    std::vector<std::array<int,2> > &checkedFilaSubs);


    bool check_bending_excVolGrid(const std::vector<Actin> &actinVec,
                                  const std::vector<ExcZone> &excZones,
                                  const std::vector<MembraneWall> &memWalls,
                                  const std::vector<Membrane> &membranes,
                                  Cortex &cortex,
                                  const bool steric, StericGrid &stericGrid,
                                  std::vector<std::array<int,2> > &checkedFilaSubs);


    bool checkBranchExcVolGrid(std::vector<Actin> &actinVec,
                               const std::vector<ExcZone> &excZones,
                               const std::vector<Membrane> &membranes,
                               const std::vector<MembraneWall> &memWalls,
                               Cortex &cortex,
                               const bool steric, StericGrid &stericGrid);

    bool check_min_dist_Grid(const std::vector<Actin> &actinVec,
                             std::vector<int> &cellsToCheck,
                             StericGrid &stericGrid);

    bool check_min_dist_Diff_Grid(int subid, const std::vector<Actin> &actinVec,
                                  std::vector<int> &cellsToCheck,
                                  StericGrid &stericGrid, const bool diffBool);

    bool check_min_dist_Diff_Grid2(int subid, const std::vector<Actin> &actinVec,
                                  std::vector<int> &cellsToCheck,
                                  StericGrid &stericGrid, const bool diffBool,
                                  std::vector<int> &checkedFilas);

    bool check_min_dist_BD_Grid(int subid, const std::vector<Actin> &actinVec,
                                std::vector<int> &cellsToCheck,
                                StericGrid &stericGrid,
                                std::vector<std::array<int,2> > &checkedFilaSubs);

    bool check_min_dist_branch_Grid(const std::vector<Actin> &actinVec,
                                   std::vector<int> &cellsToCheck,
                                   StericGrid &stericGrid);
    double FindeffL_box(const ProteinRegion &region);

    bool checkBranchSpacingPar(const Actin &parent);
    bool checkBranchSpacingSis(const Actin &sister, std::vector<Actin> &actinVec);
    bool checkBranchInRegion(const ProteinRegion &region);
    void checkAndMoveBranchDown(std::vector<Actin> &actinVec, int i);
    void checkAndMoveBranchUp(std::vector<Actin> &actinVec, int i);
    void checkAndMoveBranch(std::vector<Actin> &actinVec, int i);
    void moveBranchSub_poly_barb(std::vector<Actin> &actinVec);
    void moveBranchSub_poly_point(std::vector<Actin> &actinVec);
    void moveBranchSub_depoly_barb(std::vector<Actin> &actinVec);
    void moveBranchSub_depoly_point(std::vector<Actin> &actinVec);
    void moveBranchSub_down(std::vector<Actin> &actinVec);
    void moveBranchSub_up(std::vector<Actin> &actinVec);
    void moveBranchSub_Add_barb(std::vector<Actin> &actinVec);
    void moveBranchSub_Add_point(std::vector<Actin> &actinVec);
    void moveBranchSub_Rem_barb(std::vector<Actin> &actinVec);
    void moveBranchSub_Rem_point(std::vector<Actin> &actinVec);

    void moveBranchSubUnit_2Subs(std::vector<Actin> &actinVec);

    void regen_branch_pos(Actin &parent);
    void regen_branch_posRegions(Actin &parent, double branchSubPos,
                                 double branchPos, int branchSubUnit,
                                 int oldBranchSubUnit);

    void diffuse2D(double dxprime, double dyprime, std::vector<Actin> &actinVec,
                   Eigen::Matrix3d &rotationMatrix,
                   const std::vector<ExcZone> &excZones,
                   const std::vector<MembraneWall> & memWalls,
                   const std::vector<Membrane> &membranes,
                   Cortex &cortex,
                   const bool steric, StericGrid &stericGrid,
                   const bool translation,
                   bool crossLinking);

    void movebranches(double dx, double dy, std::vector<Actin> &actinVec,
                      const bool translation);
    void moveactin(double dx, double dy);
    void detach(std::vector<Actin> &actinVec, int nActin);
    void purgeBranchFromDaughterVec(int branchid);
    void changeStructure(std::vector<Actin> &actinVec, int nActin);
    void changeStructure(std::vector<Actin> &actinVec, int nActin,
                                std::vector<int> &changedFilas);

    void changeStructureALL(std::vector<Actin> &actinVec, int newStructureID);
    void changeSubStructure(std::vector<Actin> &actinVec);
    void changeCLStructure(std::vector<Actin> &actinVec);


    std::array<double,3> fitEllipsoids(const std::vector<Actin> &actinVec,
                                       const double viscosity, const double KT,
                                       Eigen::Matrix3d &rotationMatrix,
                                       std::array<double,3> &centrePoint,
                                       const bool translation,
                                       const bool rotation);

    std::array<double,3> fitEllipsoidsCLink(const std::vector<Actin> &actinVec,
                                       const double viscosity, const double KT,
                                       Eigen::Matrix3d &rotationMatrix,
                                       std::array<double,3> &centrePoint,
                                       const bool translation,
                                       const bool rotation);

    std::array<double,2> transEllipMobility(const double viscosity, double a,
                                            double b, double c, double P,
                                            double Q, double S);

    double rotEllipMobility(const double viscosity, double a,
                                           double b, double c, double P,
                                           double Q);

    std::array<double,3> kirkwoodDiff(const double viscosity,
                                      const double KT, const bool translation,
                                      const bool rotation);

    void rotateEllipsoid2D(double d_xyangle, std::vector<Actin> &actinVec,
                           Eigen::Matrix3d &rotationMatrix,
                           std::array<double,3> &centrePoint,
                           const std::vector<ExcZone> &excZones,
                           const std::vector<MembraneWall> &memWalls,
                           const std::vector<Membrane> &membranes,
                           Cortex &cortex,
                           const bool steric, StericGrid &stericGrid,
                           const bool rotation);

    void rotateEllipsoid2DTether(double d_xyangle, std::vector<Actin> &actinVec,
                                 Eigen::Matrix3d &rotationMatrix,
                                 const std::vector<ExcZone> &excZones,
                                 const std::vector<MembraneWall> &memWalls,
                                 const std::vector<Membrane> &membranes,
                                 Cortex &cortex,
                                 const bool steric, StericGrid &stericGrid,
                                 const bool rotation);

    void rotateEllipsoidALL(double d_xyangle, std::vector<Actin> &actinVec,
                                Eigen::Matrix3d &rotationMatrix,
                                std::array<double,3> &centrePoint,
                                const std::vector<ExcZone> &excZones,
                                const std::vector<MembraneWall> &memWalls,
                                const std::vector<Membrane> &membranes,
                                Cortex &cortex,
                                const bool steric, StericGrid &stericGrid,
                                const bool rotation);

    void rotateMiniEllipsoid(double d_xyangle, std::vector<Actin> &actinVec,
                                Eigen::Matrix3d &rotationMatrix,
                                std::array<double,3> &centrePoint,
                                const std::vector<ExcZone> &excZones,
                                const std::vector<MembraneWall> &memWalls,
                                const std::vector<Membrane> &membranes,
                                Cortex &cortex,
                                const bool steric, StericGrid &stericGrid,
                                const bool rotation);

    void rotateAllStructure(std::vector<Actin> &actinVec, double d_xyangle,
                            std::array<double,3> &centrePoint);

    void rotateCLStructure(std::vector<Actin> &actinVec, double d_xyangle,
                            std::array<double,3> &centrePoint);

    void updateBranchpositions(std::vector<Actin> &actinVec);
    void updateBranchSubpos(std::vector<Actin> &actinVec);
    void updateBranchSubpos_woRecursion(std::vector<Actin> &actinVec);

    void rotbranches(std::vector<Actin> &actinVec, const bool rotation);
    void rotbranchesInclFlex(std::vector<Actin> &actinVec, const bool rotation);
    void rotbranchesSub(std::vector<Actin> &actinVec, int p, const bool rotation);

    std::array<double,3> fActinDiffusion(const std::vector<Actin> &actinVec,
                                         const double viscosity, const double KT,
                                         const bool branching,
                                         Eigen::Matrix3d &rotationMatrix);
    double singlefilamentmobility(const double viscosity);
    double x_filamentmobility(const double viscosity);
    double y_filamentmobility(const double viscosity);
    std::array<double,3> ellipsoidmobility(const double viscosity, double a, double b, double c);
    std::array<double,3> fActinRotDiff(const std::vector<Actin> &actinVec,
                                       const double viscosity, const double KT);

    double rodDRot(const double viscosity);
    double rodDRotTether(const double viscosity);
    void rotate2D(double d_xyangle, std::vector<Actin> &actinVec,
                  const std::vector<ExcZone> &excZones,
                  const std::vector<MembraneWall> & memWalls,
                  const std::vector<Membrane> &membranes,
                  Cortex &cortex,
                  const bool steric, StericGrid &stericGrid,
                  const bool rotation);

    void rotate2DTether(double d_xyangle, std::vector<Actin> &actinVec,
                        const std::vector<ExcZone> &excZones,
                        const std::vector<MembraneWall> &memWalls,
                        const std::vector<Membrane> &membranes,
                        Cortex &cortex,
                        const bool steric, StericGrid &stericGrid,
                        const bool rotation);

    void rotate2DCLink(double d_xyangle, std::vector<Actin> &actinVec,
                       const std::vector<ExcZone> &excZones,
                       const std::vector<MembraneWall> &memWalls,
                       const std::vector<Membrane> &membranes,
                       Cortex &cortex,
                       const bool steric, StericGrid &stericGrid,
                       const bool rotation);

    std::array<double,2> findCentrePoint(int &centreSubUnit);
    std::array<double,2> findCentrePoint(int &centreSubUnit) const;

    void sever(std::vector<Actin> &actinVec, double currTime, int severMonomer,
               int &nActin, const bool steric, StericGrid &stericGrid,
               GactinGrid &gActinGrid,
               const bool arpPool, ArpGrid &arpGrid);

    void changeAngles(int index, double angle);
    void changeToNewAngle(int index, int centreSub, double newAngle);
    void initialiseLength(int numMonomers,std::vector<Actin> &actinVec,
                          const bool steric, StericGrid &stericGrid,
                          const std::vector<Membrane> &membranes,
                          Cortex &cortex,
                          const std::vector<ExcZone> &excZones,
                          const std::vector<MembraneWall> &memWalls,
                          const bool tethering);

    void addPoint(std::array<double,3> point) { m_points.push_back(point); }
    void addActualSubLength(double subLength) { m_actualSubLengths.push_back(subLength); }
    void addPresSubLength(double subLength) { m_prescribedSubLengths.push_back(subLength); }
    void addEmptyVecToDaughterID() { m_daughterIDs.push_back(std::vector<int>()); };
    void addBranchtoParentVector(int subunit, int ID) { m_daughterIDs[subunit].push_back(ID); }
    void popBranchoffParentVector(int subunit) { m_daughterIDs[subunit].pop_back(); }

    // For steric grid
    void addCellToActinSub(int subunit, int cell) { m_stericCells[subunit].push_back(cell); }
    std::vector<int> getStericCells(int subunit) const { return m_stericCells[subunit]; }
    void clearStericCellsSub(int subunit) { m_stericCells[subunit].clear(); }

    void updateUnitVecs();
    void updateUnitVecsNoCheck();
    void updateUnitVec(int p);
    Eigen::VectorXd calcBendingForces(double temp, int N);

    Eigen::MatrixXd buildDMatrixAniso(double temp, double viscosity, int N);
    Eigen::VectorXd getRand_dR(double dt, double temp, double viscosity, int N);

    Eigen::MatrixXd getBMatrix(int N);
    void addTetherPointNuc();
    void addTetherPointBarb();
    void addTetherPointPointed();
    void remTetherB(std::vector<Actin> &actinVec);
    void remTetherP(std::vector<Actin> &actinVec);
    void calcTetherPos();
    void updateTetherLoc();
    void checkAndMoveTetherDown(int subOne=-1, int subTwo=-2, bool checkAll=true);
    void checkAndMoveTetherUp(int subOne=-1, int subTwo=-2, bool checkAll=true);
    void moveTetherSub_down();
    void moveTetherSub_up();
    void moveTetherPoint_poly_barb();
    void moveTetherPoint_poly_point();
    void moveTetherPoint_depoly_barb();
    void moveTetherPoint_depoly_point();
    void moveTetherPoint_Add_barb();
    void moveTetherPoint_Add_point();
    void moveTetherPoint_Rem_barb();
    void moveTetherPoint_Rem_point();
    Eigen::VectorXd calcTetherForces(int N, bool tethering);
    double calcDistancetoTether(int subunitID, double projPos, bool down);

    Eigen::VectorXd calcBranchForces(std::vector<Actin> &actinVec, int N);
    Eigen::VectorXd calcMotherForces(std::vector<Actin> &actinVec, int N);

    // Definitions of trivial member functions
    double getLength() const { return m_length; }
    double getStericRadius() const { return m_stericRadius; }
    double getRadius() const { return m_radius; }
    int getNumMonomers() const { return m_num_monomers; }
    int getID() const { return m_id; }

    std::vector< std::array<double,3> > getPoints() const { return m_points; }
    std::array<double,3> getPoint(int pointID) const { return m_points[pointID]; }
    void setPointedEndPoints(double x, double y) { m_points[0] = { x, y, 0 }; }
    void setBarbedEndPoints(double x, double y) { m_points[m_num_subunits] = { x, y, 0 }; }
    std::array<double,3> getBarbedEnd() const { return m_points[m_num_subunits]; }
    std::array<double,3> getPreBarbedEnd() const { return m_points[m_num_subunits-1]; }
    std::array<double,3> getPointedEnd() const { return m_points[0]; }
    std::array<double,3> getPrePointedEnd() const { return m_points[1]; }


    std::vector< std::array<double,2> > getUnitVecs() const { return m_subunit_unitVecs; }
    std::array<double,2> getUnitVec(int subUnit) const { return m_subunit_unitVecs[subUnit]; }
    void setUnitVec(int subunit, std::array<double,2> unitVec) { m_subunit_unitVecs[subunit] = unitVec; }



    double getbirthTime() const { return m_birthtime; }
    int getParentID() const { return m_parent_id; }
    void setParentID(int parentID) { m_parent_id = parentID; }
    int getBranchDir() const { return m_branchdir; }
    int getStructureID() const { return m_structure_id; }
    void setStructureID(int structure_id) { m_structure_id = structure_id; }
    int getSubStructureID() const { return m_subStructure_id; }
    void setSubStructureID(int subStructure_id) { m_subStructure_id = subStructure_id; }
    int getCLStructureID() const { return m_crossStructure_id; }
    void setCLStructureID(int structure_id) { m_crossStructure_id = structure_id; }
    int getDaughternum() const { return m_daughter_num; }
    void incrementDaughternum() { m_daughter_num++; }
    void decrementDaughternum() { m_daughter_num--; }
    bool getBarbedCapped() const { return m_capped_b; }
    bool getPointedCapped() const { return m_capped_p; }
    void setBarbedCapped() { m_capped_b = true; }
    void setPointedCapped() { m_capped_p = true; }
    void setBarbedUncapped() { m_capped_b = false; }
    void setPointedUncapped() { m_capped_p = false; }

    //void setTetherPos(int i, double tetherpos) { m_tetherPos[i] = tetherpos; }
    //std::vector<double> getTetherPos() const { return m_tetherPos; }
    std::vector<int> getTetherMono() const { return m_tetherMono; }
    unsigned int getNumTethers() const { return m_tetherMono.size(); }
    void setTetherSub(int i, int subid) { m_tetherSubunits[i] = subid; }
    int getTetherSub(int i) const { return m_tetherSubunits[i]; }
    std::array<double,4> getTetherPoints(int i) const { return m_tetherPoints[i]; }
    int getTetherDistB() const { return m_distToTetBarb; }
    int getTetherDistP() const { return m_distToTetPoint; }
    void addToTetherDistB(int dl) { m_distToTetBarb += dl; }
    void addToTetherDistP(int dl) { m_distToTetPoint += dl; }
    void setTetherDistB(int l) { m_distToTetBarb = l; }
    void setTetherDistP(int l) { m_distToTetBarb = l; }
    void setPreDetTethDistB(int d) { m_preDetDistToTetBarb = d; }
    void setPreDetTethDistP(int d) { m_preDetDistToTetPoint = d; }
    int getPreDetTethDistB() const { return m_preDetDistToTetBarb; }
    int getPreDetTethDistP() const { return m_preDetDistToTetPoint; }
    void shiftTetherMonosUp();
    void shiftTetherMonosDown();


    int getNumSubs() const { return m_num_subunits; }
    void setNumSubs(int numSubs) { m_num_subunits = numSubs; }
    std::vector<double> getActualSubLengths() const { return m_actualSubLengths; }
    void setActualSubLength(int subid, double len) { m_actualSubLengths[subid] = len; }
    std::vector<double> getPresSubLengths() const { return m_prescribedSubLengths; }
    void setPresSubLength(int subid, double len) { m_prescribedSubLengths[subid] = len; }
    std::vector<int> getBranchIDSubVector(int subunit) const { return m_daughterIDs[subunit]; }
    int getRegionID() const { return m_regionID; }
    void setChangedBoolTrue() { m_changed = true; }
    void setChangedBoolFalse() { m_changed = false; }
    bool getChanged() const { return m_changed; }
    void setDarray(const std::array<double,3> &D_array) { m_Darray = D_array; }
    std::array<double,3> getDarray() const { return m_Darray; }
    void setRotMat(const Eigen::Matrix3d &rotMatrix) { m_rotationMatrix = rotMatrix; }
    Eigen::Matrix3d getRotMat() const { return m_rotationMatrix; }
    void setEllipCentre(const std::array<double,3> &centrePoint) { m_EllipCentrePoint = centrePoint; }

    std::array<double,3> getEllipCentre() const { return m_EllipCentrePoint; }
    void turnStericOff() { m_steric = false; }
    bool getStericBool() const { return m_steric; }

    void setPoints(int index, std::array<double,2> point)
    {
        m_points[index][0] = point[0];
        m_points[index][1] = point[1];
    }

    int findSubunit(int monoID, double &lenAlongSub);
    std::array<double,2> findMonoCoord(int monoID);
    std::array<double,4> findMonoStartEndCoord(int monoID);
    void shiftAvailMonosUp();
    void addAvailMono();
    void shiftAvailMonosDown();
    void shiftMonosCompatUp();
    void shiftMonosCompatDown();
    void shiftMonosBirthTimeUp();
    void shiftMonosBirthTimeDown();

    int getDistToLeadBR() const { return m_distToLeadBR; }
    bool getMonoCompat(int monoID) const { return m_monosCompat[monoID]; }
    std::vector<int> getAvailMonoVec() const { return m_availMonos; };
    int getAvailMonoID(int idx) const { return m_availMonos[idx]; }
    void adjustBranchCompat(int monoID);
    void shiftMotherMonosUp(std::vector<Actin> &actinVec);
    void shiftMotherMonosDown(std::vector<Actin> &actinVec);
    void incrementMotMono() { ++m_motherMonomer; }
    void decrementMotMono() { --m_motherMonomer; }
    int getMotherMonoID() const { return m_motherMonomer; }
    void setMotherMonoID(int newID) { m_motherMonomer = newID; }
    int getBranchSubunit() const { return m_branchSubUnit; }
    void setBranchSubunit(int BrSub) { m_branchSubUnit = BrSub; }
    void incrementBRsub() { ++m_branchSubUnit; }
    void decrementBRsub() { --m_branchSubUnit; }
    void recalcDistToBR(std::vector<Actin> &actinVec, const int nActin);
    void recalcDistToBR2(std::vector<Actin> &actinVec);
    std::vector<int> getUnavailMonos(std::vector<Actin> &actinVec);
    void recalcMonoCompat(std::vector<Actin> &actinVec, int monoID);
    void recalcMonoCompatALL(std::vector<Actin> &actinVec);
    void sortAvailMonos();
    int findMonoInReg(int monoAlongEffL, const ProteinRegion &region);

    double getMonoBirthTime(int monoID) const { return m_monoBirthTime[monoID]; }

    void reduceID(std::vector<Actin> &actinVec);
    void reduceParentID() { m_parent_id--; };
    void checkandReduceDaughtID(int doomedID);
    void savePointsToPrev() { m_prevPoints = m_points; }
    void resetPointsToPrev() { m_points = m_prevPoints; }
    void saveOldPos(std::vector<Actin> &actinVec);
    void saveOldPos(std::vector<Actin> &actinVec, std::vector<int> &calledFilas);
    void resetToOldPos(std::vector<Actin> &actinVec);
    void resetToOldPos(std::vector<Actin> &actinVec, std::vector<int> &calledFilas);

    void saveBranchPos(std::vector<Actin> &actinVec, const bool rotation);
    void saveBranchPosInclFlex(std::vector<Actin> &actinVec, const bool rotation);
    void saveBranchPosSub(std::vector<Actin> &actinVec, int p, const bool rotation);
    void resetBranchPos(std::vector<Actin> &actinVec, const bool rotation);
    void resetBranchPosInclFlex(std::vector<Actin> &actinVec, const bool rotation);
    void resetBranchPosSub(std::vector<Actin> &actinVec, int p, const bool rotation);

    void createCrossLink(Actin &otherFila, int monoID, int subID,
                         int otherMonoID, int otherSubID,
                         std::vector<Actin> &actinVec);

    void setMasterBoolFalse(std::vector<Actin> &actinVec);
    bool checkForMasterInStructure(std::vector<Actin> &actinVec);
    bool checkForCrossLinksInStructure(std::vector<Actin> &actinVec);
    bool checkForCrossLinksInStructure(std::vector<Actin> &actinVec,
                                       std::vector<int> &checkedFilas);

    bool findNewMaster(std::vector<Actin> &actinVec);
    bool findNewMaster(std::vector<Actin> &actinVec,
                       std::vector<int> &checkedFilas);

    bool checkForSameStructure(std::vector<Actin> &actinVec, int otherFilaID);

    bool checkForSameStructure(std::vector<Actin> &actinVec, int otherFilaID,
                               std::vector<int> &checkedFilas);

    int countCrossLinksInCLStruct(std::vector<Actin> &actinVec);
    std::array<double,3> getCLCoordInCLStruct(std::vector<Actin> &actinVec);

    void addToLinkVect(int mono, int subid, int otherID, int otherCLIDX,
                       int otherMono, int otherSubID);

    void setMonoToUnavailable(int mono);
    void remMonoFromAvail(int mono);
    Eigen::VectorXd calcCrossLinkForces(int N, bool crosslinking,
                                        std::vector<Actin> &actinVec);


    void updateCrossLinkLoc(std::vector<Actin> &actinVec);
    void updateCrossLinkLocs(std::vector<Actin> &actinVec);

    void setCrossLinkPoint(int cLID, double pX, double pY);
    void setOtherCrossLinkPoint(int cLID, double pX, double pY);
    void addEmptyToCLPoints() { m_cLinkPoints.push_back(std::array<double,4>()); }
    void recalcDistToLeadCL(int cLMono);
    void availMonosAddCLSp();

    // For crosslinking
    int getDistToLeadCL() const { return m_distToLeadCL; }
    void setDistToLeadCL(int newDist) { m_distToLeadCL = newDist; }
    //std::vector<int> getAvailCLSites() const { return m_availMonoCLink; }
    //int getAvailCLSiteID(int id) const { return m_availMonoCLink[id]; }
    int getNumCLinks() const { return m_cLinkActinAndSites.size(); }
    std::vector< std::array<int,6> > getCLinkActinAndSites() const { return m_cLinkActinAndSites; }
    std::array<int,6> getCLinkActinAndSite(int clID) const { return m_cLinkActinAndSites[clID]; }

    std::vector< std::array<double,4> > getCLinksPoints() const { return m_cLinkPoints; }
    std::array<double,4> getCLinksPoint(int i) const { return m_cLinkPoints[i]; }
    //bool getCLinkFree(int monoID) const { return m_cLinkFree[monoID]; }
    //std::vector<bool> getCLinkFreeVec() const { return m_cLinkFree; }
    void eraseCLink(int idx);



    bool getCLinkMaster() const { return m_cLinkMaster; }
    void turnOffCLinkMaster() { m_cLinkMaster = false; }
    void setClSub(int i, int subid, std::vector<Actin> &actinVec);
    void setClOther(int i, int filaID) { m_cLinkActinAndSites[i][2] = filaID; }
    void setClOtherIDX(int i, int idx) { m_cLinkActinAndSites[i][3] = idx; }
    void setClMonoOther(int i, int monoID) { m_cLinkActinAndSites[i][4] = monoID; }
    void setClSubOther(int i, int subID) { m_cLinkActinAndSites[i][5] = subID; }
    void setSelfMasterBoolTrue() { m_cLinkMaster = true; }
    void setSelfMasterBoolFalse() { m_cLinkMaster = false; }


    void shiftClMonosUp(std::vector<Actin> &actinVec);
    void shiftClMonoUp(int clIDX);
    void shiftClMonosDown(std::vector<Actin> &actinVec);
    void shiftClMonoDown(int clIDX);
    void checkAndMoveClDown(std::vector<Actin> &actinVec, int subOne=-1, int subTwo=-2, bool checkAll=true);
    void checkAndMoveClUp(std::vector<Actin> &actinVec, int subOne=-1, int subTwo=-2, bool checkAll=true);
    void checkAndMoveCl(std::vector<Actin> &actinVec);
    void moveClSub_down(std::vector<Actin> &actinVec);
    void moveClSub_up(std::vector<Actin> &actinVec);
    void moveClPoint_poly_barb(std::vector<Actin> &actinVec);
    void moveClPoint_poly_point(std::vector<Actin> &actinVec);
    void moveClPoint_depoly_barb(std::vector<Actin> &actinVec);
    void moveClPoint_depoly_point(std::vector<Actin> &actinVec);
    void moveClPoint_Add_barb(std::vector<Actin> &actinVec);
    void moveClPoint_Add_point(std::vector<Actin> &actinVec);
    void moveClPoint_Rem_barb(std::vector<Actin> &actinVec);
    void moveClPoint_Rem_point(std::vector<Actin> &actinVec);


    std::array<int,2> findMonoIDLimitsOfSub(int subid);
    void rotRigidCLinked(std::vector<Actin> &actinVec, const bool rotation);
    void rotRigidCLinked(std::vector<Actin> &actinVec, const bool rotation,
                         std::vector<int> &checkedFilas);

    bool getFlexibleBool() const { return m_flexible; }
    void setFlexibleBool(bool flex) { m_flexible = flex; }

};

#endif
