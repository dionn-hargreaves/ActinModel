

#ifndef MEMBRANE_H
#define MEMBRANE_H

#include "ExcZone.h"
#include "StericGrid.h"
#include "MembraneWall.h"
#include "ArpGrid.h"

class Membrane
{
private:
    bool m_exist;
    bool m_cellMem; // only true for cell membrane, false for phagosome
    int m_id;
    double m_ypos;
    double m_length;
    int m_numPoints;
    double m_segLength;
    double m_bendingMod;
    double m_memTension;
    double m_stiffness;
    double m_thickness;
    double m_zWidth;
    double m_effPersLen; // effective persistence length Lp = K_bend*2R / kbT

    bool m_activate;
    double m_activeDist;
    bool m_fusable; // can this membrane undergo fusion?
    bool m_initCortex;

    double m_tetherStiff; // for tethering of cortex to membrane
    std::vector<int> m_tetherSubIDs;
    std::vector<double> m_tetherRelDists;

    std::array<double,2> m_COM; // centre of mass

    std::vector< std::array<double,2> > m_points;
    std::vector< std::array<double,2> > m_subunit_unitVecs; // unit vectors of the subunits

    std::vector<double> m_actualSubLengths;
    std::vector<double> m_prescribedSubLengths;
    std::vector<std::vector<int>> m_coupBranchRegIDs;
    std::vector<std::vector<int>> m_coupNucRegIDs;
    std::vector<std::vector<int>> m_coupCapRegIDs;
    std::vector<std::vector<int>> m_coupAntiCapRegIDs;
    std::vector<std::vector<int>> m_coupSevRegIDs;

    // For steric grid
    // vector for each subunit, vector for each cell that sub is in
    std::vector< std::vector<int> > m_stericCells;


    Eigen::VectorXd calcBendingForces(double temp);
    Eigen::MatrixXd buildDMatrix(double temp, double viscosity);
    Eigen::VectorXd getRand_dR(double dt, Eigen::MatrixXd D);

    Eigen::MatrixXd getBMatrix();
    bool checkExVol(const std::vector<Actin> &actinvec);
    bool checkExVolSUB(const std::vector<Actin> &actinvec, int pointID);
    bool checkExVolSUBGrid(const std::vector<Actin> &actinvec, int pointID,
                           StericGrid &stericGrid);
    bool checkExVolSUBGridMembrane(const std::vector<Membrane> &membranes,
                                   int pointID, StericGrid &stericGrid);

    bool checkExVolSUBGridCortex(const Cortex &cortex,
                                 int pointID, StericGrid &stericGrid);

    bool check_min_dist_Grid(int subid, const std::vector<Actin> &actinvec,
                                                  std::vector<int> &cellsToCheck,
                                                  StericGrid &stericGrid);

    std::array<int,2> check_min_dist_Grid_Identify(int subid, const std::vector<Actin> &actinvec,
                                                   std::vector<int> &cellsToCheck,
                                                   StericGrid &stericGrid);

    Eigen::VectorXd calcBSForce(double temp, double dt, double eta);
    Eigen::VectorXd calcConsForce(double temp, double dt, double eta,
                                                    Eigen::VectorXd f_guess);


public:
    Membrane();

    Membrane(double ypos, double length, double temp, double bendMod,
             bool fusable = true, bool active = false, double activeDist = 1E-8,
             bool initCortex = false, double memTension = 1E-5);

    Membrane(int numPoints, double temp, std::array<double,2> fusePoint,
             Membrane &membrane, std::array<int,3> fusePointsToAv,
             int fuseType, int ID);

    ~Membrane(); // destructor

    bool getExist() const { return m_exist; }
    double getYpos() const { return m_ypos; }
    double getRadius() const { return (m_length/(2*M_PI)); }
    std::vector<double> initCircle();
    void calcSubunitVecs();
    void calcSubLengths();


    std::vector< std::array<double,2> > getPoints() const { return m_points; }
    double getNumPoints() const { return m_numPoints; };
    double getLength() const { return m_length; }
    double getThickness() const { return m_thickness; }
    double getSegLength() const { return m_segLength; }
    double get1DBendMod() const { return (m_bendingMod*m_zWidth); }
    double getEffLp() const { return m_effPersLen; }
    double getStiffness() const { return m_stiffness; }
    double getMemTension() const { return m_memTension; }
    std::vector<std::array<double,2>> getUnitVecs() const { return m_subunit_unitVecs; }
    std::vector<double> getActualSubLengths() const { return m_actualSubLengths; }
    std::vector<double> getPresSubLengths() const { return m_prescribedSubLengths; }
    bool getActive() const { return m_activate; }
    std::array<double,2> getCOM() const { return m_COM; }

    int getID() const { return m_id; }
    bool getFusable() const { return m_fusable; }
    bool getInitCortex() const { return m_initCortex; }

    void fluctuate(double temp, double viscosity, double dt,
                   std::vector<Actin> &actinvec,
                   std::vector<ExcZone> &excZones,
                   const bool steric, unsigned int &timesFused,
                   StericGrid &stericGrid, double currTime,
                   std::vector<ProteinRegion> &branchRegions,
                   bool &branching, double arpConc,
                   std::vector<ProteinRegion> &capRegions,
                   std::vector<ProteinRegion> &antiCapRegions,
                   std::vector<ProteinRegion> &nucRegions,
                   std::vector<ProteinRegion> &sevRegions,
                   std::vector<Membrane> &membranes,
                   int &nActin, const bool tether,
                   bool arpPool, ArpGrid &arpGrid,
                   std::vector<MembraneWall> &memWalls,
                   std::ofstream &fusionTime,
                   const bool bDynamics, Cortex &cortex, bool COMadjust,
                   const bool memSprings,
                   const bool activeBranch,
                   const double activeBranchHeight,
                   const bool activeNuc, const double activeNucHeight,
                   const bool activeCap, const double activeCapHeight,
                   const bool activeAntiCap,
                   const double activeAntiCapHeight,
                   const bool activeSever,
                   const double activeSeverHeight);

    void adjustCOM(std::vector<Actin> &actinvec,
                   std::vector<ExcZone> &excZones,
                   std::vector<Membrane> &membranes, Cortex &cortex,
                   const bool steric, StericGrid &stericGrid,
                   std::vector<ProteinRegion> &branchRegions,
                   std::vector<ProteinRegion> &capRegions,
                   std::vector<ProteinRegion> &antiCapRegions,
                   std::vector<ProteinRegion> &nucRegions,
                   std::vector<ProteinRegion> &sevRegions);


    static bool s_checkExVolBarbPoly(const Actin &filament,
                                     std::vector<Membrane> &membranes,
                                     StericGrid &stericGrid);

    bool checkExVolBarbPoly(const Actin &filament) const;
    bool checkExVolBarbPolyGrid(const Actin &filament,
                                StericGrid &stericGrid) const;

    static bool s_checkExVolPointPoly(const Actin &filament,
                                      std::vector<Membrane> &membranes,
                                      StericGrid &stericGrid);
    bool checkExVolPointPoly(const Actin &filament) const;
    bool checkExVolPointPolyGrid(const Actin &filament,
                                 StericGrid &stericGrid) const;

    bool checkExVolNuc(const Actin &actin) const;

    static bool s_checkExVolNuc(const Actin &actin,
                                const std::vector<Membrane> &membranes,
                                std::vector<int> cellsToCheck,
                                StericGrid &stericGrid);

    bool checkExVolNuc_Grid(const Actin &actin,
                            std::vector<int> cellsToCheck,
                            StericGrid &stericGrid) const;

    bool checkExVolBM(const Actin &actin) const;
    bool checkInOut(const Actin &actin) const;
    bool checkInOutAll(const std::vector<Actin> &actinVec) const;
    bool checkInOutPoly(gte::Segment<3,double> monomer) const;
    bool checkPointIn(std::array<double,2> point) const;

    std::vector<int> getStericCells(int subunit) const { return m_stericCells[subunit]; }
    void clearStericCellsSub(int subunit) { m_stericCells[subunit].clear(); }
    void addCellToMemSub(int subunit, int cell) { m_stericCells[subunit].push_back(cell); }

    static bool s_checkExVolBM(int subid, const Actin &actin,
                               const std::vector<Membrane> &membranes,
                               std::vector<int> cellsToCheck,
                               StericGrid &stericGrid);

    bool check_min_dist_Actin_Grid(int subid, const Actin &actin,
                                   std::vector<int> cellsToCheck,
                                    StericGrid &stericGrid) const;

    std::vector<int> getBranchRegID(int sub) const { return m_coupBranchRegIDs[sub]; }
    std::vector<int> getNucRegID(int sub) const { return m_coupNucRegIDs[sub]; }
    std::vector<int> getCapRegID(int sub) const { return m_coupCapRegIDs[sub]; }
    std::vector<int> getACapRegID(int sub) const { return m_coupAntiCapRegIDs[sub]; }

    void addToBranchReg(int sub, int id) { m_coupBranchRegIDs[sub].push_back(id); }
    void addToNucReg(int sub, int id) { m_coupNucRegIDs[sub].push_back(id); }
    void addToCapReg(int sub, int id) { m_coupCapRegIDs[sub].push_back(id); }
    void addToAntiCapReg(int sub, int id) { m_coupAntiCapRegIDs[sub].push_back(id); }
    void addToSevReg(int sub, int id) { m_coupSevRegIDs[sub].push_back(id); }


    bool calcDistToTarsForActive(const std::vector<ExcZone> &excZones,
                                 int subID);

    bool checkFusion(unsigned int &timesFused,
                     unsigned int numTars, int &fuseType,
                     std::array<double,2> &fusePoint,
                     std::array<int,3> &fusePointsToAv);

    void doFusion(double temp, std::vector<Membrane> &membranes,
                  StericGrid &stericGrid, unsigned int &timesFused,
                  unsigned int numTars,
                  std::vector<ProteinRegion> &branchRegions,
                  bool &branching,
                  std::vector<ProteinRegion> &capRegions,
                  std::vector<ProteinRegion> &nucRegions,
                  std::vector<ProteinRegion> &antiCapRegions,
                  std::vector<ProteinRegion> &sevRegions,
                  std::vector<Actin> &actinVec,
                  int &nActin, const bool steric,
                  bool arpPool, ArpGrid &arpGrid,
                  std::vector<ExcZone> &excZones,
                  std::vector<MembraneWall> &memWalls,
                  Cortex &cortex,
                  const bool tether, double currTime,
                  std::ofstream &fusionTime,
                  const bool bDynamics, int &fuseType,
                  std::array<double,2> &fusePoint,
                  std::array<int,3> &fusePointsToAv);

    bool getCellMemBool() const { return m_cellMem; }


    void initActRegion(std::vector<ProteinRegion> &branchRegions,
                       std::vector<ProteinRegion> &capRegions,
                       std::vector<ProteinRegion> &antiCapRegions,
                       std::vector<ProteinRegion> &nucRegions,
                       std::vector<ProteinRegion> &sevRegions,
                       std::vector<ExcZone> &excZones,
                       const bool activeBranch,
                       const double activeBranchHeight,
                       double arpConc,
                       const bool activeNuc, const double activeNucHeight,
                       const bool activeCap, const double activeCapHeight,
                       const bool activeAntiCap,
                       const double activeAntiCapHeight,
                       const bool activeSever,
                       const double activeSeverHeight);

    void initCortexRegions(std::vector<ProteinRegion> &branchRegions,
                           std::vector<ProteinRegion> &nucRegions,
                           double arpConc);

    void setCOM();
    std::array<double,2> calcCOM();

    void exocytosis(double dt, StericGrid &stericGrid, double temp,
                    double k_exo, double exo_mean, double exo_stdev);


    void forceDepol(int subid, std::vector<Actin> &actinVec,
                    const bool steric, StericGrid &stericGrid,
                    int &nActin, const bool tether,
                    bool arpPool, ArpGrid &arpGrid,
                    std::vector<ExcZone> &excZones,
                    std::vector<MembraneWall> &memWalls,
                    std::vector<Membrane> &membranes,
                    Cortex &cortex,
                    const bool bDynamics,
                    std::vector<ProteinRegion> &branchRegions,
                    std::vector<ProteinRegion> &nucRegions,
                    std::vector<ProteinRegion> &capRegions,
                    std::vector<ProteinRegion> &sevRegions);

    void forceDissociate(std::vector<Actin> &actinVec,
                         const bool steric, StericGrid &stericGrid,
                         int &nActin, const bool tether,
                         bool arpPool, ArpGrid &arpGrid,
                         std::vector<ExcZone> &excZones,
                         std::vector<MembraneWall> &memWalls,
                         std::vector<Membrane> &membranes,
                         Cortex &cortex,
                         const bool bDynamics,
                         std::vector<ProteinRegion> &branchRegions,
                         std::vector<ProteinRegion> &nucRegions,
                         std::vector<ProteinRegion> &capRegions,
                         std::vector<ProteinRegion> &sevRegions);

    void createTetherPointsCortex(Membrane &membrane);
    Eigen::VectorXd calcTetherForces(Cortex &cortex);
    void setTetherSubIDs(std::vector<int> subids) { m_tetherSubIDs = subids; }
    std::vector<int> getTetherSubIDs() const { return m_tetherSubIDs; }
    void eraseTether(int pos)
    {
        m_tetherSubIDs.erase(m_tetherSubIDs.begin()+pos);
        m_tetherRelDists.erase(m_tetherRelDists.begin()+pos);
    }
    void setTetherDists(std::vector<double> relDists) { m_tetherRelDists = relDists; }
    std::vector<double> getTetherDists() const { return m_tetherRelDists; }


    Eigen::VectorXd calcBSForces();
    void moveMem(std::array<double,2> dR, const bool steric,
                           StericGrid &stericGrid);

};

#endif
