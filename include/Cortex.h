

#ifndef CORTEX_H
#define CORTEX_H

class Cortex
{
private:
    bool m_exist;
    int m_id;
    double m_ypos;
    double m_length;
    int m_numPoints;
    double m_bendingMod;
    double m_segLength;
    double m_thickness;
    double m_zWidth;
    double m_effPersLen; // effective persistence length Lp = K_bend*2R / kbT

    double m_tetherStiff; // for tethering of cortex to membrane
    std::vector<int> m_tetherSubIDs;
    std::vector<double> m_tetherRelDists;

    std::array<double,2> m_COM; // centre of mass

    std::vector< std::array<double,2> > m_points;
    std::vector< std::array<double,2> > m_subunit_unitVecs; // unit vectors of the subunits

    std::vector<double> m_actualSubLengths;
    std::vector<double> m_prescribedSubLengths;

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

    bool check_min_dist_Grid(int subid, const std::vector<Actin> &actinvec,
                                                  std::vector<int> &cellsToCheck,
                                                  StericGrid &stericGrid);

    std::array<int,2> check_min_dist_Grid_Identify(int subid, const std::vector<Actin> &actinvec,
                                                   std::vector<int> &cellsToCheck,
                                                   StericGrid &stericGrid);


public:
    Cortex();

    Cortex(double temp, Membrane &membrane);


    ~Cortex(); // destructor

    bool getExist() const { return m_exist; }
    double getYpos() const { return m_ypos; }
    double getRadius() const { return (m_length/(2*M_PI)); }
    std::vector<double> initCircle();
    void calcSubunitVecs();
    void calcSubLengths();
    int getID() const { return m_id; }

    std::vector< std::array<double,2> > getPoints() const { return m_points; }
    double getNumPoints() const { return m_numPoints; };
    double getLength() const { return m_length; }
    double getThickness() const { return m_thickness; }
    double getSegLength() const { return m_segLength; }
    double get1DBendMod() const { return (m_bendingMod*m_zWidth); }
    double getEffLp() const { return m_effPersLen; }
    std::vector<std::array<double,2>> getUnitVecs() const { return m_subunit_unitVecs; }
    std::vector<double> getActualSubLengths() const { return m_actualSubLengths; }
    std::vector<double> getPresSubLengths() const { return m_prescribedSubLengths; }
    std::array<double,2> getCOM() const { return m_COM; }

    void fluctuate(double temp, double viscosity, double dt,
                   std::vector<Actin> &actinvec,
                   std::vector<ExcZone> &excZones,
                   const bool steric,
                   StericGrid &stericGrid, double currTime,
                   std::vector<ProteinRegion> &branchRegions,
                   bool &branching,
                   std::vector<ProteinRegion> &capRegions,
                   std::vector<ProteinRegion> &antiCapRegions,
                   std::vector<ProteinRegion> &nucRegions,
                   std::vector<Membrane> &membranes,
                   int &nActin, const bool tether,
                   bool arpPool, ArpGrid &arpGrid,
                   std::vector<MembraneWall> &memWalls,
                   std::ofstream &fusionTime,
                   const bool bDynamics);

    static bool s_checkExVolBarbPoly(const Actin &filament, Cortex &cortex,
                                     StericGrid &stericGrid);

    bool checkExVolBarbPoly(const Actin &filament) const;
    bool checkExVolBarbPolyGrid(const Actin &filament,
                                StericGrid &stericGrid) const;

    static bool s_checkExVolPointPoly(const Actin &filament, Cortex &cortex,
                                      StericGrid &stericGrid);
    bool checkExVolPointPoly(const Actin &filament) const;
    bool checkExVolPointPolyGrid(const Actin &filament,
                                 StericGrid &stericGrid) const;

    bool checkExVolNuc(const Actin &actin) const;

    static bool s_checkExVolNuc(const Actin &actin, Cortex &cortex,
                                std::vector<int> cellsToCheck,
                                StericGrid &stericGrid);

    bool checkExVolNuc_Grid(const Actin &actin,
                            std::vector<int> cellsToCheck,
                            StericGrid &stericGrid) const;

    bool checkExVolBM(const Actin &actin) const;
    bool checkInOut(const Actin &actin) const;
    bool checkInOutPoly(gte::Segment<3,double> monomer) const;
    bool checkPointIn(std::array<double,2> point) const;

    std::vector<int> getStericCells(int subunit) const { return m_stericCells[subunit]; }
    void clearStericCellsSub(int subunit) { m_stericCells[subunit].clear(); }
    void addCellToCortexSub(int subunit, int cell) { m_stericCells[subunit].push_back(cell); }

    static bool s_checkExVolBM(int subid, const Actin &actin, Cortex &cortex,
                               std::vector<int> cellsToCheck,
                               StericGrid &stericGrid);

    bool check_min_dist_Actin_Grid(int subid, const Actin &actin,
                                   std::vector<int> cellsToCheck,
                                    StericGrid &stericGrid) const;

    void setCOM();
    std::array<double,2> calcCOM();

    void exocytosis(double dt, StericGrid &stericGrid, double temp);

    void createTetherPointsCortex(Membrane &membrane);
    Eigen::VectorXd calcTetherForces(Membrane &other);
    void setTetherSubIDs(std::vector<int> subids) { m_tetherSubIDs = subids; }
    std::vector<int> getTetherSubIDs() const { return m_tetherSubIDs; }
    void eraseTether(int pos)
    {
        m_tetherSubIDs.erase(m_tetherSubIDs.begin()+pos);
        m_tetherRelDists.erase(m_tetherRelDists.begin()+pos);
    }
    void setTetherDists(std::vector<double> relDists) { m_tetherRelDists = relDists; }
    std::vector<double> getTetherDists() const { return m_tetherRelDists; }

    void moveCortex(std::array<double,2> dR, const bool steric,
                    StericGrid &stericGrid);

};

#endif
