

#ifndef MEMBRANEWALL_H
#define MEMBRANEWALL_H



class MembraneWall
{
private:
    int m_id;
    bool m_exist;
    double m_ypos;
    double m_length;
    double m_extForce;
    double m_origYPos;

    std::vector<int> m_coupledBranchRegionsIDs;
    std::vector<int> m_coupledNucRegionsIDs;
    std::vector<int> m_coupledCapRegionsIDs;
    std::vector<int> m_coupledSevRegionsIDs;

    std::vector<int> m_stericCells;

public:
    MembraneWall(double ypos, double length, double extForce);
    int getID() const { return m_id; }
    double getForce() const { return m_extForce; }
    double getYpos() const { return m_ypos; }
    double getOrigYpos() const { return m_origYPos; }
    double getXlength() const { return m_length; }
    //void moveWall(double deltaD) { m_ypos += deltaD; }
    void moveWall(double deltaD, StericGrid &stericGrid);

    double distToMoveWall_barbed(Actin &actin);
    double distToMoveWall_pointed(Actin &actin);
    double distToMoveWallBack(const std::vector<Actin> &actinvec);
    double distToMoveWallBack2(const std::vector<Actin> &actinvec,
                               StericGrid &stericGrid);

    double calc_min_dist_Grid(const std::vector<Actin> &actinvec,
                              std::vector<int> &cellsToCheck,
                              StericGrid &stericGrid);

    void moveCoupledRegions(double dx, double dy,
                            std::vector<ProteinRegion> &branchRegions,
                            std::vector<ProteinRegion> &nucRegions,
                            std::vector<ProteinRegion> &capRegions,
                            std::vector<ProteinRegion> &sevRegions);

    bool checkExVolNuc(const Actin &actin);
    bool checkExVolBM(const Actin &actin) const;

    void appendToBranchRegions(int id) { m_coupledBranchRegionsIDs.push_back(id); }
    void appendToNucRegions(int id) { m_coupledNucRegionsIDs.push_back(id); }
    void appendToCapRegions(int id) { m_coupledCapRegionsIDs.push_back(id); }
    void appendToSevRegions(int id) { m_coupledSevRegionsIDs.push_back(id); }

    std::vector<int> getStericCells() const { return m_stericCells; }
    void clearStericCells() { m_stericCells.clear(); }
    void addCellToMemWall(int cell) { m_stericCells.push_back(cell); }


};

#endif
