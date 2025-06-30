




#ifndef STERICGRID_H
#define STERICGRID_H

typedef gte::FIQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> Find_I;

class ExcZone;
class Membrane;
class Cortex;
class MembraneWall;
class StericGrid
{
private:
    double m_min_xy;
    double m_max_xy;
    double m_cellSize;
    int m_numCellsAcross;
    // m_cellContents[0] is the actin id (index in actin vector), m_cellContents[1] is the subunit id on that actin
    std::vector< std::vector< std::array<int,2> > > m_cellContents;
    std::vector< std::vector< std::array<int,2> > > m_cellContentsMembrane;
    std::vector< std::vector< std::array<int,2> > > m_cellContentsCortex;
    std::vector< std::vector<int> > m_cellContentsMemWall;

public:
    StericGrid();

    StericGrid(double min_xy, double max_xy);

    void updateCellsBarbPoly(Actin &filament);
    void updateCellsPointPoly(Actin &filament);
    void updateCellsAddPointBarb(Actin &filament);
    void updateCellsAddPointPoint(Actin &filament);
    void updateCellsRemPointBarb(Actin &filament);
    void updateCellsRemPointPoint(Actin &filament);
    void updateCellsAll(Actin &filament);
    void moveCellsUp(Actin &filament);
    void moveCellsDown(Actin &filament);
    void dissociate(Actin &filament);
    void resetAndUpdateAllDaughters(Actin &filament, std::vector<Actin> &actinVec);
    void resetAndUpdateAllCLinks(Actin &filament, std::vector<Actin> &actinVec);
    void resetSubAndUpdateAllDaughters(Actin &filament, int p, std::vector<Actin> &actinVec);
    void resetAndUpdateAllCLinksAndDaughters(Actin &filament,
                                             std::vector<Actin> &actinVec);

    void resetAndUpdateAllCLinksAndDaughters(Actin &filament,
                                             std::vector<Actin> &actinVec,
                                             std::vector<int> &calledFilas);

    void resetCells(Actin &filament, int subid);
    void resetSubandUpdate(Actin &filament, int subid);
    void checkSubandUpdate(Actin &filament, int subid);

    int getNumCellsAcross() const { return m_numCellsAcross; }
    std::vector<std::array<int,2>> getCellContents(int cellID) const { return m_cellContents[cellID]; }
    std::vector<std::array<int,2>> getCellContentsMem(int cellID) const { return m_cellContentsMembrane[cellID]; }
    std::vector<std::array<int,2>> getCellContentsCortex(int cellID) const { return m_cellContentsCortex[cellID]; }
    std::vector<int> getCellContentsMemWall(int cellID) const { return m_cellContentsMemWall[cellID]; }

    std::vector<int> getCellsToCheck(std::vector<int> occCells);
    std::vector<int> getCellsToCheckBIG(std::vector<int> occCells, double extraDist);
    std::vector<int> getCellsToCheckMemWall(std::vector<int> occCells);

    void resetMemandUpdate(Membrane &memFila, int subid);
    void resetMem(Membrane &memFila, int subid);
    void updateMembrane(Membrane &memFila, int subid);

    void resetCortexandUpdate(Cortex &cortex, int subid);
    void resetCortex(Cortex &cortex, int subid);
    void updateCortex(Cortex &cortex, int subid);

    void resetMemWallandUpdate(MembraneWall &memWall);
    void resetMemWall(MembraneWall &memWall);
    void updateMemWall(MembraneWall &memWall);


    std::vector<int> findOccCells(int row1, int col1, int row2, int col2,
                                double p1X, double p1Y, double p2X, double p2Y);

    std::vector<int> findOccCellsMemWall(int row1, int col1, int row2, int col2);

    int getCellFromCoord(std::array<double,2> coord);


};

#endif
