/*
 *  StericGrid.cpp
 *
 *  C++ file containing the definition of the StericGrid class.
 *  This contains member variable initialisations and non-trivial member
 *  functions.
 *  The declaration of this class (including the definitions of the trivial
 *  member functions) is found in the associated header file "StericGrid.h"
 *
 *  James Bradford
 *  University of Sheffield
 *  Dec 2018
 */


#include "configHeader.h"
#include "geometry.h"
#include "Membrane.h"
#include "Cortex.h"
#include "StericGrid.h"
#include "Actin.h"
#include "globals.h"
#include "RNG.h"
#include <unordered_set>

extern RNG rng;

typedef gte::TIQuery<double, gte::Segment<2, double>, gte::Segment<2, double> > Check_I;

StericGrid::StericGrid()
{}

StericGrid::StericGrid(double min_xy, double max_xy)
    : m_min_xy { min_xy },
    m_max_xy { max_xy }
{
    // Cell size should be diameter of an actin filament!


    m_cellSize = Actin::s_stericRadius*2;

    double width = m_max_xy - m_min_xy;
    m_numCellsAcross = floor(width/m_cellSize);
    m_cellSize = width/m_numCellsAcross;

    std::cout << "StericGrid cell size: " <<  m_cellSize << std::endl;
    for (int i = 0; i < m_numCellsAcross*m_numCellsAcross; ++i)
    {
        // Loop over every cell, and the (outside) cell 0
        std::vector< std::array<int,2> > empty;
        m_cellContents.push_back(empty);
        m_cellContentsMembrane.push_back(empty);
        m_cellContentsCortex.push_back(empty);
        std::vector<int> emptyInt;

        m_cellContentsMemWall.push_back(emptyInt);
    }

    // m_cellContents will contain the actin ids in each cell
}

void StericGrid::updateCellsBarbPoly(Actin &filament)
{

    // Called after barbed end polymerisation

    // DO on a subunit level

    // new end point

    if (filament.getLength() < 3*Actin::s_segmentationLength)
    {
        for (int i = 0; i < filament.getNumSubs(); ++i)
        {
            resetSubandUpdate(filament, i);
        }
    }
    else
    {
        int endSub = filament.getNumSubs()-1;
        resetSubandUpdate(filament, endSub);
        resetSubandUpdate(filament, endSub-1);
    }
}

void StericGrid::updateCellsPointPoly(Actin &filament)
{

    // Called after pointed end polymerisation

    // DO on a subunit level

    // new end point

    if (filament.getLength() < 3*Actin::s_segmentationLength)
    {
        for (int i = 0; i < filament.getNumSubs(); ++i)
        {
            resetSubandUpdate(filament, i);
        }
    }
    else
    {
        resetSubandUpdate(filament, 0);
        resetSubandUpdate(filament, 1);
    }
}

void StericGrid::updateCellsAddPointBarb(Actin &filament)
{
    // Called after adding a point
    int newEnd = filament.getNumSubs()-1;
    checkSubandUpdate(filament,newEnd);
    resetSubandUpdate(filament, newEnd-1);
    resetSubandUpdate(filament, newEnd-2);

}

void StericGrid::updateCellsAddPointPoint(Actin &filament)
{
    // Called after adding a point
    checkSubandUpdate(filament,0);
    resetSubandUpdate(filament, 1);
    resetSubandUpdate(filament, 2);
}

void StericGrid::updateCellsAll(Actin &filament)
{
    for (int i = 0; i < filament.getNumSubs(); ++i)
    {
        resetSubandUpdate(filament, i);
    }
}
void StericGrid::updateCellsRemPointBarb(Actin &filament)
{
    // Called after removing a point
    int newEnd = filament.getNumSubs()-1;
    resetSubandUpdate(filament, newEnd);
    resetSubandUpdate(filament, newEnd-1);

}

void StericGrid::updateCellsRemPointPoint(Actin &filament)
{
    // Called after removing a point
    resetSubandUpdate(filament, 0);
    resetSubandUpdate(filament, 1);

}

void StericGrid::moveCellsUp(Actin &filament)
{
    // Given adding a pointed end point, need to account for subids incrementing
    // by one
    int id = filament.getID();
    std::vector<int> cellIDs;
    for (int subid = 0; subid < filament.getNumSubs()-1; ++subid)
    {
        std::vector<int> stericCells = filament.getStericCells(subid);

        for (unsigned int i = 0; i < stericCells.size(); ++i)
        {
            int cellID = stericCells[i];
            // Check to see if we've already looked in this cell
            if (std::find(cellIDs.begin(), cellIDs.end(),cellID)==cellIDs.end())
            {
                cellIDs.push_back(cellID);

                for (unsigned int j = 0; j < m_cellContents[cellID].size(); ++j)
                {
                    if (m_cellContents[cellID][j][0] == id)
                    {
                        m_cellContents[cellID][j][1] += 1;
                    }
                }
            }
        }
    }
}

void StericGrid::moveCellsDown(Actin &filament)
{
    // Given removing a pointed end point, need to account for subids decrementing
    // by one
    int id = filament.getID();
    std::vector<int> cellIDs;
    for (int subid = 0; subid < filament.getNumSubs(); ++subid)
    {
        std::vector<int> stericCells = filament.getStericCells(subid);

        for (unsigned int i = 0; i < stericCells.size(); ++i)
        {
            int cellID = stericCells[i];
            // Check to see if we've already looked in this cell
            if (std::find(cellIDs.begin(), cellIDs.end(),cellID)==cellIDs.end())
            {
                cellIDs.push_back(cellID);

                for (unsigned int j = 0; j < m_cellContents[cellID].size(); ++j)
                {
                    if (m_cellContents[cellID][j][0] == id)
                    {
                        m_cellContents[cellID][j][1] -= 1;
                    }
                }
            }
        }
    }
}

void StericGrid::resetCells(Actin &filament, int subid)
{
    // Reset the trackers for this particular subunit only!
    std::vector<int> stericCells = filament.getStericCells(subid);
    filament.clearStericCellsSub(subid);
    int filaID = filament.getID();

    for (unsigned int i = 0; i < stericCells.size(); ++i)
    {
        int cellID = stericCells[i];
        for (unsigned int j = 0; j < m_cellContents[cellID].size(); ++j)
        {
            // Remove it, careful about changing what you are iterating over
            if (m_cellContents[cellID][j][0] == filaID && m_cellContents[cellID][j][1] == subid)
            {
                // Should be removed
                m_cellContents[cellID].erase(m_cellContents[cellID].begin() + (j));
                break;
            }

        }
    }
}

void StericGrid::dissociate(Actin &filament)
{
    // Called after dissociation on every filament with id number greater than that
    // which dissociated
    int id = filament.getID() + 1; // old id
    std::vector<int> cellIDs;
    for (int subid = 0; subid < filament.getNumSubs(); ++subid)
    {
        std::vector<int> stericCells = filament.getStericCells(subid);

        for (unsigned int i = 0; i < stericCells.size(); ++i)
        {
            int cellID = stericCells[i];
            // Check to see if we've already looked in this cell
            if (std::find(cellIDs.begin(), cellIDs.end(),cellID)==cellIDs.end())
            {
                cellIDs.push_back(cellID);

                for (unsigned int j = 0; j < m_cellContents[cellID].size(); ++j)
                {
                    if (m_cellContents[cellID][j][0] == id)
                    {
                        m_cellContents[cellID][j][0] -= 1;
                    }
                }
            }
        }
    }
}


void StericGrid::resetAndUpdateAllDaughters(Actin &filament, std::vector<Actin> &actinVec)
{
    for (int i = 0; i < filament.getNumSubs(); ++i)
    {
        // update the steric grid for the test seed
        resetSubandUpdate(filament, i);
        // Any small branches that have moved need to be updated too
        for (int j : filament.getBranchIDSubVector(i))
        {
            // only consider short branches as longer branches will have not moved
            if (actinVec[j].getLength() < 2*Actin::s_segmentationLength)
            {
                resetAndUpdateAllDaughters(actinVec[j], actinVec);
            }
        }
    }
}

void StericGrid::resetAndUpdateAllCLinks(Actin &filament, std::vector<Actin> &actinVec)
{
    for (int i = 0; i < filament.getNumCLinks(); ++i)
    {
        int cLinkID = filament.getCLinkActinAndSite(i)[2];
        if (actinVec[cLinkID].getLength() < 2*Actin::s_segmentationLength)
        {
            for (int j = 0; j < actinVec[cLinkID].getNumSubs(); ++j)
            {
                resetSubandUpdate(actinVec[cLinkID], j);
            }

        }
    }
    // Check branches
    for (int i = 0; i < filament.getNumSubs(); ++i)
    {
        for (int j : filament.getBranchIDSubVector(i))
        {
            // only consider short branches as longer branches will have not moved
            if (actinVec[j].getLength() < 2*Actin::s_segmentationLength)
            {
                resetAndUpdateAllCLinks(actinVec[j], actinVec);
            }
        }
    }

}

void StericGrid::resetAndUpdateAllCLinksAndDaughters(Actin &filament, std::vector<Actin> &actinVec)
{
    // Going to pass through a vector containing ids that the function has already been called on
    std::vector<int> calledFilas;
    resetAndUpdateAllCLinksAndDaughters(filament, actinVec, calledFilas);
}

void StericGrid::resetAndUpdateAllCLinksAndDaughters(Actin &filament,
                                                     std::vector<Actin> &actinVec,
                                                     std::vector<int> &calledFilas)
{
    // first check if this function has been called on filament before
    if (std::find(calledFilas.begin(), calledFilas.end(), filament.getID()) != calledFilas.end())
    {
            // We have done this filament already, so return
            return;
    }

    // Then add "filament" to calledFilas
    calledFilas.push_back(filament.getID());

    for (int i = 0; i < filament.getNumSubs(); ++i)
    {
        // update the self filament
        resetSubandUpdate(filament, i);

        // Any small branches that have moved need to be updated too
        for (int j : filament.getBranchIDSubVector(i))
        {
            // only consider short branches as longer branches will have not moved
            if (actinVec[j].getLength() < 2*Actin::s_segmentationLength)
            {
                resetAndUpdateAllCLinksAndDaughters(actinVec[j], actinVec, calledFilas);
            }
        }
    }

    // Now look at clinks
    for (int i = 0; i < filament.getNumCLinks(); ++i)
    {
        // update the steric grid for the test seed
        int cLinkID = filament.getCLinkActinAndSite(i)[2];
        if (actinVec[cLinkID].getLength() < 2*Actin::s_segmentationLength)
        {
            resetAndUpdateAllCLinksAndDaughters(actinVec[cLinkID], actinVec, calledFilas);

        }
    }

}

void StericGrid::resetSubAndUpdateAllDaughters(Actin &filament, int p, std::vector<Actin> &actinVec)
{
    // for sub p and all daughters on that sub

    // update the steric grid for the test seed
    resetSubandUpdate(filament, p);
    // Any small branches that have moved need to be updated too
    for (int j : filament.getBranchIDSubVector(p))
    {
        // only consider short branches as longer branches will have not moved
        if (actinVec[j].getLength() < 2*Actin::s_segmentationLength)
        {
            resetAndUpdateAllDaughters(actinVec[j], actinVec);
        }
    }

}

void StericGrid::resetSubandUpdate(Actin &filament, int subid)
{

    resetCells(filament, subid);
    checkSubandUpdate(filament, subid);
}

void StericGrid::checkSubandUpdate(Actin &filament, int subid)
{

    // Refind all the cells that this sub belongs to
    //std::cout << "Testing placement cells across: " << m_numCellsAcross << std::endl;


    int col1 = floor((filament.getPoints()[subid][0] - m_min_xy)/m_cellSize);
    // col1 = 0 for left
    int row1 = floor((filament.getPoints()[subid][1] - m_min_xy)/m_cellSize);
    // row1 = 0 for bottom
    if (col1 < 0 || row1 < 0 || col1 >= m_numCellsAcross || row1 >= m_numCellsAcross)
    {
        // Outside the grid
        std::cout << "A filament is outside the grid (1), ending sim" << std::endl;
        std::cout << col1 << ", " << row1 << std::endl;
        std::cout << filament.getPoints()[subid][0] << " : " << subid << std::endl;
        exit(0);
    }
    int col2 = floor((filament.getPoints()[subid+1][0] - m_min_xy)/m_cellSize);
    int row2 = floor((filament.getPoints()[subid+1][1] - m_min_xy)/m_cellSize);

    if (col2 < 0 || row2 < 0 || col2 >= m_numCellsAcross || row2 >= m_numCellsAcross)
    {
        // Outside the grid
        std::cout << "A filament is outside the grid (2), ending sim" << std::endl;
        std::cout << col2 << ", " << row2 << std::endl;
        exit(0);
    }

    std::vector<int> occCells = findOccCells(row1, col1, row2, col2, filament.getPoints()[subid][0], filament.getPoints()[subid][1], filament.getPoints()[subid+1][0], filament.getPoints()[subid+1][1]);
    std::array<int,2> iDAndSubID = { filament.getID(), subid };
    for (unsigned int j = 0; j < occCells.size(); ++j)
    {
        m_cellContents[occCells[j]].push_back(iDAndSubID);
        filament.addCellToActinSub(subid, occCells[j]);
    }


}

void StericGrid::resetMemandUpdate(Membrane &membrane, int subid)
{

    resetMem(membrane, subid);
    updateMembrane(membrane, subid);
}


void StericGrid::resetMem(Membrane &membrane, int subid)
{
    // Reset the trackers for this particular subunit only!
    std::vector<int> stericCells = membrane.getStericCells(subid);
    membrane.clearStericCellsSub(subid);

    for (unsigned int i = 0; i < stericCells.size(); ++i)
    {
        int cellID = stericCells[i];
        for (unsigned int j = 0; j < m_cellContentsMembrane[cellID].size(); ++j)
        {
            // Remove it, careful about changing what you are iterating over

            if (m_cellContentsMembrane[cellID][j][0] == membrane.getID() && m_cellContentsMembrane[cellID][j][1] == subid)
            {
                // Should be removed
                m_cellContentsMembrane[cellID].erase(m_cellContentsMembrane[cellID].begin() + (j));
                break;
            }

        }
    }
}


std::vector<int> StericGrid::getCellsToCheck(std::vector<int> occCells)
{

    // Changed so will have repeats here, this *might* be more efficient than
    // Checking

    std::vector<int> cellsToCheck;
    cellsToCheck.reserve(9);
    //std::array<int,9> cellsToCheck;

    for (unsigned int i = 0; i < occCells.size(); ++i)
    {
        int cell = occCells[i];
        assert(cell > 0);

        // Bot row
        if (cell >= m_numCellsAcross)
        {
            // left
            if (cell % m_numCellsAcross != 0)
            {
                cellsToCheck.push_back(cell-1-m_numCellsAcross);
            }
            // mid
            cellsToCheck.push_back(cell-m_numCellsAcross);

            // right
            if (cell % m_numCellsAcross != (m_numCellsAcross-1))
            {
                cellsToCheck.push_back(cell+1-m_numCellsAcross);
            }
        }

        // Mid row
        // left
        if (cell % m_numCellsAcross != 0)
        {
            cellsToCheck.push_back(cell-1);
        }

        // mid
        cellsToCheck.push_back(cell);

        // right
        if (cell % m_numCellsAcross != (m_numCellsAcross-1))
        {
            cellsToCheck.push_back(cell+1);
        }

        // Top row
        if (cell < (m_numCellsAcross*m_numCellsAcross)-m_numCellsAcross) // If the cell is not on the top row
        {
            // left
            if (cell % m_numCellsAcross != 0) // If the cell is not on the leftmost column
            {
                // Cell to the top left of the occupied cell
                cellsToCheck.push_back(cell-1+m_numCellsAcross);
            }

            // mid
            cellsToCheck.push_back(cell+m_numCellsAcross);

            // right
            if (cell % m_numCellsAcross != (m_numCellsAcross-1))
            {
                cellsToCheck.push_back(cell+1+m_numCellsAcross);
            }
        }
    }

    return cellsToCheck;
}

std::vector<int> StericGrid::getCellsToCheckBIG(std::vector<int> occCells,
                                                double extraDist)
{
    // For the cortex! - which is thicker than our actin filaments and so
    // requires checking over more than 9 cells

    // Changed so will have repeats here, this *might* be more efficient than
    // Checking

    // extraDist is cortexThickness or crosslink distance

    // Now will have more cells to check since bigger membrane or whatever
    double sizeToCheck = extraDist + 6*Actin::s_stericRadius;
    assert(sizeToCheck < m_max_xy - m_min_xy);
    int numCells = m_cellContents.size();
    int widthNumCells = ceil(sizeToCheck/m_cellSize); // Has to be an odd number
    if (widthNumCells % 2 == 0)
    {
        // Even
        widthNumCells += 1;
    }

    std::vector<int> cellsToCheck;
    for (unsigned int i = 0; i < occCells.size(); ++i)
    {
        int cell = occCells[i];
        assert(cell > 0);
        // start bot left
        for (int relRow = -(widthNumCells-1)/2; relRow <= (widthNumCells-1)/2; ++relRow) // row relative to cell
        {
            for (int relCol = -(widthNumCells-1)/2; relCol <= (widthNumCells-1)/2; ++relCol) // col relative to cell
            {
                // Bounds check
                if ((cell + relCol + m_numCellsAcross*relRow) >= 0 && (cell + relCol + m_numCellsAcross*relRow) < numCells )
                {
                    cellsToCheck.push_back(cell + relCol + m_numCellsAcross*relRow);
                }
                else
                {
                  std::cout << "Grid out of bounds?" << std::endl;
                }
            }
        }

    }

    return cellsToCheck;
}

std::vector<int> StericGrid::getCellsToCheckMemWall(std::vector<int> occCells)
{

    // Changed so will have repeats here, this *might* be more efficient than
    // Checking

    std::vector<int> cellsToCheck;
    cellsToCheck.reserve(2);

    for (unsigned int i = 0; i < occCells.size(); ++i)
    {
        int cell = occCells[i];
        assert(cell > 0);

        // Bot row
        if (cell >= m_numCellsAcross)
        {
            // mid
            cellsToCheck.push_back(cell-m_numCellsAcross);
        }

        // Mid row

        // mid
        cellsToCheck.push_back(cell);


    }

    return cellsToCheck;
}

void StericGrid::updateMembrane(Membrane &membrane, int subid)
{

    // Refind all the cells that this sub belongs to

    int col1 = floor((membrane.getPoints()[subid][0] - m_min_xy)/m_cellSize);
    int row1 = floor((membrane.getPoints()[subid][1] - m_min_xy)/m_cellSize);
    if (col1 < 0 || row1 < 0 || col1 >= m_numCellsAcross || row1 >= m_numCellsAcross)
    {
        // Outside the grid
        std::cout << "A part of membrane is outside the grid 1, ending sim" << std::endl;
        std::cout << membrane.getPoints()[subid][0] << ", " << membrane.getPoints()[subid][1] << std::endl;
        exit(0);
    }

    int np = subid + 1; // next point
    if (subid == membrane.getNumPoints()-1)
    {
        // end sub
        np = 0;
    }
    int col2 = floor((membrane.getPoints()[np][0] - m_min_xy)/m_cellSize);
    int row2 = floor((membrane.getPoints()[np][1] - m_min_xy)/m_cellSize);


    if (col2 < 0 || row2 < 0 || col2 >= m_numCellsAcross || row2 >= m_numCellsAcross)
    {
        // Outside the grid
        std::cout << "A part of membrane is outside the grid 2, ending sim" << std::endl;
        std::cout << membrane.getPoints()[np][0] << ", " << membrane.getPoints()[np][1] << std::endl;
        exit(0);
    }

    std::vector<int> occCells = findOccCells(row1, col1, row2, col2, membrane.getPoints()[subid][0], membrane.getPoints()[subid][1], membrane.getPoints()[np][0], membrane.getPoints()[np][1]);
    std::array<int,2> iDAndSubID = { membrane.getID(), subid };
    for (unsigned int j = 0; j < occCells.size(); ++j)
    {
        m_cellContentsMembrane[occCells[j]].push_back(iDAndSubID);
        membrane.addCellToMemSub(subid, occCells[j]);
    }

}

void StericGrid::resetCortexandUpdate(Cortex &cortex, int subid)
{

    resetCortex(cortex, subid);
    updateCortex(cortex, subid);
}

void StericGrid::resetCortex(Cortex &cortex, int subid)
{
    // Reset the trackers for this particular subunit only!
    std::vector<int> stericCells = cortex.getStericCells(subid);
    cortex.clearStericCellsSub(subid);

    for (unsigned int i = 0; i < stericCells.size(); ++i)
    {
        int cellID = stericCells[i];
        for (unsigned int j = 0; j < m_cellContentsCortex[cellID].size(); ++j)
        {
            // Remove it, careful about changing what you are iterating over

            if (m_cellContentsCortex[cellID][j][0] == cortex.getID() && m_cellContentsCortex[cellID][j][1] == subid)
            {
                // Should be removed
                m_cellContentsCortex[cellID].erase(m_cellContentsCortex[cellID].begin() + (j));
                break;
            }

        }
    }
}

void StericGrid::updateCortex(Cortex &cortex, int subid)
{

    // Refind all the cells that this sub belongs to

    int col1 = floor((cortex.getPoints()[subid][0] - m_min_xy)/m_cellSize);
    int row1 = floor((cortex.getPoints()[subid][1] - m_min_xy)/m_cellSize);
    if (col1 < 0 || row1 < 0 || col1 >= m_numCellsAcross || row1 >= m_numCellsAcross)
    {
        // Outside the grid
        std::cout << "A part of cortex is outside the grid 1, ending sim" << std::endl;
        std::cout << cortex.getPoints()[subid][0] << ", " << cortex.getPoints()[subid][1] << std::endl;
        exit(0);
    }

    int np = subid + 1; // next point
    if (subid == cortex.getNumPoints()-1)
    {
        // End sub
        np = 0;
    }
    int col2 = floor((cortex.getPoints()[np][0] - m_min_xy)/m_cellSize);
    int row2 = floor((cortex.getPoints()[np][1] - m_min_xy)/m_cellSize);


    if (col2 < 0 || row2 < 0 || col2 >= m_numCellsAcross || row2 >= m_numCellsAcross)
    {
        // Outside the grid
        std::cout << "A part of cortex is outside the grid 2, ending sim" << std::endl;
        std::cout << cortex.getPoints()[np][0] << ", " << cortex.getPoints()[np][1] << std::endl;
        exit(0);
    }

    std::vector<int> occCells = findOccCells(row1, col1, row2, col2, cortex.getPoints()[subid][0], cortex.getPoints()[subid][1], cortex.getPoints()[np][0], cortex.getPoints()[np][1]);
    std::array<int,2> iDAndSubID = { cortex.getID(), subid };
    for (unsigned int j = 0; j < occCells.size(); ++j)
    {
        m_cellContentsCortex[occCells[j]].push_back(iDAndSubID);
        cortex.addCellToCortexSub(subid, occCells[j]);
    }
}


void StericGrid::resetMemWallandUpdate(MembraneWall &memWall)
{

    resetMemWall(memWall);
    updateMemWall(memWall);
}

void StericGrid::resetMemWall(MembraneWall &memWall)
{
    // Reset the trackers for this particular subunit only!
    std::vector<int> stericCells = memWall.getStericCells();
    memWall.clearStericCells();

    for (unsigned int i = 0; i < stericCells.size(); ++i)
    {
        int cellID = stericCells[i];
        for (unsigned int j = 0; j < m_cellContentsMemWall[cellID].size(); ++j)
        {
            // Remove it, careful about changing what you are iterating over

            if (m_cellContentsMemWall[cellID][j] == memWall.getID())
            {
                // Should be removed
                m_cellContentsMemWall[cellID].erase(m_cellContentsMemWall[cellID].begin() + (j));
                break;
            }

        }
    }
}

void StericGrid::updateMemWall(MembraneWall &memWall)
{

    // Refind all the cells that this sub belongs to

    int col1 = floor(((-memWall.getXlength()/2) - m_min_xy)/m_cellSize);
    int row1 = floor((memWall.getYpos() - m_min_xy)/m_cellSize);
    if (col1 < 0 || row1 < 0 || col1 >= m_numCellsAcross || row1 >= m_numCellsAcross)
    {
        // Outside the grid
        std::cout << "A part of memWall is outside the grid 1, ending sim" << std::endl;
        exit(0);
    }

    int col2 = floor(((memWall.getXlength()/2) - m_min_xy)/m_cellSize);
    int row2 = row1;


    if (col2 < 0 || row2 < 0 || col2 >= m_numCellsAcross || row2 >= m_numCellsAcross)
    {
        // Outside the grid
        std::cout << "A part of memWall is outside the grid 2, ending sim" << std::endl;
        exit(0);
    }

    std::vector<int> occCells = findOccCellsMemWall(row1, col1, row2, col2);
    for (unsigned int j = 0; j < occCells.size(); ++j)
    {
        m_cellContentsMemWall[occCells[j]].push_back(memWall.getID());
        memWall.addCellToMemWall(occCells[j]);
    }
}


std::vector<int> StericGrid::findOccCells(int row1, int col1, int row2, int col2,
                                          double p1X, double p1Y, double p2X, double p2Y)
{
  // Use intersections to see which cells are occupied and which are not
    std::vector<int> occCells;


    gte::Segment<2,double> filament;
    filament.p[1][0] = p2X;
    filament.p[0][0] = p1X;
    filament.p[1][1] = p2Y;
    filament.p[0][1] = p1Y;


    int ystep, xstep;    // the step on y and x axi
    int y = row1, x = col1;  // the line points
    int dx = col2 - col1; // positive
    int dy = row2 - row1; // positive
    int cellID = y*m_numCellsAcross + x;
    occCells.push_back(cellID);
    // NB the last point can't be here, because of its previous point (which has to be verified)
    if (dy < 0)
    {
      ystep = -1;
      dy = -dy; // kept positive
    }
    else
      ystep = 1;

    if (dx < 0)
    {
      xstep = -1;
      dx = -dx; // kept positive
    }
    else
      xstep = 1;

    if (dx == 0 && dy != 0)
    {
        // Just moves in y, no need to check intersections
        for (int i = 0; i < dy; ++i)
        {
          y += ystep;
          cellID = y*m_numCellsAcross + x;
          occCells.push_back(cellID);
        }

        assert ((y == row2) && (x == col2));  // the last point (row2,col2) has to be the same with the last point of the algorithm
        return occCells;
    }
    else if (dy == 0 && dx != 0)
    {
      // Just moves in x, no need to check intersections
      for (int i = 0; i < dx; ++i)
      {
        x += xstep;
        cellID = y*m_numCellsAcross + x;
        occCells.push_back(cellID);
      }

      assert ((y == row2) && (x == col2));  // the last point (row2,col2) has to be the same with the last point of the algorithm
      return occCells;
    }
    else if (dy == 0 && dx == 0)
    {
        // Both in same cell, already have it
        assert ((row1 == row2) && (col1 == col2));
        return occCells;
    }

    for (int i = 0 ; i < dx+dy ; i++)
    {  // do not use the first point (already done)

      // check intersection of either right or left grid line
      if (xstep > 0)
      {
        // check right
        // x is column, y is row
        gte::Segment<2,double> gridLine;
        gridLine.p[1][0] = m_min_xy + m_cellSize*(x+1); // right
        gridLine.p[0][0] = m_min_xy + m_cellSize*(x+1); // right
        gridLine.p[1][1] = m_min_xy + m_cellSize*(y+1); // top
        gridLine.p[0][1] = m_min_xy + m_cellSize*y; // bot
        Check_I intersect;
        auto result = intersect(gridLine, filament);


        if (result.intersect)
        {
            x += xstep;
        }
        else
        {
            y += ystep;
        }

      }
      else
      {
          // check left
          gte::Segment<2,double> gridLine;
          gridLine.p[1][0] = m_min_xy + m_cellSize*x; // left
          gridLine.p[0][0] = m_min_xy + m_cellSize*x; // left
          gridLine.p[1][1] = m_min_xy + m_cellSize*(y+1); // top
          gridLine.p[0][1] = m_min_xy + m_cellSize*y; // bot
          Check_I intersect;
          auto result = intersect(gridLine, filament);

          if (result.intersect)
          {
              x += xstep;
          }
          else
          {
              y += ystep;
          }

      }
      cellID = (y)*m_numCellsAcross + (x);
      occCells.push_back(cellID);
    }

    assert ((y == row2) && (x == col2));  // the last point (row2,col2) has to be the same with the last point of the algorithm

    return occCells;
}


std::vector<int> StericGrid::findOccCellsMemWall(int row1, int col1, int row2, int col2)
{
  // Use intersections to see which cells are occupied and which are not
    std::vector<int> occCells; // JAMES: in cellid form

    int initCellID = row1*m_numCellsAcross + col1;
    int endCellID = row2*m_numCellsAcross + col2;
    int dCell = endCellID - initCellID;

    for (int i = 0; i <= dCell; ++i)
    {
      occCells.push_back(initCellID + i);
    }

    return occCells;
}

int StericGrid::getCellFromCoord(std::array<double,2> coord)
{
    // Given a 2d coordinate, return the cell it is in

    int col = floor((coord[0] - m_min_xy)/m_cellSize);
    int row = floor((coord[1] - m_min_xy)/m_cellSize);

    if (col < 0 || row < 0 || col >= m_numCellsAcross || row >= m_numCellsAcross)
    {
        // Outside the grid
        std::cout << "Coordinate point is outside the steric Grid" << std::endl;
        exit(0);
    }

    return (row*m_numCellsAcross + col);

}
