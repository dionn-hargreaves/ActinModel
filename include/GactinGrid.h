




#ifndef GACTINGRID_H
#define GACTINGRID_H

typedef gte::FIQuery<double, gte::Segment<2, double>, gte::Segment<2, double>> Find_I;

class GactinGrid
{
private:
    bool m_exist;
    double m_x_min;
    double m_y_min;
    double m_x_max;
    double m_y_max;
    double m_x_width;
    double m_y_height;
    double m_cellWidth;
    double m_cellHeight;
    double m_cellDepth;
    double m_convert;
    int m_grids_y;
    int m_grids_x;
    bool m_flow;
    std::vector< std::vector<double> > m_conc_vec;
    std::vector< std::vector<int>  > m_num_vec;

    // For latrunculin A
    std::vector< std::vector<double> > m_LatA_conc_vec;
    std::vector< std::vector<int>  > m_LatA_num_vec;
    std::vector< std::vector<int>  > m_boundLatA_num_vec;
public:
    GactinGrid();

    GactinGrid(double x_min, double y_min, double x_max, double y_max,
               double cellWidth, double cellHeight, double cellDepth,
               double gActinConc, double latAConc = 0, bool flow = false);


    std::array<int, 2> findGrid(double xpos, double ypos);
    void decrement(double xpos, double ypos);
    void decrementPyrene(double xpos, double ypos);
    void increment(double xpos, double ypos);
    void incrementPyrene(double xpos, double ypos);
    void dissociate(double xpos, double ypos);
    void dissociateBR(double xpos, double ypos);
    void nucleation(double xpos, double ypos);
    void branching(double xpos, double ypos);
    void branchingCell(int cellX, int cellY);
    double getConcentration(int xcoord, int ycoord) const;
    int getNumber(int xcoord, int ycoord) const;
    int getPyreneNumber(int xcoord, int ycoord) const;
    int getLatANumber(int xcoord, int ycoord) const;
    int getLatActinNumber(int xcoord, int ycoord) const;
    std::array<int,2> findNearestCell(double xpos, double ypos);
    void calcFlow(double dt);

    bool getExist() const { return m_exist; }
    bool getFlow() const { return m_flow; }
    int getNumCellsX() const { return m_grids_x; }
    int getNumCellsY() const { return m_grids_y; }

    double getCellWidth() const { return m_cellWidth; }
    double getCellHeight() const { return m_cellHeight; }
    double getCellDepth() const { return m_cellDepth; }

    void changeRegion(double newXmin, double newYmin, double newXmax,
                      double newYmax, double newCellDepth);

    void latBindingMC(double dt);
    void latBinding(int cellX, int cellY, int amount);
    void latUnbindingMC(double dt);



};

#endif
