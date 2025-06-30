




#ifndef EXCZONE_H
#define EXCZONE_H

#include "configHeader.h"
#include "StericGrid.h"

class Actin;
class Membrane;

class ExcZone
{
private:
    int m_id;
    double m_x;
    double m_y;
    double m_radius;
    double m_x_min;
    double m_y_min;
    double m_x_max;
    double m_y_max;
    bool m_isCircle;
    bool m_isRectangle;
    bool m_motion;
    bool m_directedmotion; // in for eg AFM experiment
    double m_totalMoveDown; // the total movement directed

    //Triangle
    double m_len1; // shared length
    double m_len2; // third length
    std::array<double,2> m_p1; // lowest point
    std::array<double,2> m_p2;
    std::array<double,2> m_p3;


public:
    ExcZone(double x, double y, double radius, int ID, bool motion=false, bool directedmotion=false);
    ExcZone(double x_min, double y_min, double x_max, double y_max, int ID, bool dummy1, bool dummy2);
    ExcZone(double p1X, double p1Y, double len1, double len2, int ID, bool directedmotion, bool dummy1, bool dummy2);

    bool check_ex_vol(const Actin &actin) const;
    bool circ_check_ex_vol(const Actin &actin) const;
    bool circ_check_ex_vol(int subid, const Actin &actin) const;

    bool rect_check_ex_vol(const Actin &actin) const;
    bool rect_check_ex_vol(int subid, const Actin &actin) const;
    bool rect_checkPointwithin(std::array<double,2> Point) const;

    static bool s_checkExVolMemFila(const Membrane &memFila,
                                    const std::vector<ExcZone> &excZones,
                                    int pointID);

    static bool s_checkExVol(int subid, const Actin &actin,
                             const std::vector<ExcZone> excZones);


    bool checkExVolMemFilaCirc(const Membrane &memFila);
    bool checkExVolMemFilaCirc(const std::vector<Membrane> &membranes);
    bool checkExVolMemFilaCircSUB(const Membrane &memFila, int pointID);

    bool checkExVolMemFilaRectSUB(const Membrane &membrane, int pointID);

    bool checkExVolMemFilaTri(const std::vector<Membrane> &membranes);
    bool checkExVolMemFilaTriSUB(const Membrane &membrane, int pointID);

    static bool s_check_polymerisation_barbed(const Actin &actin, const std::vector<ExcZone> &excZones);
    bool circ_check_polymerisation_barbed(const Actin &actin);
    bool rect_check_polymerisation_barbed(const Actin &actin);

    static bool s_check_polymerisation_pointed(const Actin &actin, const std::vector<ExcZone> &excZones);
    bool circ_check_polymerisation_pointed(const Actin &actin);
    bool rect_check_polymerisation_pointed(const Actin &actin);

    static bool s_checkExVolCortex(const Cortex &cortex,
                                   const std::vector<ExcZone> &excZones,
                                   int pointID);

    bool checkExVolCortexCircSUB(const Cortex &cortex, int pointID);
    bool checkExVolCortexRectSUB(const Cortex &cortex, int pointID);
    bool checkExVolCortexTriSUB(const Cortex &cortex, int pointID);


    double getX() const { return m_x; }
    double getY() const { return m_y; }
    double getR() const { return m_radius; }
    double getXmin() const { return m_x_min; }
    double getYmin() const { return m_y_min; }
    double getXmax() const { return m_x_max; }
    double getYmax() const { return m_y_max; }
    bool getCircBool() const { return m_isCircle; }
    bool getRectBool() const { return m_isRectangle; }
    int getID() const { return m_id; }
    bool getMotion() const { return m_motion; }
    bool getDirMove() const { return m_directedmotion; }
    std::array<double,2> getP1() const { return m_p1; }
    std::array<double,2> getP2() const { return m_p2; }
    std::array<double,2> getP3() const { return m_p3; }
    double getTotalMoveDown() const { return m_totalMoveDown; }

    void addX(double dx) { m_x += dx; }
    void addY(double dy) { m_y += dy; }
    void addTriX(double dx) { m_p1[0] += dx; m_p2[0] += dx; m_p3[0] += dx; }
    void addTriY(double dy) { m_p1[1] += dy; m_p2[1] += dy; m_p3[1] += dy; }

    void addTotalMoveDown(double ds) { m_totalMoveDown += ds; }



    bool checkExVolCircTar(const std::vector<ExcZone> excZones);
    void targetBM(double temp, double viscosity, double dt,
                  const std::vector<Membrane> &membranes,
                  std::vector<ExcZone> &excZones,
                  std::vector<ProteinRegion> &branchRegions,
                  std::vector<ProteinRegion> &nucRegions,
                  std::vector<ProteinRegion> &aCapRegions,
                  std::vector<Actin> &actinVec);

    void targetdirectedMovement(double dt,
                                const std::vector<Membrane> &membranes,
                                std::vector<ExcZone> &excZones,
                                std::vector<ProteinRegion> &branchRegions,
                                std::vector<ProteinRegion> &nucRegions,
                                std::vector<ProteinRegion> &aCapRegions);

    static void s_targetdirectedMovementCircANDTri(double dt, bool &firstTouch,
                                                   std::array<double,2> velVector,
                                                   double indentDepth,
                                                   const std::vector<Membrane> &membranes,
                                                   std::vector<ExcZone> &excZones);

    void move(double dx, double dy);
};

#endif
