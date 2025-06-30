/*
 *  main.cpp
 *
 *  Main C++ file containing the body of the simulation.
 *
 *  James Bradford
 *  University of Sheffield
 *  March 2021
 */
#include "configHeader.h"
#include "Actin.h"
#include "print.h"
#include "fileoutput.h"
#include "globals.h"
#include "RNG.h"
//#include "pythoncaller.h"
//#include "createVid.h"
#include "ProteinRegion.h"
#include "geometry.h"
#include "nucleation.h"
#include "polymerisation.h"
#include "branching.h"
#include "capping.h"
#include "GactinGrid.h"
#include "ArpGrid.h"
#include "ExcZone.h"
#include "MembraneWall.h"
#include "Membrane.h"
#include "bDynamics.h"
#include "phagoAnalysis.h"
#include "StericGrid.h"
#include "severing.h"
#include "Cortex.h"
#include "crosslinking.h"


// -----------------------------------------------------------------------------
// For boost program options config file input

namespace po = boost::program_options;
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
        copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
        return os;
}
// -----------------------------------------------------------------------------
RNG rng; // random number generator is global

// Below for steric hindrance checks if needed
typedef gte::DCPQuery<double, gte::Segment<2, double>, gte::Segment<2, double> > RobustQuery;

int main(int argc, char *argv[])
{
        // VARIABLES

        // setup of sim and output files
        std::string config_file;
        // setup
        int nActin; // number of actin filaments
        std::uint32_t rn_seed; // rng seed
        int nActinLimit; // Maximum number of filaments
        bool nActinLimitBool = false; // has the sim reached the limit?
        double runTime; // Time (sec) at which sim ends
        bool printToFile; // can choose not to output a results file
        bool printLen; // can choose to print out total length of F-actin each frame
        bool copyCFG; // can choose to copy the config file used and output it
        bool printLogFile; // can choose to print log file containing details of sim
        bool printMonoTime; // can choose to print out positions of all monomers and their birthtimes
        int nFrames; // number of "frames" to print to file

        // Below: whether or not the sim is a "test", test sims will dump output
        // files to a "test" directory which then can be deleted later
        bool test;
        bool phagoAnalysis; // choose to output a file with measurements of phagocytosis


        // actin parameters
        bool steric; // choose to observe steric hindrance
        double segLength; // subunit length of actin filaments

        // kinetics
        bool poly_b; // barbed end polymerisation on/off
        bool poly_p; // pointed end polymerisation on/off
        bool depoly_b; // barbed end depolymerisation on/off
        bool depoly_p; // pointed end depolymerisation on/off
        double gActinConc; // (Initial) concentration of g-actin monomers
        double k_on_b; // barbed end polymerisation rate (per uM per sec)
        double k_off_b; // barbed end depolymerisation rate (per sec)
        double k_on_p; // pointed end polymerisation rate (per uM per sec)
        double k_off_p; // pointed end depolymerisation rate (per sec)

        // Below: dissociation which is if a trimer undergoes depolymerisation
        // it will be deleted, if dissociate is false then trimers
        // cannot depolymerise
        bool dissociate;

        // branching
        bool branching;
        double arp23Conc;
        bool branch_detach; // Choose to have "biochemical" branch detachment
        double k_db;
        int branch_space; // spacing in monomers between branches on filament
        // Below 2: When considering Brownian dynamics branch points are treated
        // as springs connecting the branch filament to its mother
        double branchStiffness; // spring constant (N/m) of branch point
        double torsionCoeff; // torsion spring constant (Nm/rad) of branch point

        // nucleation
        bool nucleation;
        // Nucleation rate is a rate per area per g-actin cubed
        // per second
        double k_nuc; // nucleation rate (per um**2 per uM**3 per sec)
        // Can choose to have directionality of nucleated filaments which involves
        // choosing angle of direction from a Gaussian distribution using RNG
        bool nucdirection;
        double nucAngleMean;
        double nucAngleStDev;

        // capping
        double capConc; // Concentration of capping protein (uM)

        bool capping_b; // Barbed end capping
        double k_cap_b; // Barbed end capping rate (per uM)

        bool capping_p; // Pointed end capping
        double k_cap_p; // Pointed end capping rate (per uM)

        bool uncapping_b; // Barbed end uncapping
        double k_uncap_b; // Barbed end uncapping rate

        bool uncapping_p; // Pointed end uncapping
        double k_uncap_p; // Pointed end uncapping rate

        // severing
        bool severing;
        double k_sever; // severing rate (per um (length of filament) per second)

        // Crosslinking
        bool crossLinking;
        double k_cLink; // crosslinking rate per potential link per second
        int cLSpacing; // spacing in monomers between crosslinks on filament
        double cLStiff; // stiffness of crosslink springs
        double cLDist; // distance threshold of crosslinks and rest length of springs

        // Unlinking
        bool unLinking;
        double k_unLink;

        // Brownian motion
        bool bDynamics; // Filaments use Brownian Dynamics algorithm instead
        double viscosity; // (Pa s)
        double temperature; // Kelvin

        // Tethering
        // Allows for the appearance of tether points along the filament
        // modelling what is seen in the TIRFM system
        // tether points are modelled as springs
        bool tether;
        double tetherDistMean; // mean distance between tether points along fila
        double tetherDistStdDev;
        double tetherStiffness; // stiffness of tether spring (N/m)

        // Option to have starting filaments with a different length
        bool initLen;
        int initLenMonos; // (initial length in monomers)
        bool initBent; // option to have filament bent if long enough

        // regions
        // These empty vectors MUST exist, regardless of if they get filled or not

        std::vector<std::string> nucRegionInput;
        std::vector<ProteinRegion> nucleationRegions;
        std::vector<std::string> nucRegionRInput;

        std::vector<std::string> branchRegionInput;
        std::vector<ProteinRegion> branchingRegions;

        std::vector<std::string> capRegionInput;
        std::vector<ProteinRegion> cappingRegions;

        std::vector<std::string> antiCapRegionInput;
        std::vector<ProteinRegion> antiCapRegions;

        std::vector<std::string> sevRegionInput;
        std::vector<ProteinRegion> sevRegions;

        // -------------------- Membranes ---------------------------------
        std::vector<std::string> memWallInput;
        std::vector<MembraneWall> memWalls;

        std::vector<std::string> memInput;
        Membrane membrane = Membrane();
        std::vector<Membrane> membranes;
        membranes.reserve(10);
        unsigned int timesFused = 0; // count of fusion events
        bool fusionEnd = false; // True when ending sim early due to fusion
        // set below to true if do not want fusion
        bool COMadjust; // If true then after the membrane moves everything is moved such that the Centre of (mass) membrane is unchanged
        bool memSprings;
        bool memExocytosis;
        double k_exo;
        double exo_mean;
        double exo_stdev;
        bool activeBranch;
        double activeBranchHeight;
        bool activeNuc;
        double activeNucHeight;
        bool activeCap;
        double activeCapHeight;
        bool activeAntiCap;
        double activeAntiCapHeight;
        bool activeSever;
        double activeSeverHeight;
        // -------------------- Cortex ---------------------------------
        bool cortexBool;
        Cortex cortex = Cortex();
        // -------------------- Zones of exclusion ---------------------------------
        std::vector<std::string> circInput;
        std::vector<std::string> rectInput;
        std::vector<std::string> triInput;

        double directMoveTime;
        double directMoveVelX, directMoveVelY;
        double indentDepth;
        bool firstTouch = false;
        // Theis empty vector MUST exist, regardless of if it gets filled or not
        std::vector<ExcZone> excZones { };
        // ----------------------- Steric Grid -------------------------------------
        std::vector<std::string> stericGridInput;
        StericGrid stericGrid = StericGrid();
        bool stericGridBool = false;
        double stericGridMinXY = 0;
        double stericGridMaxXY = 0;
        // ------------------ Monomeric actin grid (for limited pool) --------------
        std::vector<std::string> gGridInput;
        GactinGrid gActinGrid; // declare the grid here
        bool latA;
        // ----------------- Arp2/3 grid (for limited pool) ------------------------
        std::vector<std::string> arpGridInput;
        ArpGrid arpGrid; // declare the grid here
        bool arpPool = false;


        // Plotting options
        bool pythonPlot;
        double pythonTPF;
        std::string pythonFileType;
        std::string pythonFrameLimits;

        bool makeVid;
        int fps;

        std::random_device seed_generator;
        // Declare the supported options.
        po::options_description generic("Generic options");
        generic.add_options()
                ("help,h", "Produces this help message and exits")
        ;

        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        po::options_description hidden("Hidden options");
        hidden.add_options()
                ("config-file", po::value<std::string>(&config_file)->default_value("configExamples/default.cfg"), "Name of a configuration file")
        ;

        // Config file options
        po::options_description config("Configuration file options (=default_value)");
        config.add_options()

                ("seed", po::value<std::uint32_t>(&rn_seed)->default_value(seed_generator()," generated from std::random_device (entropy pool)"), "Random number seed" )
                ("nActin",po::value<int>(&nActin)->default_value(10), "Number of starting trimers")
                ("nActinLimit",po::value<int>(&nActinLimit)->default_value(500), "Maximum number of actin filaments")
                ("runtime",po::value<double>(&runTime)->default_value(10), "Runtime of simulation (seconds)")

                ("test",po::value<bool>(&test)->default_value(0, "off"), "Choose to output results to test directory")
                ("printRes",po::value<bool>(&printToFile)->default_value(1, "on"), "Choose to print results to file")
                ("printLen",po::value<bool>(&printLen)->default_value(0, "off"), "Choose to total lengths to file")
                ("nFrames",po::value<int>(&nFrames)->default_value(1001), "Number of frames to print out to results file, it will always print out the initial frame.")
                ("copyConfig", po::value<bool>(&copyCFG)->default_value(1,"on"),"Choose to copy the input config file to the output directory")
                ("printLog", po::value<bool>(&printLogFile)->default_value(1,"on"), "Choose to print the log file (strongly recommended)")
                ("printPhago", po::value<bool>(&phagoAnalysis)->default_value(0,"off"), "Choose to print measurements of engulfment during phagocytosis")
                ("printMonoTime", po::value<bool>(&printMonoTime)->default_value(0,"off"), "Choose to print positions and creation times of all monomers")

                ("steric", po::value<bool>(&steric)->default_value(true), "Steric hindrance bool")
                ("stericGrid", po::value< std::vector<std::string> >(&stericGridInput), "Steric Grid (minimum x,y pos (m), maximum x,y pos(m))")

                ("poly_b", po::value<bool>(&poly_b)->default_value(true), "Barbed end polymerisation bool")
                ("poly_p", po::value<bool>(&poly_p)->default_value(true), "Pointed end polymerisation bool")
                ("depoly_b", po::value<bool>(&depoly_b)->default_value(true), "Barbed end depolymerisation bool")
                ("depoly_p", po::value<bool>(&depoly_p)->default_value(true), "Pointed end depolymerisation bool")

                ("dissociate", po::value<bool>(&dissociate)->default_value(true), "Dissociation of trimers bool")

                ("gActinConc",po::value<double>(&gActinConc)->default_value(1), "Concentration of G-actin (uM)")
                ("k_poly_b",po::value<double>(&k_on_b)->default_value(12), "Barbed end polymerisation rate (monomers per uM G-actin per second)")
                ("k_poly_p",po::value<double>(&k_on_p)->default_value(1.3), "Pointed end polymerisation rate (monomers per uM G-actin per second)")
                ("k_depoly_b",po::value<double>(&k_off_b)->default_value(1.4, "1.4"), "Barbed end depolymerisation rate (monomers per second)")
                ("k_depoly_p",po::value<double>(&k_off_p)->default_value(0.8, "0.8"), "Pointed end depolymerisation rate (monomers per second)")

                ("branching", po::value<bool>(&branching)->default_value(0, "off"), "Branching bool")
                ("arpConc", po::value<double>(&arp23Conc)->default_value(0), "Concentration of Arp2/3 complex (uM)")
                ("branchDetach", po::value<bool>(&branch_detach)->default_value(0), "Branch detachment bool")
                ("k_db", po::value<double>(&k_db)->default_value(0.0018), "Branch detachment rate (/s)")
                ("branchSpacing", po::value<int>(&branch_space)->default_value(13), "Min allowed distance between branches (num of monomers)")
                ("branchStiff", po::value<double>(&branchStiffness)->default_value(5E-5,"5E-5"), "Stiffness of the branch springs (N/m)")
                ("torsionCoeff", po::value<double>(&torsionCoeff)->default_value(5E-18, "5E-18"), "Angular stiffness of the branch springs (Nm/rad)")

                ("nucleation", po::value<bool>(&nucleation)->default_value(0, "off"), "Nucleation bool")
                ("k_nuc", po::value<double>(&k_nuc)->default_value(0), "Nucleation rate density (per um^3 per second)")
                ("nucDir",po::value<bool>(&nucdirection)->default_value(0, "off"), "Nucleation direction bool")
                ("nucDirMean", po::value<double>(&nucAngleMean), "Gaussian mean of nucleation direction (factor of pi, (rad))")
                ("nucDirStDev", po::value<double>(&nucAngleStDev), "Gaussian standard deviation of nucleation direction (factor of pi, (rad))")

                ("capping_barb", po::value<bool>(&capping_b)->default_value(0, "off"), "Barbed end capping bool")
                ("capping_point", po::value<bool>(&capping_p)->default_value(0, "off"), "Pointed end capping bool")
                ("uncapping_barb", po::value<bool>(&uncapping_b)->default_value(0, "off"), "Barbed end uncapping bool")
                ("uncapping_point", po::value<bool>(&uncapping_p)->default_value(0, "off"), "Pointed end uncapping bool")
                ("cpConc", po::value<double>(&capConc)->default_value(0), "Concentration of capping protein (uM)")
                ("k_barb_cap",po::value<double>(&k_cap_b)->default_value(3.5), "Barbed end capping rate (per uM CP per second)")
                ("k_point_cap",po::value<double>(&k_cap_p)->default_value(0.8, "0.8"), "Pointed end capping rate (per uM CP per second)")
                ("k_barb_uncap",po::value<double>(&k_uncap_b)->default_value(4.0E-4, "4E-4"), "Barbed end uncapping rate (per second)")
                ("k_point_uncap",po::value<double>(&k_uncap_p)->default_value(1.8E-3, "1.8E-3"), "Pointed end uncapping rate (per second)")

                ("severing", po::value<bool>(&severing)->default_value(0, "off"), "Severing bool")
                ("k_sever", po::value<double>(&k_sever)->default_value(5E-3, "5E-3"), "Severing rate (per um Filament per second)")

                ("crosslinking", po::value<bool>(&crossLinking)->default_value(0, "off"), "Crosslinking bool")
                ("k_cLink", po::value<double>(&k_cLink)->default_value(1E-3, "1E-3"), "Crosslinking rate (per potential link per second)")
                ("cLSpacing", po::value<int>(&cLSpacing)->default_value(13), "Min allowed distance between crosslinks (num of monomers)")
                ("cLDist", po::value<double>(&cLDist)->default_value(1E-8, "1E-8"), "Crosslinking threshold distance (m)")
                ("cLStiff", po::value<double>(&cLStiff)->default_value(5E-5, "5E-5"), "Stiffness of crosslink springs (N/m)")

                ("unlinking", po::value<bool>(&unLinking)->default_value(0, "off"), "Unlinking of crosslinks bool")
                ("k_unLink", po::value<double>(&k_unLink)->default_value(1E-3, "1E-3"), "Unlinking rate (per crosslink per second)")

                ("bDynamics", po::value<bool>(&bDynamics)->default_value(0, "off"), "Filaments move using Brownian Dynamics algorithm bool")
                ("segLength", po::value<double>(&segLength)->default_value(1.0), "Desired segmentation length of a single filament (m)")

                ("tethering", po::value<bool>(&tether)->default_value(0, "off"), "Tethering bool")
                ("tetherDistMean", po::value<double>(&tetherDistMean)->default_value(1E-6, "1E-6"), "Mean distance between tethering points (m)")
                ("tetherDistStdev", po::value<double>(&tetherDistStdDev)->default_value(2.5E-7, "2.5E-7"), "Standard dev of distance between tethering points (m)")
                ("tetherStiff", po::value<double>(&tetherStiffness)->default_value(1E-5, "1E-5"), "Stiffness of the tether springs (N/m)")

                ("initLen", po::value<bool>(&initLen)->default_value(0, "off"), "Bool for initialising to a different length")
                ("initBent", po::value<bool>(&initBent)->default_value(0, "off"), "Bool for initialising to a bent configuration")
                ("initLenLength", po::value<int>(&initLenMonos)->default_value(3), "Length (in number of monomers) to initialise to")

                ("viscosity", po::value<double>(&viscosity)->default_value(0.02, "0.02"), "Viscosity of medium (Pascal seconds)")
                ("temp", po::value<double>(&temperature)->default_value(300), "Absolute temperature (Kelvin)")

                ("nucRegion", po::value< std::vector<std::string> >(&nucRegionInput), "Nucleation region (angle (factor of pi radians), bottom left x coord, bottom left y coord, width, height (metres), opt: nucleation coefficient, coupled to membrane bool, membrane subunit to couple to)")
                ("nucRegionR", po::value< std::vector<std::string> >(&nucRegionRInput), "Ringed Nucleation region (angle (factor of pi radians), bottom left x coord, bottom left y coord, width, height (metres), opt: nucleation coefficient, coupled to membrane bool, membrane subunit to couple to)")
                ("branchRegion", po::value< std::vector<std::string> >(&branchRegionInput), "Branching region (angle (factor of pi radians), bottom left x coord, bottom left y coord, width, height (metres)) opt: Arp2/3 conc, coupled to membrane bool, membrane subunit to couple to)")
                ("capRegion", po::value< std::vector<std::string> >(&capRegionInput), "Capping region (angle (factor of pi radians), bottom left x coord, bottom left y coord, width, height (metres), opt: capping coefficient, coupled to membrane bool, membrane subunit to couple to)")
                ("antiCapRegion", po::value< std::vector<std::string> >(&antiCapRegionInput), "Anti-Capping region (angle (factor of pi radians), bottom left x coord, bottom left y coord, width, height (metres), opt: coupled to membrane bool, membrane subunit to couple to)")
                ("sevRegion", po::value< std::vector<std::string> >(&sevRegionInput), "Severing region (angle (factor of pi radians), bottom left x coord, bottom left y coord, width, height (metres), opt: coupled to membrane bool, membrane subunit to couple to)")

                ("circle", po::value< std::vector<std::string> >(&circInput), "Exclusion zone; Circle: (centre in x, centre in y, radius (metres)) ")
                ("triangle", po::value< std::vector<std::string> >(&triInput), "Exclusion zone; Triangle: (xmin, ymin, xmax, ymax (metres))")
                ("rectangle", po::value< std::vector<std::string> >(&rectInput), "Exclusion zone; Rectangle: (xmin, ymin, xmax, ymax (metres))")
                ("directMoveTime", po::value<double>(&directMoveTime)->default_value(1000), "Time at which directed movement of circles and/or triangles start")
                ("directMoveVelX", po::value<double>(&directMoveVelX)->default_value(0, "0"), "Velocity in x direction (m/s)")
                ("directMoveVelY", po::value<double>(&directMoveVelY)->default_value(-4E-7, "-4E-7"), "Velocity in y direction (m/s)")
                ("indentDepth", po::value<double>(&indentDepth)->default_value(2E-7, "2E-7"), "How far to move before stopping")


                ("gActinGrid", po::value< std::vector<std::string> >(&gGridInput), "G actin grid (xmin, ymin, xmax, ymax, cellwidth, cellheight, celldepth (metres), opt:flow (bool))")
                ("arpGrid", po::value< std::vector<std::string> >(&arpGridInput), "Arp2/3 grid (xmin, ymin, xmax, ymax, cellwidth, cellheight, celldepth (metres), opt:flow (bool))")
                ("latA", po::value<bool> (&latA)->default_value(0, "off"), "Put latrunculin in the system")

                ("membraneWall", po::value< std::vector<std::string> >(&memWallInput), "Membrane wall (y-position, x-length (metres), external Force (Newtons))")

                ("membrane", po::value< std::vector<std::string> >(&memInput), "Membrane (y-position, x-length (metres), membrane bending mod (units of Kb*T), opt: membrane fusion bool, opt:membrane activation bool, opt: membrane activation distance threshold, opt: init regions replicating cortex), opt: membrane tension")
                ("memSprings", po::value<bool> (&memSprings)->default_value(0, "off"), "bool to have membrane with springs rather than constraints")
                ("memExocytosis", po::value<bool> (&memExocytosis)->default_value(0, "off"), "bool to have membrane grow with exocytosis")
                ("k_exo", po::value<double> (&k_exo)->default_value(0.1, "0.1"), "Exocytosis rate (per second)")
                ("exo_mean", po::value<double> (&exo_mean)->default_value(220E-9, "220E-9"), "Exocytosis length mean (metres)")
                ("exo_stdev", po::value<double> (&exo_stdev)->default_value(45E-9, "45E-9"), "Exocytosis length standard deviation (metres)")

                ("activeBranch", po::value<bool> (&activeBranch)->default_value(0, "off"), "bool to have activated branch regions under membrane")
                ("activeBranchHeight", po::value<double> (&activeBranchHeight)->default_value(200E-9, "200E-9"), "Height of activated branch region (m)")
                ("activeNuc", po::value<bool> (&activeNuc)->default_value(0, "off"), "bool to have activated nucleation regions under membrane")
                ("activeNucHeight", po::value<double> (&activeNucHeight)->default_value(200E-9, "200E-9"), "Height of activated nucleation region (m)")
                ("activeCap", po::value<bool> (&activeCap)->default_value(0, "off"), "bool to have activated capping regions under membrane")
                ("activeCapHeight", po::value<double> (&activeCapHeight)->default_value(200E-9, "200E-9"), "Height of activated capping region (m)")
                ("activeAntiCap", po::value<bool> (&activeAntiCap)->default_value(0, "off"), "bool to have activated antiCap regions under membrane")
                ("activeAntiCapHeight", po::value<double> (&activeAntiCapHeight)->default_value(200E-9, "200E-9"), "Height of activated antiCap region (m)")
                ("activeSever", po::value<bool> (&activeSever)->default_value(0, "off"), "bool to have activated sever regions under membrane")
                ("activeSeverHeight", po::value<double> (&activeSeverHeight)->default_value(200E-9, "200E-9"), "Height of activated sever region (m)")

                ("memRefFrame", po::value<bool> (&COMadjust)->default_value(0, "off"), "bool to set reference frame to be centre of membrane")
                ("cortex", po::value<bool> (&cortexBool)->default_value(0, "off"), "Cortex bool")

                ("pythonPlot", po::value<bool> (&pythonPlot)->default_value(1, "on"), "Python plotting bool")
                ("pythonTimeperFrame", po::value<double> (&pythonTPF)->default_value(0.2, "0.2"), "Length of time between frames in seconds")
                ("pythonFileType", po::value<std::string> (&pythonFileType)->default_value("png"), "File type for python plots, please only give png, eps or svg")
                ("pythonFrameLimits", po::value<std::string> (&pythonFrameLimits)->default_value("-5E-6, -5E-6, 5E-6, 5E-6"), "Frame ranges (xmin,ymin,xmax,ymax) (metres)")

                ("video", po::value<bool> (&makeVid)->default_value(0, "off"), "Automatically make video")
                ("videoFPS", po::value<int> (&fps)->default_value(20, "20"), "Frames per second of video")


        ;

        po::options_description cmdline_options;
        cmdline_options.add(generic).add(config).add(hidden);

        po::options_description visible("Allowed options");
        visible.add(generic).add(config);

        po::positional_options_description p;
        p.add("config-file", -1);


        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help"))
        {
                std::cout << "Usage: main [config_file.cfg]" << "\n";
                std::cout << visible << "\n";
                return 0;
        }

        std::ifstream ifs;
        ifs.open(config_file);
        if (!ifs.is_open())
        {
                std::cout << "can not open config file: " << config_file << "\n";
                return 0;
        }
        else
        {
                store(parse_config_file(ifs, config), vm);
                notify(vm);
        }

        if (vm.count("nucRegion"))
        {
                for (unsigned int i = 0; i < nucRegionInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(nucRegionInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 3)
                        {
                                // Circular region
                                // Arguments 0 and 1 are the x,y coords of centre
                                // Argument 2 is the radius
                                std::array<double,2> pointC = { vect[0], vect[1] };
                                nucleationRegions.push_back(ProteinRegion(vect[2], pointC, 0));

                        }
                        else if (vect.size() == 4)
                        {
                                // Circ region with coefficient
                                // Argument 3 is the nucleation coefficient
                                // So if you want nucleation in this region to be
                                // double the "base rate" this should be 2.0
                                std::array<double,2> pointC = { vect[0], vect[1] };
                                nucleationRegions.push_back(ProteinRegion(vect[2], pointC, 0, 0, vect[3]));
                        }
                        else if (vect.size() == 5)
                        {
                                //Rectangular region
                                // Argument 0 is the orientation angle in units of pi radians
                                // Arguments 1,2 is the coordinate of the bottom left point of the rectangle
                                // Argument 3 is the width
                                // Argument 4 is the height
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                nucleationRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, 0)));
                        }
                        else if (vect.size() == 6)
                        {
                                //Rectangular region with coefficient
                                // Argument 5 is the nucleation coefficient
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                nucleationRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, 0, 0, vect[5])));

                        }
                        else if (vect.size() == 8)
                        {
                                //Rectangular region with coefficient and stuck to a membrane
                                // Argument 6 is the bool for sticking to the membrane 1-yes, 0-no
                                // Argument 7 is the subunit of the membrane it should be coupled to
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                nucleationRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, 0, 0, vect[5], 0, vect[6], vect[7])));
                                if (vect[6] == 1)
                                {
                                        std::cout << "Coupling nucleation region to a membrane, subunit: " << vect[7] << std::endl;
                                }
                        }
                        else
                        {
                                std::cout << "Nucleation region must only have 5, 6 or 8 arguments" << std::endl;
                                return 0;
                        }
                }
        }

        if (vm.count("nucRegionR"))
        {
                for (unsigned int i = 0; i < nucRegionRInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(nucRegionRInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 4)
                        {
                                // Nucleation Ring region
                                // Argument 0,1 is the centre point
                                // Argument 2 is the outer radius of the ring
                                // Argument 3 is the inner radius of the ring
                                std::array<double,2> pointC = { vect[0], vect[1] };
                                nucleationRegions.push_back(ProteinRegion(pointC, vect[2], vect[3]));
                        }
                        else if (vect.size() == 5)
                        {
                            // Nucleation Ring region with nucleation coefficient
                            // Argument 4 is the nucleation coefficient
                            std::array<double,2> pointC = { vect[0], vect[1] };
                            nucleationRegions.push_back(ProteinRegion(pointC, vect[2], vect[3], vect[4]));
                        }
                        else if (vect.size() == 7)
                        {
                                // Nucleation Ring region with nucleation coefficient and coupling to a membrane sub
                                // Argument 5 is the bool for sticking to the membrane 1-yes, 0-no
                                // Argument 6 is the subunit of the membrane it should be coupled to
                                std::array<double,2> pointC = { vect[0], vect[1] };
                                nucleationRegions.push_back(ProteinRegion(pointC, vect[2], vect[3], vect[4], vect[5], vect[6]));
                        }
                        else
                        {
                                std::cout << "Ring Nucleation region must only have 4, 5 or 7 arguments" << std::endl;
                                return 0;
                        }
                }
        }

        if (vm.count("branchRegion"))
        {
                for (unsigned int i = 0; i < branchRegionInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(branchRegionInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }


                        if (vect.size() == 3)
                        {
                                // Circular region
                                std::array<double,2> pointC = { vect[0], vect[1] };
                                branchingRegions.push_back(ProteinRegion(vect[2], pointC, arp23Conc));
                        }
                        else if (vect.size() == 4)
                        {
                                // Circular region with local arp2/3 conc
                                // Argument 3 is a CONCENTRATION not a coefficient as with other regions
                                std::array<double,2> pointC = { vect[0], vect[1] };
                                branchingRegions.push_back(ProteinRegion(vect[2], pointC, vect[3]));
                        }
                        else if (vect.size() == 5)
                        {
                                //Rectangular region
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                branchingRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, arp23Conc)));
                        }
                        else if (vect.size() == 6)
                        {
                                //Rectangular region with local arp2/3 conc
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                branchingRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, vect[5])));
                        }
                        else if (vect.size() == 8)
                        {
                                //Rectangular region with local arp2/3 conc and coupling to a membrane
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                branchingRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, vect[5], 0, 0, 0, vect[6], vect[7])));
                                if (vect[6] == 1)
                                {
                                        std::cout << "Coupling branching region to a membrane, subunit: " << vect[7] << std::endl;
                                }
                        }
                        else
                        {
                                std::cout << "Branching region must only have 3, 4, 5, 6 or 8 arguments" << std::endl;
                                return 0;
                        }
                }
        }

        if (vm.count("capRegion"))
        {
                for (unsigned int i = 0; i < capRegionInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(capRegionInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 3)
                        {
                                // Circular region
                                std::array<double,2> pointC = { vect[0], vect[1] };
                                cappingRegions.push_back(ProteinRegion(vect[2], pointC, 0));
                        }
                        else if (vect.size() == 4)
                        {
                                // Circular region with coefficient
                                std::array<double,2> pointC = { vect[0], vect[1] };
                                cappingRegions.push_back(ProteinRegion(vect[2], pointC, 0, vect[3]));
                        }
                        else if (vect.size() == 5)
                        {
                                //Rectangular region
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                cappingRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, 0)));
                        }
                        else if (vect.size() == 6)
                        {
                                //Rectangular region with coefficient
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                cappingRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, 0, vect[5])));
                        }
                        else if (vect.size() == 8)
                        {
                                //Rectangular region withcoefficient and coupling to membrane
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                cappingRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, 0, vect[5], 0, 0, vect[6], vect[7])));
                                if (vect[6] == 1)
                                {
                                        std::cout << "Coupling capping region to a membrane, subunit: " << vect[7] << std::endl;
                                }
                        }
                        else
                        {
                                std::cout << "Capping region must only have 5, 6 or 8 arguments" << std::endl;
                                return 0;
                        }
                }
        }
        if (vm.count("antiCapRegion"))
        {
                for (unsigned int i = 0; i < antiCapRegionInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(antiCapRegionInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }
                        if (vect.size() == 3)
                        {
                                // Circular region
                                std::array<double,2> pointC = { vect[0], vect[1] };
                                antiCapRegions.push_back(ProteinRegion(vect[2], pointC, 0));
                        }

                        else if (vect.size() == 5)
                        {
                                // Rectangular region
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                antiCapRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, 0)));
                        }
                        else if (vect.size() == 7)
                        {
                                // Rectangular region allow for coupling to membrane
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                antiCapRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, 0, 0, 0, 0, vect[5], vect[6])));
                                if (vect[5] == 1)
                                {
                                        std::cout << "Coupling Anti-capping region to a membrane, subunit: " << vect[6] << std::endl;
                                }
                        }
                        else
                        {
                                std::cout << "Anti-Capping region must only have 3, 5 or 7 arguments" << std::endl;
                                return 0;
                        }
                }
        }

        if (vm.count("sevRegion"))
        {
                for (unsigned int i = 0; i < sevRegionInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(sevRegionInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 3)
                        {
                                // Circular region
                                std::array<double,2> pointC = { vect[0], vect[1] };
                                sevRegions.push_back(ProteinRegion(vect[2], pointC, 0));
                        }
                        else if (vect.size() == 4)
                        {
                                // Circular region with sever coefficient
                                std::array<double,2> pointC = { vect[0], vect[1] };
                                sevRegions.push_back(ProteinRegion(vect[2], pointC, 0, 0, 0, vect[3]));

                        }
                        else if (vect.size() == 5)
                        {
                                // Rectangular region
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                sevRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, 0)));
                        }
                        else if (vect.size() == 6)
                        {
                                // Rectangular region with sever coefficient
                                std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                                std::array<double,2> pointBL = { vect[1], vect[2] };
                                sevRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, 0, 0, 0, vect[5])));

                        }
                        else if (vect.size() == 8)
                        {
                            // Rectangular region with sever coefficient and allowing coupling to membrane
                            std::array<double,2> uVec = { cos(vect[0]*M_PI), sin(vect[0]*M_PI) };
                            std::array<double,2> pointBL = { vect[1], vect[2] };
                            sevRegions.push_back((ProteinRegion(uVec, vect[3], vect[4], pointBL, 0, 0, 0, vect[5], vect[6], vect[7])));

                            if (vect[6] == 1)
                            {
                                    std::cout << "Coupling severing region to a membrane" << std::endl;
                            }
                        }
                        else
                        {
                                std::cout << "Severing region must only have 3, 4, 5, 6 or 8 arguments" << std::endl;
                                return 0;
                        }
                }
        }

        if (vm.count("circle"))
        {
                for (unsigned int i = 0; i < circInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(circInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 3)
                        {
                                // 3 arguments given: Circle, no movement
                                // Argument 0,1 x,y coordinates of circle centre
                                // Argument 2: radius of circle
                                excZones.push_back(ExcZone(vect[0], vect[1], vect[2], excZones.size()));
                        }
                        else if (vect.size() == 4)
                        {
                                // 4 arguments given: Circle, with/without brownian motion
                                // Argument 3: Movement bool 1 - BM, 0 - no BM
                                excZones.push_back(ExcZone(vect[0], vect[1], vect[2], excZones.size(), vect[3]));
                        }
                        else if (vect.size() == 5)
                        {
                            // 5 arguments given: Circle, with/without brownian motion, with/without directed movement
                            // Argument 4: Directed movement bool 1- Yes, 0 - No
                            excZones.push_back(ExcZone(vect[0], vect[1], vect[2], excZones.size(), vect[3], vect[4]));
                        }
                        else
                        {
                            std::cout << "Circular exclusion zone must only have 3, 4 or 5 arguments" << std::endl;
                            return 0;
                        }
                }
        }

        if (vm.count("rectangle"))
        {
                for (unsigned int i = 0; i < rectInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(rectInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 4)
                        {
                                // 4 arguments given: Rectangular exclusion zone
                                // Argument 0: Min X coordinate (left side)
                                // Argument 1: Min Y coordinate (bottom side)
                                // Argument 2: Max x coordinate (right side)
                                // Argument 3: Max Y coordinate (top side)
                                excZones.push_back(ExcZone(vect[0], vect[1], vect[2], vect[3], excZones.size(), 0, 0));
                        }
                        else
                        {
                                std::cout << "Rectangular exclusion zone must only have 4 arguments" << std::endl;
                                return 0;
                        }
                }
        }

        if (vm.count("triangle"))
        {
                for (unsigned int i = 0; i < triInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(triInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 4)
                        {
                                // 4 arguments given: triangular exclusion zone
                                // Argument 0, 1: x,y coordinates of bottom point
                                // Argument 2: Height of triangle
                                // Argument 3: Top width of triangle
                                excZones.push_back(ExcZone(vect[0], vect[1], vect[2], vect[3], excZones.size(), 0, 0, 0));
                        }
                        else if (vect.size() == 5)
                        {
                            // 5 arguments given: triangle as above, plus with/without directed movement
                            // Argument 4: Directed movement bool
                            excZones.push_back(ExcZone(vect[0], vect[1], vect[2], vect[3], excZones.size(), vect[4], 0, 0));
                        }
                        else
                        {
                                std::cout << "Triangular exclusion zone must only have 4 or 5 arguments" << std::endl;
                                return 0;
                        }
                }
        }

        if (vm.count("gActinGrid"))
        {
                for (unsigned int i = 0; i < gGridInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(gGridInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 7)
                        {
                                // 7 arguments given: G-actin grid
                                // Argument 0: Min X coordinate (left side)
                                // Argument 1: Min Y coordinate (bottom side)
                                // Argument 2: Max x coordinate (right side)
                                // Argument 3: Max Y coordinate (top side)
                                // Argument 4: Width of one cell
                                // Argument 5: Height of one cell
                                // Argument 6: Depth of one cell (needed for calculating conc -> number)

                                gActinGrid = GactinGrid(vect[0], vect[1], vect[2], vect[3], vect[4], vect[5], vect[6], gActinConc);
                        }
                        else if (vect.size() == 8)
                        {
                                // 8 arguments given: G-actin grid with latrunculin
                                // Argument 7: latrunculin concentration
                                gActinGrid = GactinGrid(vect[0], vect[1], vect[2], vect[3], vect[4], vect[5], vect[6], gActinConc, vect[7]);
                        }
                        else if (vect.size() == 9)
                        {
                                // 9 arguments given: G-actin grid with latrunculin and flow between cells
                                // Argument 8: Flow bool 1 - Flow, 0- No Flow
                                gActinGrid = GactinGrid(vect[0], vect[1], vect[2], vect[3], vect[4], vect[5], vect[6], gActinConc, vect[7], vect[8]);
                        }
                        else
                        {
                                std::cout << "gActin grids must only have 7, 8 or 9 arguments" << std::endl;
                                return 0;
                        }
                }
        }

        if (vm.count("arpGrid"))
        {
                arpPool = true;
                for (unsigned int i = 0; i < arpGridInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(arpGridInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 7)
                        {
                                // 7 arguments given: Arp2/3 grid
                                // Argument 0: Min X coordinate (left side)
                                // Argument 1: Min Y coordinate (bottom side)
                                // Argument 2: Max x coordinate (right side)
                                // Argument 3: Max Y coordinate (top side)
                                // Argument 4: Width of one cell
                                // Argument 5: Height of one cell
                                // Argument 6: Depth of one cell (needed for calculating conc -> number)
                                arpGrid = ArpGrid(vect[0], vect[1], vect[2], vect[3], vect[4], vect[5], vect[6], arp23Conc);
                        }
                        else if (vect.size() == 8)
                        {
                                // 8 arguments given: Arp2/3 grid with flow
                                // Argument 7: Flow bool
                                arpGrid = ArpGrid(vect[0], vect[1], vect[2], vect[3], vect[4], vect[5], vect[6], arp23Conc, vect[7]);
                        }
                        else
                        {
                                std::cout << "arp2/3 grids must only have 7 or 8 arguments" << std::endl;
                                return 0;
                        }
                }
        }

        if (vm.count("membraneWall"))
        {
                for (unsigned int i = 0; i < memWallInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(memWallInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 3)
                        {
                                // 3 arguments given
                                // Argument 0: Y coordinate
                                // Argument 1: Length of wall (it is finite)
                                // Argument 2: External force pushing down (Newtons)
                                memWalls.push_back(MembraneWall(vect[0], vect[1], vect[2]));
                        }
                        else
                        {
                                std::cout << "MembraneWall must only have 3 arguments" << std::endl;
                                return 0;
                        }
                }
        }

        if (vm.count("membrane"))
        {
                for (unsigned int i = 0; i < memInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(memInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 3)
                        {
                                // 3 arguments given: Membrane
                                // Argument 0 y coordinate of centre of membrane
                                // Argument 1: Length of membrane
                                // Argument 2: Bending modulus of membrane (multiples of KT)
                                membrane = Membrane(vect[0], vect[1], temperature, vect[2]);
                        }
                        else if (vect.size() == 4)
                        {
                                // 4 Arguments given: Membrane with fusable bool
                                // Argument 3: Fusable bool 1 - Can fuse, 0 - Can't fuse
                                membrane = Membrane(vect[0], vect[1], temperature, vect[2], vect[3]);
                        }
                        else if (vect.size() == 5)
                        {
                                // 5 Arguments given: Membrane with fusable bool and active regions
                                // Argument 4: Bool allowing for active regions to turn on
                                // If within 10nm from target
                                membrane = Membrane(vect[0], vect[1], temperature, vect[2], vect[3], vect[4]);
                        }
                        else if (vect.size() == 6)
                        {
                            // 6 Arguments given: Membrane with fusable bool and active regions and active region distance
                            // Argument 5: Active region distance (if want to change from 10nm as above)
                            membrane = Membrane(vect[0], vect[1], temperature, vect[2], vect[3], vect[4], vect[5]);
                        }
                        else if (vect.size() == 7)
                        {
                            // 7 Arguments given: Membrane with fusable bool and active regions and active region distance and initialising with a cortex
                            // Argument 6: Initialise cortex bool
                            membrane = Membrane(vect[0], vect[1], temperature, vect[2], vect[3], vect[4], vect[5], vect[6]);
                        }
                        else if (vect.size() == 8)
                        {
                            // 8 Arguments given: Membrane with fusable bool and active regions and active region distance and initialising with a cortex and membrane tension
                            // Argument 6: Membrane tension
                            assert(memSprings);
                            membrane = Membrane(vect[0], vect[1], temperature, vect[2], vect[3], vect[4], vect[5], vect[6], vect[7]);

                        }
                        else
                        {
                                std::cout << "Membrane must only have 3, 4, 5, 6, 7 or 8 arguments" << std::endl;
                                return 0;
                        }
                }

        }
        // Push back to vector, even if there is no membrane!
        membranes.push_back(membrane);

        if (vm.count("stericGrid"))
        {
                for (unsigned int i = 0; i < stericGridInput.size(); ++i)
                {
                        // vect will store our arguments
                        std::vector<double> vect;
                        // take in a line
                        std::stringstream ss(stericGridInput[i]);
                        // j is our value
                        double j;

                        while (ss >> j)
                        {
                                vect.push_back(j);
                                // ignore white space and commas
                                if (ss.peek() == ',' || ss.peek() == ' ')
                                        ss.ignore();
                        }

                        if (vect.size() == 2)
                        {
                                // 2 arguments given: correct
                                // Argument 0: Min xy (left and bottom length)
                                // Argument 1: Max xy (top and right length)
                                stericGridBool = true;
                                stericGridMinXY = vect[0];
                                stericGridMaxXY = vect[1];
                        }
                        else
                        {
                                std::cout << "Steric Grid must only have 2 arguments" << std::endl;
                                return 0;
                        }
                }
        }
        // COUPLING OF REGIONS TO EITHER MEMBRANE WALL OR MEMBRANE
        for (unsigned int i = 0; i < nucleationRegions.size(); ++i)
        {
                if (nucleationRegions[i].getCoupledBool())
                {
                        if (memWalls.size() == 0 && !membranes[0].getExist())
                        {
                                std::cout << "There is no membrane to couple the nucleation region to!";
                                std::cout << " Please check the config file supplied" << std::endl;
                                exit(0);
                        }
                        else if (!membranes[0].getExist())
                        {
                                memWalls[0].appendToNucRegions(i);
                        }
                        else
                        {
                                nucleationRegions[i].stickToMem(membranes[0]);
                                membranes[0].addToNucReg(nucleationRegions[i].getCoupledSub(), i);
                        }
                }
        }
        for (unsigned int i = 0; i < branchingRegions.size(); ++i)
        {
                if (branchingRegions[i].getCoupledBool())
                {
                        if (memWalls.size() == 0 && !membranes[0].getExist())
                        {
                                std::cout << "There is no membrane to couple the branch region to!";
                                std::cout << " Please check the config file supplied" << std::endl;
                                exit(0);
                        }
                        else if (!membranes[0].getExist())
                        {
                                memWalls[0].appendToBranchRegions(i);
                        }
                        else
                        {
                                branchingRegions[i].stickToMem(membranes[0]);
                                membranes[0].addToBranchReg(branchingRegions[i].getCoupledSub(), i);
                        }
                }
        }
        for (unsigned int i = 0; i < cappingRegions.size(); ++i)
        {
                if (cappingRegions[i].getCoupledBool())
                {
                        if (memWalls.size() == 0 && !membranes[0].getExist())
                        {
                                std::cout << "There is no membrane to couple the capping region to!";
                                std::cout << " Please check the config file supplied" << std::endl;
                                exit(0);
                        }
                        else if (!membranes[0].getExist())
                        {
                                memWalls[0].appendToCapRegions(i);
                        }
                        else
                        {
                                cappingRegions[i].stickToMem(membranes[0]);
                                membranes[0].addToCapReg(cappingRegions[i].getCoupledSub(), i);
                        }
                }
        }
        for (unsigned int i = 0; i < antiCapRegions.size(); ++i)
        {
                if (antiCapRegions[i].getCoupledBool())
                {
                        if (memWalls.size() == 0 && !membranes[0].getExist())
                        {
                                std::cout << "There is no membrane to couple the capping region to!";
                                std::cout << " Please check the config file supplied" << std::endl;
                                exit(0);
                        }
                        else if (!membranes[0].getExist())
                        {
                                std::cout << "Anti cap regions not currently compatible with memWall" << std::endl;
                                //memWalls[0].appendToCapRegions(i);
                        }
                        else
                        {
                                antiCapRegions[i].stickToMem(membranes[0]);
                                membranes[0].addToAntiCapReg(antiCapRegions[i].getCoupledSub(), i);
                        }
                }
        }
        for (unsigned int i = 0; i < sevRegions.size(); ++i)
        {
                if (sevRegions[i].getCoupledBool())
                {
                        if (memWalls.size() == 0 && !membranes[0].getExist())
                        {
                                std::cout << "There is no membrane to couple the severing region to!";
                                std::cout << " Please check the config file supplied" << std::endl;
                                exit(0);
                        }
                        else if (!membranes[0].getExist())
                        {
                            memWalls[0].appendToSevRegions(i);
                        }
                        else
                        {
                            sevRegions[i].stickToMem(membranes[0]);
                            membranes[0].addToSevReg(sevRegions[i].getCoupledSub(), i);
                        }
                }
        }


        // Cortex option

        if (membranes[0].getExist() && cortexBool)
        {
            // Create a cortex under the membrane
            cortex = Cortex(temperature, membranes[0]);
        }
        else if (cortexBool)
        {
            std::cout << "Error: Cannot have a cortex without a membrane, check config file" << std::endl;
            exit(1);
        }

        // -------------------------- Timestep -------------------------------------
        double t_step = 1E-2; // max timestep is 1E-2
        double min_dt_bend, D_stokes, a;
        double Yama_corr_X = 0.044;

        double fracErr = 0.2;
        // Set the biochemical probabilitys to max of 0.1 per timestep
        if (k_on_b != 0 && poly_b)
        {
                t_step = std::min((0.1 / (k_on_b * gActinConc)), t_step);
        }

        if (k_on_p != 0 && poly_p)
        {
                t_step = std::min((0.1 / (k_on_p * gActinConc)), t_step);
        }

        if (k_off_b != 0 && depoly_b)
        {
                t_step = std::min((0.1 / k_off_b), t_step);
        }

        if (k_off_p != 0 && depoly_p)
        {
                t_step = std::min((0.1 / k_off_p), t_step);
        }

        std::cout << "Actin biochem tstep: " << t_step << std::endl;
        // Any springs
        if (branching && bDynamics)
        {
                min_dt_bend = fracErr*((2*M_PI*viscosity*segLength)/(log(segLength/(Actin::s_trueRadius*2))+Yama_corr_X) / branchStiffness);
                t_step = std::min(min_dt_bend, t_step);
                std::cout << "Actin branching tstep: " << min_dt_bend << std::endl;
        }

        if (tether)
        {
                min_dt_bend = fracErr*((2*M_PI*viscosity*segLength)/(log(segLength/(Actin::s_trueRadius*2))+Yama_corr_X) / tetherStiffness);
                t_step = std::min(min_dt_bend, t_step);
                std::cout << "Actin tether tstep: " << min_dt_bend << std::endl;

        }

        if (crossLinking && bDynamics)
        {
            min_dt_bend = fracErr*((2*M_PI*viscosity*segLength)/(log(segLength/(Actin::s_trueRadius*2))+Yama_corr_X) / cLStiff);
            t_step = std::min(min_dt_bend, t_step);
            std::cout << "Actin crosslinking tstep: " << min_dt_bend << std::endl;
        }


        if (bDynamics)
        {
                D_stokes = (g_Kb*temperature)*((log(segLength/(Actin::s_trueRadius*2))+Yama_corr_X)/(2*M_PI*segLength*viscosity));
                a = (g_Kb * temperature) / (6*M_PI*viscosity*D_stokes);
                min_dt_bend = ((fracErr*segLength)*(fracErr*segLength)) / (2*D_stokes);
                std::cout << "Actin bDynamics 1 tstep: " << min_dt_bend << std::endl;

                t_step = std::min(min_dt_bend, t_step);
                min_dt_bend = (fracErr*(M_PI*viscosity*a*segLength*segLength*segLength)/(Actin::s_persistenceLength*g_Kb*temperature));
                t_step = std::min(min_dt_bend, t_step);
                std::cout << "Actin bDynamics 2 tstep: " << min_dt_bend << std::endl;

        }

        if (membranes[0].getExist())
        {
                // Bead rod

                D_stokes = (g_Kb*temperature)*((log(membranes[0].getSegLength()/(membranes[0].getThickness()))+Yama_corr_X)/(2*M_PI*membranes[0].getSegLength()*viscosity));

                a = (g_Kb * temperature) / (6*M_PI*viscosity*D_stokes);
                min_dt_bend = ((fracErr*membranes[0].getSegLength())*(fracErr*membranes[0].getSegLength())) / (2*D_stokes);
                std::cout << "Membrane bDynamics 1 tstep: " << min_dt_bend << std::endl;

                t_step = std::min(min_dt_bend, t_step);
                min_dt_bend = (fracErr*(M_PI*viscosity*a*membranes[0].getSegLength()*membranes[0].getSegLength()*membranes[0].getSegLength())/(membranes[0].get1DBendMod()));
                t_step = std::min(min_dt_bend, t_step);
                std::cout << "Membrane bDynamics 2 tstep: " << min_dt_bend << std::endl;

        }


        if (membranes[0].getExist() && memSprings)
        {
                // Bead spring

                double min_dt_spr = fracErr*((2*M_PI*viscosity*membranes[0].getSegLength())/(log(membranes[0].getSegLength()/(membranes[0].getThickness()))+0.044) / (membranes[0].getStiffness()));
                t_step = std::min(min_dt_spr, t_step);
                std::cout << "Membrane spring tstep: " << min_dt_spr << std::endl;
        }

        if (cortex.getExist())
        {
            // Cortex
            D_stokes = (g_Kb*temperature)*((log(cortex.getSegLength()/(cortex.getThickness()))+Yama_corr_X)/(2*M_PI*cortex.getSegLength()*viscosity));

            a = (g_Kb * temperature) / (6*M_PI*viscosity*D_stokes);
            min_dt_bend = ((fracErr*cortex.getSegLength())*(fracErr*cortex.getSegLength())) / (2*D_stokes);
            std::cout << "Cortex bDynamics 1 tstep: " << min_dt_bend << std::endl;

            t_step = std::min(min_dt_bend, t_step);
            min_dt_bend = (fracErr*(M_PI*viscosity*a*cortex.getSegLength()*cortex.getSegLength()*cortex.getSegLength())/(cortex.get1DBendMod()));
            t_step = std::min(min_dt_bend, t_step);
            std::cout << "Cortex bDynamics 2 tstep: " << min_dt_bend << std::endl;

        }

        Actin::s_branchSpacing = (branching) ? branch_space : 0;
        Actin::s_segmentationLength = segLength;
        Actin::s_branchStiff = branchStiffness;
        Actin::s_tetherStiff = tetherStiffness;
        Actin::s_torsionCoeff = torsionCoeff;
        Actin::s_cLSpacing = (crossLinking) ? cLSpacing : 0;
        Actin::s_cLStiff = cLStiff;
        Actin::s_cLDist = cLDist;
        Actin::s_maxSpacing = (Actin::s_branchSpacing > Actin::s_cLSpacing ? Actin::s_branchSpacing : Actin::s_cLSpacing);

        std::array<double,2> AFMVelVector = { directMoveVelX, directMoveVelY };

        assert(t_step > 0);
        //t_step = 1E-4;
        std::cout << "Timestep: " << t_step << std::endl;

        // -------------------- Set RNG Seed -----------------------------------
        rng.setSeed(rn_seed);
        // ------------------------ Setup --------------------------------------

        clock_t walltime;
        double starttime = omp_get_wtime();
        // Set number of threads to use, if you change this you need to completely
        // rebuild, so "make clean" then "make"
        omp_set_num_threads(2);
        std::cout << "Number of OMP threads: " << omp_get_max_threads() << std::endl;

        // ---------------- Create results output file -----------------------------
        // This is the number of timesteps between printing to file, cannot be lower
        // than 1 (Which would print each timestep). (Always print the 0th frame)
        nFrames -= 1;
        int n_bw_print = std::max(int(runTime/(t_step*double(nFrames))), 1);
        double time_bw_print = n_bw_print*t_step;

        // Create a new directory if one doesn't exist
        std::string outdir = newDirectory(test);

        // Create a results file in this directory for writing to
        std::ofstream file_results;
        if (printToFile)
                file_results = openResults(outdir);

        std::ofstream file_lengths;
        if (printLen)
                file_lengths = openLengths(outdir);

        std::ofstream fileLog;
        if (printLogFile)
                fileLog = openLog(outdir);

        // Create a copy of the input configuration file, and dump it in the directory
        if (copyCFG)
                copyConfig(outdir, ifs);

        std::ofstream file_mono;
        if (printMonoTime)
            file_mono = openMono(outdir);

        std::ofstream file_concMap;
        if (gActinGrid.getExist())
                file_concMap = openConcMap(outdir);

        std::ofstream file_LatAconcMap;
        if (gActinGrid.getExist() && latA)
                file_LatAconcMap = openLatAConcMap(outdir);

        std::ofstream file_LatActinconcMap;
        if (gActinGrid.getExist() && latA)
                file_LatActinconcMap = openLatActinConcMap(outdir);

        std::ofstream file_ArpConcMap;
        if (arpPool)
                file_ArpConcMap = openArpConcMap(outdir);

        std::ofstream file_memLengths;
        std::ofstream file_fusionTime;
        if (membranes[0].getExist() && membranes[0].getFusable())
        {
            file_memLengths = openMemLengths(outdir);
            if (excZones.size() > 0)
            {
                file_fusionTime = openFusionTime(outdir);
            }
            else
            {
                std::cout << "Membrane is fusable but no targets to fuse!";
                std::cout << " Ending" << std::endl;
                exit(1);
            }
        }
        else
        {
            timesFused = excZones.size();
        }

        std::ofstream file_engulf;
        if (phagoAnalysis)
        {
                file_engulf = openEngulfment(outdir);
        }


        // --------------- Steric grid setup -----------------------------------

        if (steric && stericGridBool)
        {
                stericGrid = StericGrid(stericGridMinXY, stericGridMaxXY);
                if (membranes[0].getExist())
                {
                        for (int i = 0; i < membranes[0].getNumPoints(); ++i)
                        {
                                stericGrid.updateMembrane(membranes[0], i);
                        }
                }

                if (cortex.getExist())
                {
                        for (int i = 0; i < cortex.getNumPoints(); ++i)
                        {
                                stericGrid.updateCortex(cortex, i);
                        }
                }

                if (memWalls.size() != 0)
                {
                    stericGrid.updateMemWall(memWalls[0]);
                }

        }
        else if (steric)
        {
                std::cout << "You must define a steric grid! Or turn off steric";
                std::cout << " hindrance. Exiting..." << std::endl;
                exit(1);
        }


        // -------------- Init regions from membrane ---------------------------
        // Create any "active" regions based on initcondition
        if (membranes[0].getExist() && membranes[0].getActive())
        {
                membranes[0].initActRegion(branchingRegions, cappingRegions,
                                           antiCapRegions,
                                           nucleationRegions, sevRegions,
                                           excZones, activeBranch,
                                           activeBranchHeight, arp23Conc,
                                           activeNuc,
                                           activeNucHeight, activeCap,
                                           activeCapHeight, activeAntiCap,
                                           activeAntiCapHeight, activeSever,
                                           activeSeverHeight);
        }

        if (membranes[0].getExist() && membranes[0].getInitCortex())
        {
            membranes[0].initCortexRegions(branchingRegions, nucleationRegions, arp23Conc);
        }

        //----------------------------------------------------------------------
        // Number of starter trimers
        std::vector<Actin> actinVec {  };

        // Generate distribution ranging from 0 to nRegions-1 (inclusive)
        // TO DO: WRITE FUNCTIONS FOR THESE
        double totalRegionAreaNuc { 0.0 };
        for (unsigned int i = 0; i < nucleationRegions.size(); ++i)
        {
                totalRegionAreaNuc += nucleationRegions[i].getArea();
        }

        // ----- free nucleation angle distribution --------------------------------

        //default is a uniform distribution from 0 to 2pi but can change this
        nucAngleMean *= M_PI;
        nucAngleStDev *= M_PI;

        // Split up nucleation regions if they are too big
        for (unsigned int i = 0; i < nucleationRegions.size(); ++i)
        {
                double nucProb = k_nuc * t_step * gActinConc * gActinConc * gActinConc * nucleationRegions[i].getNucCoeff() * (nucleationRegions[i].getArea() * 1E12);
                while (nucProb > 0.1)
                {
                        nucleationRegions[i].splitRegion(nucleationRegions,
                                                         memWalls);
                        nucProb /= 2;
                }
        }

        // Used for initial nucleation
        for (unsigned int i = 0; i < nucleationRegions.size(); ++i)
        {
                nucleationRegions[i].setNormalisedArea(totalRegionAreaNuc);
        }

        // ----------- Tethering distribution --------------------------------------
        auto tetherDistri = std::normal_distribution<double> (tetherDistMean, tetherDistStdDev);

        // ---------------- Nucleate starter trimers at time 0 ---------------------

        initNuc(actinVec, nActin, nucleationRegions, excZones, membranes,
                memWalls, cortex, steric, stericGrid, tether,
                tetherDistri, nucdirection, nucAngleMean, nucAngleStDev);

        // -------------------- initialisation -------------------------------------

        if (initLen)
        {

                for (unsigned int i = 0; i < actinVec.size(); ++i)
                {
                        actinVec[i].initialiseLength(initLenMonos, actinVec,
                                                     steric, stericGrid,
                                                     membranes, cortex,
                                                     excZones, memWalls, tether);
                }

                if (initBent)
                {
                        initConfiguration(actinVec, excZones, memWalls,
                                          membranes, cortex, steric,
                                          stericGrid);
                }

        }

        // reseed RNG here if needed
        std::uint32_t origSeed = rng.getSeed();
        //rng.setSeed(seed_generator()); // reseed

        // --------------------- Simulation conditions--------------------------

        double currentTime { 0.0 };

        // -------------------- Branching --------------------------------------
        // Branching rate per monomer on the mother
        // We now use Carlsson 2004 model 5.3x10^-4 (end plus side model)
        // This is now done on a per monomer level, not per length
        // Carlsson (2004) models
        // HardCode in the rate, this is dependent on C_arp*C_gactin*C_gactin
        const double k_branch { 5.3E-4 };

        // Dislocation of branches
        // Carlsson (0.0018/s) Pollard (0.0013/s) use greater of the two
        // If Cofilin present increases to 0.023/s (Pollard 2000)
        const double p_detach { k_db * t_step };

        // ------------------- Set maximum number of filaments --------------
        actinVec.reserve(nActinLimit);
        // ------------------- Capping and Uncapping -------------------------------
        // Capping for barbed and/or pointed ends
        double p_cap_b { k_cap_b * capConc * t_step};

        double p_cap_p { k_cap_p * capConc * t_step };

        double p_uncap_b { k_uncap_b * t_step };

        double p_uncap_p { k_uncap_p * t_step };

        // --------------- Print out parameters and initial frame ------------------

        if (printLogFile)
        {
                printLog(fileLog, t_step, origSeed, fracErr);
        }
        if (printToFile)
        {
                printActinHeaders(file_results, time_bw_print);
                printActinframes(file_results, nActin, actinVec, 0.0, t_step, excZones,
                                 nucleationRegions, branchingRegions, cappingRegions,
                                 antiCapRegions, sevRegions,
                                 memWalls, membranes, tether, crossLinking,
                                 cortex);
        }
        if (printMonoTime)
        {
            printActinHeaders(file_mono, time_bw_print);
            printMonoTimeFrames(file_mono, nActin, actinVec, 0.0, t_step);
        }

        if (gActinGrid.getExist())
        {
                printGActinGridHeader(file_concMap, time_bw_print, gActinGrid);
                if (latA)
                {
                  printGActinGridHeader(file_LatAconcMap, time_bw_print, gActinGrid);
                  printGActinGridHeader(file_LatActinconcMap, time_bw_print, gActinGrid);
                }
                printGActinGrid(file_concMap, gActinGrid, 0);
        }

        if (arpPool)
        {
                printArpGridHeader(file_ArpConcMap, time_bw_print, arpGrid);
                printArpGrid(file_ArpConcMap, arpGrid, 0);
        }

        if (membranes[0].getExist())
        {
                printMemLenHeader(file_memLengths, time_bw_print);
                printMemLen(file_memLengths, membranes, 0);
        }

        //-------------------- Main time loop of simulation ------------------------
        int t_stepCount { 0 };
        while (currentTime < runTime)
        {
                currentTime += t_step;
                ++t_stepCount;

                // Membrane exocytosis and fluctuation
                if (membranes[0].getExist() && memExocytosis)
                {
                    membranes[0].exocytosis(t_step, stericGrid,
                                            temperature, k_exo, exo_mean,
                                            exo_stdev);
                }

                for (unsigned int i = 0; i < membranes.size(); ++i)
                {
                    if (membranes[i].getExist())
                    {
                        membranes[i].fluctuate(temperature, viscosity,
                                               t_step, actinVec,
                                               excZones, steric,
                                               timesFused, stericGrid,
                                               currentTime,
                                               branchingRegions,
                                               branching, arp23Conc,
                                               cappingRegions,
                                               antiCapRegions,
                                               nucleationRegions,
                                               sevRegions,
                                               membranes, nActin, tether,
                                               arpPool,
                                               arpGrid, memWalls,
                                               file_fusionTime,
                                               bDynamics, cortex, COMadjust,
                                               memSprings, activeBranch,
                                               activeBranchHeight, activeNuc,
                                               activeNucHeight, activeCap,
                                               activeCapHeight, activeAntiCap,
                                               activeAntiCapHeight, activeSever,
                                               activeSeverHeight);
                    }
                }

                // Cortex object fluctuation
                if (cortex.getExist())
                {
                    cortex.fluctuate(temperature, viscosity,
                                           t_step, actinVec,
                                           excZones, steric,
                                           stericGrid,
                                           currentTime,
                                           branchingRegions,
                                           branching, cappingRegions,
                                           antiCapRegions,
                                           nucleationRegions,
                                           membranes, nActin, tether,
                                           arpPool,
                                           arpGrid, memWalls,
                                           file_fusionTime,
                                           bDynamics);
                }

                // Target Brownian motion
                for (unsigned int i = 0; i < excZones.size(); ++i)
                {
                    if (excZones[i].getCircBool() && excZones[i].getMotion())
                    {
                        excZones[i].targetBM(temperature, viscosity, t_step,
                                                 membranes, excZones,
                                                 branchingRegions, nucleationRegions,
                                                 antiCapRegions, actinVec);
                    }
                }

                // Exczones directed movement
                if (currentTime >= directMoveTime)
                {
                    ExcZone::s_targetdirectedMovementCircANDTri(t_step,
                                                                firstTouch,
                                                                AFMVelVector,
                                                                indentDepth,
                                                                membranes,
                                                                excZones);
                }

                // Flow between cells in g-actin and arp2/3 grids

                if (gActinGrid.getExist() && gActinGrid.getFlow())
                {
                        gActinGrid.calcFlow(t_step);
                }

                if (arpPool && arpGrid.getFlow())
                {
                        arpGrid.calcFlow(t_step);
                }

                // latrunculin binding and unbinding
                if (latA)
                {
                        gActinGrid.latBindingMC(t_step);
                        gActinGrid.latUnbindingMC(t_step);
                }


                if (poly_b)
                {
                    if (gActinGrid.getExist())
                    {
                            polymerisationBarbedGrid(actinVec, gActinGrid, k_on_b,
                                                     t_step, currentTime, excZones,
                                                     steric, memWalls, membranes, temperature,
                                                     tether, tetherDistri, branchingRegions,
                                                     nucleationRegions, cappingRegions, sevRegions,
                                                     stericGrid, bDynamics, cortex,
                                                     nActin, arpPool, arpGrid);
                    }
                    else
                    {
                            polymerisationBarbed(actinVec, gActinConc, k_on_b,
                                                 t_step, currentTime, excZones,
                                                 steric, memWalls, membranes, temperature,
                                                 tether, tetherDistri, branchingRegions,
                                                 nucleationRegions, cappingRegions, sevRegions,
                                                 stericGrid, bDynamics, cortex,
                                                 nActin, arpPool, arpGrid);
                    }
                }


                if (poly_p)
                {
                    if (gActinGrid.getExist())
                    {
                            polymerisationPointedGrid(actinVec, gActinGrid,
                                                      k_on_p, t_step, currentTime,
                                                      excZones,
                                                      steric, memWalls, membranes, temperature,
                                                      tether, tetherDistri, branchingRegions,
                                                      nucleationRegions, cappingRegions, sevRegions,
                                                      stericGrid, bDynamics, cortex,
                                                      nActin, arpPool, arpGrid);
                    }
                    else
                    {
                            polymerisationPointed(actinVec, gActinConc, k_on_p,
                                                  t_step, currentTime, excZones,
                                                  steric, memWalls, membranes, temperature,
                                                  tether, tetherDistri, branchingRegions,
                                                  nucleationRegions, cappingRegions, sevRegions,
                                                  stericGrid, bDynamics, cortex,
                                                  nActin, arpPool, arpGrid);
                    }
                }


                if (depoly_b)
                {
                    if (gActinGrid.getExist())
                    {
                            depolymerisationBarbedGrid(actinVec, k_off_b,
                                                       t_step, tether,
                                                       nActin, dissociate,
                                                       gActinGrid, arpPool,
                                                       arpGrid, steric,
                                                       stericGrid,
                                                       excZones, memWalls,
                                                       membranes, cortex,
                                                       bDynamics,
                                                       branchingRegions,
                                                       nucleationRegions,
                                                       cappingRegions,
                                                       sevRegions);
                    }
                    else
                    {
                            depolymerisationBarbed(actinVec, k_off_b,
                                                   t_step, tether,
                                                   nActin, dissociate,
                                                   arpPool, arpGrid,
                                                   steric, stericGrid,
                                                   excZones,
                                                   memWalls, membranes,
                                                   cortex, bDynamics,
                                                   branchingRegions,
                                                   nucleationRegions,
                                                   cappingRegions,
                                                   sevRegions);
                    }
                }


                if (depoly_p)
                {
                    if (gActinGrid.getExist())
                    {
                            depolymerisationPointedGrid(actinVec, k_off_p,
                                                        t_step, tether,
                                                        nActin, dissociate,
                                                        gActinGrid, arpPool,
                                                        arpGrid, steric,
                                                        stericGrid,
                                                        excZones, memWalls,
                                                        membranes, cortex,
                                                        bDynamics,
                                                        branchingRegions,
                                                        nucleationRegions,
                                                        cappingRegions,
                                                        sevRegions);
                    }
                    else
                    {
                            depolymerisationPointed(actinVec, k_off_p,
                                                    t_step, tether,
                                                    nActin, dissociate,
                                                    arpPool, arpGrid,
                                                    steric, stericGrid,
                                                    excZones, memWalls,
                                                    membranes, cortex,
                                                    bDynamics,
                                                    branchingRegions,
                                                    nucleationRegions,
                                                    cappingRegions,
                                                    sevRegions);
                    }
                }


                if (nucleation && nucleationRegions.size() > 0)
                {
                    if (gActinGrid.getExist())
                    {
                            nucGrid(actinVec, nActin, k_nuc, t_step, currentTime,
                                    nActinLimit, nucleationRegions,
                                    excZones, membranes, memWalls, cortex,
                                    steric, stericGrid, tether, tetherDistri,
                                    gActinGrid, nucdirection, nucAngleMean,
                                    nucAngleStDev);
                    }
                    else
                    {
                            nuc(actinVec, nActin, k_nuc, t_step, gActinConc,
                                currentTime, nActinLimit, nucleationRegions,
                                excZones, membranes, memWalls, cortex,
                                steric, stericGrid, tether, tetherDistri,
                                nucdirection, nucAngleMean, nucAngleStDev);
                    }
                }


                if (branch_detach)
                        deBranch(actinVec, p_detach, nActin, arpPool, arpGrid, steric, stericGrid);



                if (branching)
                {
                    if (gActinGrid.getExist())
                    {
                        if (arpPool)
                        {
                            branchArpAndActinGridRoutine(actinVec, k_branch, t_step,
                                                         nActin,
                                                         nActinLimit, currentTime,
                                                         arpGrid, gActinGrid,
                                                         excZones,
                                                         membranes, memWalls,
                                                         cortex,
                                                         steric,
                                                         stericGrid);
                        }
                        else
                        {
                            // Function for branching with gActin grid but no arp grid
                            // Probably never want this
                            std::cout << "Error, function not included yet" << std::endl;
                            exit(1);

                        }
                    }
                    else
                    {
                        if (arpPool)
                        {
                            branchArpGridRoutine(actinVec, k_branch, t_step, gActinConc,
                                                 nActin, nActinLimit, currentTime,
                                                 arpGrid, excZones, membranes,
                                                 memWalls, cortex,
                                                 steric, stericGrid);
                        }
                        else
                        {
                            if (branchingRegions.size() == 0)
                            {
                                    branchingRoutine(actinVec, k_branch, t_step, arp23Conc, gActinConc,
                                                     nActin, nActinLimit, currentTime,
                                                     excZones, membranes,
                                                     memWalls, cortex,
                                                     steric, stericGrid,
                                                     tether, tetherDistri);
                            }
                            else
                            {
                                    branchRegionsRoutine(actinVec, k_branch, t_step,
                                                         gActinConc,
                                                         nActin,
                                                         nActinLimit, currentTime,
                                                         branchingRegions,
                                                         excZones, membranes,
                                                         memWalls, cortex,
                                                         steric, stericGrid);
                            }
                        }
                    }
                }


                if (capping_b)
                {
                        if (cappingRegions.size() == 0 && antiCapRegions.size() == 0)
                                capBarbed(actinVec, p_cap_b);
                        else if (antiCapRegions.size() == 0)
                                capBarbedRegion(actinVec, p_cap_b, cappingRegions);
                        else
                        {
                                capBarbedAntiRegion(actinVec, p_cap_b, antiCapRegions);
                        }

                }

                if (capping_p)
                {
                        if (cappingRegions.size() == 0 && antiCapRegions.size() == 0)
                                capPointed(actinVec, p_cap_p);
                        else if (antiCapRegions.size() == 0)
                                capPointedRegion(actinVec, p_cap_p, cappingRegions);
                        else
                                capPointedAntiRegion(actinVec, p_cap_p, antiCapRegions);
                }

                if (uncapping_b)
                        uncapBarbed(actinVec, p_uncap_b);

                if (uncapping_p)
                        uncapPointed(actinVec, p_uncap_p);


                if (severing)
                {
                        if (sevRegions.size() == 0)
                        {
                                severfilament(actinVec, k_sever, t_step, nActin, nActinLimit,
                                              severing, currentTime, steric, stericGrid,
                                              gActinGrid, arpPool, arpGrid);
                        }
                        else
                        {
                                   severfilamentRegions(actinVec, k_sever,
                                                        t_step, nActin,
                                                        nActinLimit, severing,
                                                        currentTime,
                                                        steric, stericGrid,
                                                        gActinGrid, arpPool,
                                                        arpGrid, sevRegions);
                        }
                }

                if (crossLinking)
                {
                    crosslinkingRoutine(actinVec, k_cLink, t_step, currentTime,
                                        steric, stericGrid, crossLinking);
                }


                if (unLinking)
                {
                    unLink(actinVec, k_unLink, t_step, currentTime, steric,
                           stericGrid, crossLinking);
                }

                if (bDynamics)
                {

                    BrownianDynamics(actinVec, excZones, memWalls,
                                     membranes, cortex, steric,
                                     stericGrid, t_step,
                                     currentTime, temperature,
                                     viscosity, tether, crossLinking);

                }

                // Allowing membrane wall to move down
                if (memWalls.size() != 0)
                {
                    double moveDistance = -memWalls[0].distToMoveWallBack2(actinVec,
                                                                           stericGrid);

                    if (memWalls[0].getYpos()+moveDistance > memWalls[0].getOrigYpos())
                    {
                        memWalls[0].moveWall(moveDistance, stericGrid);
                        memWalls[0].moveCoupledRegions(0,moveDistance,branchingRegions,
                                                       nucleationRegions, cappingRegions,
                                                       sevRegions);
                    }
                }

                if (t_stepCount % 50000 == 0)
                {
                        std::cout << "Simulation time: " << currentTime << std::endl;
                }


                if (t_stepCount % n_bw_print == 0)
                {
                        if (printToFile)
                        {

                                printActinframes(file_results, nActin, actinVec,
                                                 currentTime, t_step,
                                                 excZones, nucleationRegions,
                                                 branchingRegions,
                                                 cappingRegions, antiCapRegions,
                                                 sevRegions,
                                                 memWalls, membranes, tether,
                                                 true, cortex);
                        }

                        if (printMonoTime)
                        {
                            printMonoTimeFrames(file_mono, nActin, actinVec,
                                                currentTime, t_step);
                        }

                        printGActinGrid(file_concMap, gActinGrid, currentTime);

                        if (latA)
                        {
                                printLatAGrid(file_LatAconcMap, gActinGrid, currentTime);
                                printLatActinGrid(file_LatActinconcMap, gActinGrid, currentTime);
                        }

                        printArpGrid(file_ArpConcMap, arpGrid, currentTime);

                        if (printLen)
                        {
                                printLength(file_lengths, actinVec, currentTime);
                        }

                        if (membranes[0].getExist())
                        {
                                printMemLen(file_memLengths, membranes, currentTime);
                        }

                        if (phagoAnalysis)
                        {

                                // Want to print to file, engulfment (length of membrane touching target)
                                // and density of actin tips touching membrane

                                printEngulfment(file_engulf, currentTime, calcEngulfment(membranes, excZones));
                        }

                }

                if (membranes[0].getExist() && membranes[0].getFusable() && timesFused == excZones.size() && !fusionEnd)
                {
                    // All fusion events have taken place, start to end sim
                    runTime = currentTime + 2;
                    // have sim stop 2 seconds after
                    fusionEnd = true;
                }

                if (printLogFile && !nActinLimitBool && nActin == nActinLimit)
                {
                    // Log that the limit of number of filaments has been reached
                    fileLog << "nActin limit reached at: " << currentTime << std::endl;
                    nActinLimitBool = true;
                }


                // --------------- End of timestep assertions --------------------------
                // Turn off for speed
                /*
                assert(nActin <= nActinLimit);
                int totalSubs = 0;
                for (unsigned int i = 0; i < actinVec.size(); ++i)
                {
                    // id num is positive
                    assert(actinVec[i].getID() >= 0);
                    // id num is position in vector

                    assert(actinVec[i].getID() == i);
                    assert(actinVec[i].getNumSubs() >= 3);
                    assert(actinVec[i].getbirthTime() >= 0);

                    assert(actinVec[i].getDaughternum() >= 0);
                    if (actinVec[i].getDistToLeadBR() > actinVec[i].getNumMonomers())
                    {
                      std::cout << i << std::endl;
                    }
                    assert(actinVec[i].getDistToLeadBR() <= actinVec[i].getNumMonomers());

                    assert(actinVec[i].getDistToLeadCL() <= actinVec[i].getNumMonomers());

                    if (actinVec[i].getParentID() == -1)
                    {
                        // If it is not a branch...

                        // min length is trimer
                        assert(actinVec[i].getLength() >= Actin::s_monomerLength*Actin::s_seedSize-1e-14);
                        assert(actinVec[i].getNumMonomers() >= Actin::s_seedSize);
                    }
                    else
                    {
                        assert(actinVec[i].getBranchDir() == -1 || actinVec[i].getBranchDir() == 1);
                        assert(actinVec[i].getLength() >= Actin::s_monomerLength*Actin::s_branchSeedSize-1e-14);
                        assert(actinVec[i].getNumMonomers() >= Actin::s_branchSeedSize);
                        assert(actinVec[i].getMotherMonoID() <= actinVec[actinVec[i].getParentID()].getNumMonomers());
                    }

                    if (bDynamics && actinVec[i].getLength() > 3*Actin::s_segmentationLength)
                    {
                        if (!actinVec[i].getFlexibleBool())
                        {
                            std::cout << actinVec[i].getLength() << std::endl;
                        }
                        assert(actinVec[i].getFlexibleBool());
                    }
                    else
                    {
                        if (actinVec[i].getFlexibleBool())
                        {
                            std::cout << actinVec[i].getLength() << std::endl;
                        }
                        assert(!actinVec[i].getFlexibleBool());
                    }
                    if (actinVec[i].getLength() < 3*Actin::s_segmentationLength)
                    {
                        assert(actinVec[i].getNumSubs() == 3);
                    }
                    else
                    {
                        assert(actinVec[i].getNumSubs() > 3);
                    }
                    for (int j = 0; j < actinVec[i].getNumSubs(); ++j)
                    {
                        assert(std::isfinite(actinVec[i].getActualSubLengths()[j]));
                        assert(std::isfinite(actinVec[i].getPresSubLengths()[j]));
                        assert(std::isfinite(actinVec[i].getPoints()[j][0]));
                        assert(std::isfinite(actinVec[i].getPoints()[j][1]));
                        assert(actinVec[i].getUnitVec(j)[0] != 0 || actinVec[i].getUnitVec(j)[1] != 0);

                        totalSubs += 1;
                    }
                    assert(std::isfinite(actinVec[i].getBarbedEnd()[0]));
                    assert(std::isfinite(actinVec[i].getBarbedEnd()[1]));

                    // Assert straight ends
                    assert(fabs(actinVec[i].getUnitVec(0)[0] - actinVec[i].getUnitVec(1)[0]) < 1e-8 && fabs(actinVec[i].getUnitVec(0)[1] - actinVec[i].getUnitVec(1)[1]) < 1e-8);
                    if (fabs(actinVec[i].getUnitVec(actinVec[i].getNumSubs()-1)[0] - actinVec[i].getUnitVec(actinVec[i].getNumSubs()-2)[0]) > 1e-8 || fabs(actinVec[i].getUnitVec(actinVec[i].getNumSubs()-1)[1] - actinVec[i].getUnitVec(actinVec[i].getNumSubs()-2)[1]) > 1e-8)
                    {
                      std::cout << fabs(actinVec[i].getUnitVec(actinVec[i].getNumSubs()-1)[0] - actinVec[i].getUnitVec(actinVec[i].getNumSubs()-2)[0]) << " : " << fabs(actinVec[i].getUnitVec(actinVec[i].getNumSubs()-1)[1] - actinVec[i].getUnitVec(actinVec[i].getNumSubs()-2)[1]) << std::endl;
                      std::cout << i << std::endl;
                    }
                    assert(fabs(actinVec[i].getUnitVec(actinVec[i].getNumSubs()-1)[0] - actinVec[i].getUnitVec(actinVec[i].getNumSubs()-2)[0]) < 1e-8 && fabs(actinVec[i].getUnitVec(actinVec[i].getNumSubs()-1)[1] - actinVec[i].getUnitVec(actinVec[i].getNumSubs()-2)[1]) < 1e-8);


                    for (unsigned int j = 0; j < actinVec[i].getAvailMonoVec().size(); ++j)
                    {
                        assert(actinVec[i].getAvailMonoVec()[j] >= Actin::s_maxSpacing);
                        assert(actinVec[i].getMonoCompat(actinVec[i].getAvailMonoID(j)));
                        assert(actinVec[i].getAvailMonoID(j) < actinVec[i].getNumMonomers());
                    }

                    std::vector<int> availVec = actinVec[i].getAvailMonoVec(); // have to get this once to put it in memoery once

                    std::vector<int> unAvailMonos = actinVec[i].getUnavailMonos(actinVec);

                    std::sort(unAvailMonos.begin(), unAvailMonos.end());

                    for (int j = 0; j < actinVec[i].getNumMonomers(); ++j)
                    {
                        if (std::find(unAvailMonos.begin(), unAvailMonos.end(), j) != unAvailMonos.end())
                        {
                            if (actinVec[i].getMonoCompat(j))
                            {
                                std::cout << "i: " << i << " j1: " << j << std::endl;
                            }
                            assert(!actinVec[i].getMonoCompat(j));
                            assert(std::find(availVec.begin(), availVec.end(), j) == availVec.end());
                        }
                        else
                        {
                            //std::cout << "j2: " << j << std::endl;
                            if (!actinVec[i].getMonoCompat(j))
                            {
                                std::cout << "i: " << i << " j2: " << j << std::endl;
                            }
                            assert(actinVec[i].getMonoCompat(j));
                            assert(std::find(availVec.begin(), availVec.end(), j) != availVec.end());
                        }
                    }

                    if (actinVec[i].getNumCLinks() > 0 && actinVec[i].getLength() < 2*Actin::s_segmentationLength && actinVec[i].getParentID() == -1)
                    {
                        if (!actinVec[i].checkForMasterInStructure(actinVec))
                        {
                            std::cout << "No master in structure: " << i << std::endl;
                        }
                        //std::cout << "Checking for master: " << i << std::endl;
                        assert(actinVec[i].checkForMasterInStructure(actinVec));
                    }

                    // Structures! - branches
                    for (int j = 0; j < actinVec[i].getNumSubs(); ++j)
                    {
                        for (unsigned int k = 0; k < actinVec[i].getBranchIDSubVector(j).size(); ++k)
                        {
                            int daughtID = actinVec[i].getBranchIDSubVector(j)[k];

                            assert(actinVec[i].getStructureID() == actinVec[daughtID].getStructureID());
                            assert(actinVec[i].getCLStructureID() == actinVec[daughtID].getCLStructureID());

                            if (actinVec[daughtID].getLength() < 2*Actin::s_segmentationLength)
                            {
                                assert(actinVec[i].getSubStructureID() == actinVec[daughtID].getSubStructureID());
                            }

                            double lenAlongSub;
                            int subid = actinVec[i].findSubunit(actinVec[daughtID].getMotherMonoID(), lenAlongSub);
                            assert(subid == j);

                        }
                    }

                    // Structures - crosslinks
                    for (int j = 0; j < actinVec[i].getNumCLinks(); ++j)
                    {
                        int otherFilaID = actinVec[i].getCLinkActinAndSite(j)[2];
                        assert(actinVec[i].getStructureID() == actinVec[otherFilaID].getStructureID());
                        //assert(actinVec[i].getCLStructureID() != actinVec[otherFilaID].getCLStructureID());

                        // Check sub is correct
                        int k = actinVec[i].getCLinkActinAndSite(j)[1];
                        double lenAlongSub;
                        int subid = actinVec[i].findSubunit(actinVec[i].getCLinkActinAndSite(j)[0], lenAlongSub);
                        assert(subid == k);

                    }

                    for (unsigned int j = 0; j < membranes.size(); ++j)
                    {
                      if (membranes[j].getExist())
                      {
                        // Check filament is not outside cell membrane
                        std::array<double,2> barbedEnd;
                        barbedEnd[0] = actinVec[i].getBarbedEnd()[0];
                        barbedEnd[1] = actinVec[i].getBarbedEnd()[1];
                        if (membranes[j].getCellMemBool())
                        {
                          assert(membranes[j].checkPointIn(barbedEnd));
                        }
                        else
                        {
                          assert(!membranes[j].checkPointIn(barbedEnd));
                        }
                      }
                    }
                }
                assert (totalSubs == Actin::s_total_subunits);
                */

                // Old fashioned steric check
                /*
                if (steric)
                {
                    for (int i = 0; i < nActin-1; ++i)
                    {
                        for (int j = 0; j < actinVec[i].getNumSubs(); ++j)
                        {
                            for (int k = i; k < nActin; ++k)
                            {
                                for (int l = 0; l < actinVec[k].getNumSubs(); ++l)
                                {
                                    if (actinVec[i].getID() == actinVec[k].getID())
                                    {
                                        // Same filament, ignore same sub, and adjacent subs
                                        //std::cout << "Same filament" << std::endl;
                                        if (!actinVec[i].getFlexibleBool())
                                        {
                                            // Rigid filament, cant hinder self
                                            break;
                                        }

                                        if (j == l || j == l + 1 || j == l-1)
                                        {
                                            // ignore same or adjacent subunits while on the same filament
                                            continue;
                                        }
                                    }


                                    if (actinVec[i].getID() == actinVec[k].getParentID())
                                    {
                                         //std::cout << "Other filament parent" << std::endl;
                                         // k filament is parent, ignore sub that the branch sits on AND ADJACENT SUBS

                                        if (actinVec[i].getLength() < 2*Actin::s_segmentationLength && actinVec[k].getLength() < 2*Actin::s_segmentationLength)
                                        {
                                            // both mother and daughter rigid, cant hinder
                                            break;
                                        }


                                        if (j == actinVec[k].getBranchSubunit() || j == actinVec[k].getBranchSubunit()-1 || j == actinVec[k].getBranchSubunit()+1)
                                        {
                                            break;
                                        }
                                   }

                                    if (actinVec[i].getParentID() == actinVec[k].getID())
                                    {

                                        if (actinVec[i].getLength() < 2*Actin::s_segmentationLength && actinVec[k].getLength() < 2*Actin::s_segmentationLength)
                                        {
                                            // both mother and daughter rigid, cant hinder
                                            break;
                                        }

                                        if (l == actinVec[i].getBranchSubunit() || l == actinVec[i].getBranchSubunit()-1 || l == actinVec[i].getBranchSubunit()+1)
                                        {
                                            continue;
                                        }
                                    }

                                    //crosslinks
                                    bool cont = false;
                                    for (unsigned int m = 0; m < actinVec[i].getNumCLinks(); ++m)
                                    {
                                       if (actinVec[i].getCLinkActinAndSite(m)[2] == k && (l >= actinVec[i].getCLinkActinAndSite(m)[5]-1 && l <= actinVec[i].getCLinkActinAndSite(m)[5]+1) && (j >= actinVec[i].getCLinkActinAndSite(m)[1]-1 && j <= actinVec[i].getCLinkActinAndSite(m)[1]+1))
                                       {
                                           // Other filament is linked to this one and it is the linked sub or adjacent sub
                                           cont = true;
                                           break;
                                       }
                                   }
                                   if (cont)
                                       continue;


                                   gte::Segment<2,double> filament_1;
                                   gte::Segment<2,double> filament_2;

                                   filament_1.p[0][0] = actinVec[i].getPoints()[j][0];
                                   filament_1.p[1][0] = actinVec[i].getPoints()[j+1][0];
                                   filament_1.p[0][1] = actinVec[i].getPoints()[j][1];
                                   filament_1.p[1][1] = actinVec[i].getPoints()[j+1][1];

                                   filament_2.p[0][0] = actinVec[k].getPoints()[l][0];
                                   filament_2.p[1][0] = actinVec[k].getPoints()[l+1][0];
                                   filament_2.p[0][1] = actinVec[k].getPoints()[l][1];
                                   filament_2.p[1][1] = actinVec[k].getPoints()[l+1][1];


                                   RobustQuery min_d_GTE;
                                   auto result = min_d_GTE(filament_1, filament_2);
                                   double distance = result.distance;

                                   if (distance < (actinVec[i].getStericRadius() + actinVec[k].getStericRadius()))
                                   {
                                           std::cout << "Broke steric hind end " << distance << std::endl;
                                           std::cout << currentTime << ", " << i << " : " << j << " , " << k << " : " << l << std::endl;
                                           std::cout << actinVec[i].getLength() << ", " << actinVec[k].getLength() << std::endl;
                                           //std::cout << actinVec[i].getFlexibleBool() << ", " << actinVec[k].getFlexibleBool() << std::endl;
                                           //std::cout << actinVec[i].getBarbedEnd()[0] << ", " << actinVec[i].getBarbedEnd()[1] << std::endl;
                                           //std::cout << actinVec[k].getBarbedEnd()[0] << ", " << actinVec[k].getBarbedEnd()[1] << std::endl;
                                           //std::cout << actinVec[i].getDaughternum() << " : " << actinVec[k].getDaughternum() << std::endl;
                                           exit(1);
                                   }
                           }
                         }
                       }
                    }
                }
                */

                // ------------------------Assertions End-----------------------------
        }

        double endtime = omp_get_wtime();
        std::cout << "Simulation completed in " << (endtime - starttime);
        std::cout << " seconds" << std::endl;

        if (printLogFile)
        {
                printRunTime(fileLog, (endtime-starttime));
                if (!nActinLimitBool)
                {
                    fileLog << "nActinLimit not reached" << std::endl;
                }
        }

//        // Call python script to plot frames
//        if (pythonPlot && printToFile)
//        {
//                std::cout << "Now calling 'plotActinframes.py'" << std::endl;
//                walltime = clock();
//                calltoplotActinframes(outdir, pythonTPF, pythonFileType, pythonFrameLimits, test);
//                walltime = clock() - walltime;
//                std::cout << "Plotting frames completed in " << float(walltime)/CLOCKS_PER_SEC;
//                std::cout << " seconds" << std::endl;
//                if (makeVid)
//                {
//                        std::cout << "Now creating a video" << std::endl;
//                        callToCreateVid(outdir, fps, test);
//                }
//        }


        return 0;
}
