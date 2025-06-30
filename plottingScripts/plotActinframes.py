'''
How to run this directly from terminal given a 'results' file
python plotActinframes.py -I ../../output/yyyy-mm-dd-runnumber/results.dat -O ../../plotting/frames/yyyy-mm-dd-runnumber/direct -TPF 1 -FT png -FL "-5E-6, -5E-6, 5E-6, 5E-6"

python3 plottingScripts/plotActinframes.py -I output/yyyy-mm-dd-runX/results.dat -O output/yyyy-mm-dd-runX/plotting/frames -TPF 1 -FT png -FL "-0.75E-5, -0.5E-5, 0.75E-5, 1E-5"


'''

from __future__ import unicode_literals
from matplotlib import pyplot as plt
from os import makedirs, chdir
from scipy.spatial import ConvexHull
import numpy as np
import argparse

plt.rcParams["font.size"] = 12
plt.rcParams["font.family"] = "DejaVu Sans"

def plotframescpp(directory, time, actins, circles, rects, tris, memWalls,
                  memFilas, cortexs, tPs, cLs, nucs, brs, caps, aCaps, sevs, Ftype,
                  Fls):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.8)
    ax.set_aspect('equal', adjustable='box')

    plt.xlabel(u"x (\u03bcm)")
    plt.ylabel(u"y (\u03bcm)")

    fac = 1E6 # distance factor, so we plot in um not m
    plt.xlim([Fls[0]*fac,Fls[2]*fac])
    plt.ylim([Fls[1]*fac,Fls[3]*fac])


    plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,2))
    ax.text(0.7,1.05,
        (u"Time: %.2f s"
        % (time)),
        bbox=dict(facecolor='none', edgecolor='black'),
        transform=ax.transAxes)

    for i in range(len(circles)):
        circle = plt.Circle((circles[i].x*fac,circles[i].y*fac), circles[i].r*fac, facecolor="g", edgecolor="k", linewidth=0.5, fill=True, zorder=4)
        ax.add_artist(circle)

    for i in range(len(rects)):
        ax.plot([rects[i].pos[0]*fac,rects[i].pos[0]*fac],[rects[i].pos[1]*fac,rects[i].pos[3]*fac], color="g", linewidth=0.5, zorder=4)
        ax.plot([rects[i].pos[0]*fac,rects[i].pos[2]*fac],[rects[i].pos[3]*fac,rects[i].pos[3]*fac], color="g", linewidth=0.5, zorder=4)
        ax.plot([rects[i].pos[2]*fac,rects[i].pos[2]*fac],[rects[i].pos[3]*fac,rects[i].pos[1]*fac], color="g", linewidth=0.5, zorder=4)
        ax.plot([rects[i].pos[2]*fac,rects[i].pos[0]*fac],[rects[i].pos[1]*fac,rects[i].pos[1]*fac], color="g", linewidth=0.5, zorder=4)

    for i in range(len(tris)):
        ax.plot([tris[i].pos[0]*fac,tris[i].pos[2]*fac],[tris[i].pos[1]*fac,tris[i].pos[3]*fac], color="g", linewidth=0.5, zorder=4)
        ax.plot([tris[i].pos[2]*fac,tris[i].pos[4]*fac],[tris[i].pos[3]*fac,tris[i].pos[5]*fac], color="g", linewidth=0.5, zorder=4)
        ax.plot([tris[i].pos[4]*fac,tris[i].pos[0]*fac],[tris[i].pos[5]*fac,tris[i].pos[1]*fac], color="g", linewidth=0.5, zorder=4)


    for i in range(len(nucs)):
        if (nucs[i].__class__.__name__ == "Circle"):
            circle = plt.Circle((nucs[i].x*fac,nucs[i].y*fac), nucs[i].r*fac, color="b", fill=True, alpha=0.3, linewidth=0.0, zorder=0)
            ax.add_artist(circle)
        elif (nucs[i].__class__.__name__ == "Ring"):
            radii = [nucs[i].inR*fac, nucs[i].outR*fac]
            theta = np.linspace(0, 2*np.pi, 50, endpoint=True)
            xs = np.outer(radii, np.cos(theta)) + nucs[i].x*fac
            ys = np.outer(radii, np.sin(theta)) + nucs[i].y*fac
            xs[1,:] = xs[1,::-1]
            ys[1,:] = ys[1,::-1]

            ax.fill(np.ravel(xs), np.ravel(ys), color='b', alpha=0.3, linewidth=0.0, zorder=0)

        else:
            points = np.array([[nucs[i].corners[0]*fac, nucs[i].corners[1]*fac],
                      [nucs[i].corners[2]*fac, nucs[i].corners[3]*fac],
                      [nucs[i].corners[4]*fac, nucs[i].corners[5]*fac],
                      [nucs[i].corners[6]*fac, nucs[i].corners[7]*fac]])

            hull = ConvexHull(points)
            ax.fill(points[hull.vertices,0], points[hull.vertices,1], 'b', linewidth=0.0, alpha=0.3, zorder=0)




    for i in range(len(brs)):
        if (brs[i].__class__.__name__ == "Circle"):
            circle = plt.Circle((brs[i].x*fac,brs[i].y*fac), brs[i].r*fac, color="m", fill=True, alpha=0.3, linewidth=0.0, zorder=0)
            ax.add_artist(circle)
        else:
            points = np.array([[brs[i].corners[0]*fac, brs[i].corners[1]*fac],
                      [brs[i].corners[2]*fac, brs[i].corners[3]*fac],
                      [brs[i].corners[4]*fac, brs[i].corners[5]*fac],
                      [brs[i].corners[6]*fac, brs[i].corners[7]*fac]])

            hull = ConvexHull(points)
            ax.fill(points[hull.vertices,0], points[hull.vertices,1], 'm', linewidth=0.0, alpha=0.3, zorder=0)


    for i in range(len(caps)):
        if (caps[i].__class__.__name__ == "Circle"):
            circle = plt.Circle((caps[i].x*fac,caps[i].y*fac), caps[i].r*fac, color="g", fill=True, alpha=0.3, linewidth=0.0, zorder=0)
            ax.add_artist(circle)
        else:
            points = np.array([[caps[i].corners[0]*fac, caps[i].corners[1]*fac],
                      [caps[i].corners[2]*fac, caps[i].corners[3]*fac],
                      [caps[i].corners[4]*fac, caps[i].corners[5]*fac],
                      [caps[i].corners[6]*fac, caps[i].corners[7]*fac]])

            hull = ConvexHull(points)
            ax.fill(points[hull.vertices,0], points[hull.vertices,1], 'g', linewidth=0.0, alpha=0.3, zorder=0)

    for i in range(len(aCaps)):
        if (aCaps[i].__class__.__name__ == "Circle"):
            circle = plt.Circle((aCaps[i].x*fac,aCaps[i].y*fac), aCaps[i].r*fac, color="g", fill=True, alpha=0.3, linewidth=0.0, zorder=0)
            ax.add_artist(circle)
        else:
            points = np.array([[aCaps[i].corners[0]*fac, aCaps[i].corners[1]*fac],
                      [aCaps[i].corners[2]*fac, aCaps[i].corners[3]*fac],
                      [aCaps[i].corners[4]*fac, aCaps[i].corners[5]*fac],
                      [aCaps[i].corners[6]*fac, aCaps[i].corners[7]*fac]])
            hull = ConvexHull(points)
            ax.fill(points[hull.vertices,0], points[hull.vertices,1], 'g', linewidth=0.0, alpha=0.3, zorder=0)

    for i in range(len(sevs)):

        if (sevs[i].__class__.__name__ == "Circle"):
            circle = plt.Circle((sevs[i].x*fac,sevs[i].y*fac), sevs[i].r*fac, color="r", fill=True, alpha=0.3, linewidth=0.0, zorder=0)
            ax.add_artist(circle)
        else:
            points = np.array([[sevs[i].corners[0]*fac, sevs[i].corners[1]*fac],
                      [sevs[i].corners[2]*fac, sevs[i].corners[3]*fac],
                      [sevs[i].corners[4]*fac, sevs[i].corners[5]*fac],
                      [sevs[i].corners[6]*fac, sevs[i].corners[7]*fac]])
            hull = ConvexHull(points)
            ax.fill(points[hull.vertices,0], points[hull.vertices,1], 'r', linewidth=0.0, alpha=0.3, zorder=0)




    for i in range(len(memWalls)):
        ax.plot([(-memWalls[i].length*fac)/2, (memWalls[i].length*fac)/2], [memWalls[i].y*fac,memWalls[i].y*fac],
                color = "r", zorder=5)

    for i in range(len(memFilas)):
        ax.plot([memFilas[i].pos[0]*fac, memFilas[i].pos[2]*fac], [memFilas[i].pos[1]*fac, memFilas[i].pos[3]*fac],
                color = "r", lw=0.5, zorder=5)

        #ax.scatter([memFilas[i].pos[0]*fac, memFilas[i].pos[2]*fac], [memFilas[i].pos[1]*fac, memFilas[i].pos[3]*fac],
        #            s=20, marker="x", c="r")

    for i in range(len(cortexs)):
        ax.plot([cortexs[i].pos[0]*fac, cortexs[i].pos[2]*fac], [cortexs[i].pos[1]*fac, cortexs[i].pos[3]*fac],
                color = "#ADD8E6", lw=3, zorder=5)


    '''
    for i in range(len(tPs)):
        #Plot red x for the fixed tether point
        ax.plot(tPs[i].pos[0]*fac, tPs[i].pos[1]*fac, "rx", markersize=3, zorder=6)
        # and a blue cross for one on filament
        ax.plot(tPs[i].pos[2]*fac, tPs[i].pos[3]*fac, "bx", markersize=3, zorder=6)
    '''

    for i in range(len(cLs)):
        #Plot green line for the crosslink
        ax.plot([cLs[i].pos[0]*fac,cLs[i].pos[2]*fac],[cLs[i].pos[1]*fac,cLs[i].pos[3]*fac],
                    color = 'g', lw = 1.5, zorder=6)

        #Plot green blob at centre
        #ax.plot(cLs[i].pos[0]*fac, cLs[i].pos[1]*fac, "go", markersize=2, zorder=6)
        #ax.plot(cLs[i].pos[2]*fac, cLs[i].pos[3]*fac, "go", markersize=2, zorder=6)
        centreX = (cLs[i].pos[0]+0.5*(cLs[i].pos[2]-cLs[i].pos[0]))*fac
        centreY = (cLs[i].pos[1]+0.5*(cLs[i].pos[3]-cLs[i].pos[1]))*fac
        ax.plot(centreX, centreY, "go", markersize=2, zorder=6)


    prevActinID = -1
    for i in range(len(actins)):
        ax.plot([actins[i].pos[0]*fac,actins[i].pos[2]*fac],[actins[i].pos[1]*fac,actins[i].pos[3]*fac],
                    color = 'b', lw = 0.5, zorder=1)
        #ax.plot(actins[i].pos[0]*fac,actins[i].pos[1]*fac,'ko', markersize=1, zorder=2)


        if actins[i].pCap == 1:
            if actins[i].br == 1:
                # Plot a magenta circle at the start to represent arp2/3
                ax.plot(actins[i].pos[0]*fac,actins[i].pos[1]*fac,'mo', markersize=2, zorder=3)
            else:
                # Plot a green circle at start to represent cap
                ax.plot(actins[i].pos[0]*fac,actins[i].pos[1]*fac,'go', markersize=2, zorder=3)
        if actins[i].bCap == 1:
            # Plot a green circle at the end to represent cap
            ax.plot(actins[i].pos[2]*fac,actins[i].pos[3]*fac,'go', markersize=2, zorder=3)


        # For debugging, plot id number of filament
        '''
        if (actins[i].actinID != prevActinID):# If its a new filament
            ax.annotate(str(actins[i].actinID),xy=(actins[i].pos[0]*fac,actins[i].pos[1]*fac))
            ax.plot(actins[i].pos[0]*fac,actins[i].pos[1]*fac,'ko', markersize=1, zorder=2)
            prevActinID = actins[i].actinID
        '''

    makedirs(directory, exist_ok=True) # check to see if directory exists,
    #if it doesn't we create it here
    #plt.legend(loc='lower right',prop={'size':6})
    if (Ftype == "png"):
        plt.savefig("%s/time_%.2f_ms.png"
            %(directory,time * 1000),dpi=300)

    elif (Ftype == "svg"):
        plt.savefig("%s/time_%.2f_ms.svg"
            %(directory,time * 1000),format='svg')

    elif (Ftype == "eps"):
        plt.savefig("%s/time_%.2f_ms.eps"
            %(directory,time * 1000),format='eps')

    elif (Ftype == "pdf"):
        plt.savefig("%s/time_%.2f_ms.pdf"
            %(directory,time * 1000),format='pdf')

    else:
        print("Invalid filetype given, not plotting")

    fig.clf()
    plt.close()

class Actin_sub():
    def __init__(self,ID,subID,x1,y1,x2,y2,pCap,bCap,br):
        self.actinID = ID # NEW
        self.subID = subID
        self.pos = (x1,y1,x2,y2)
        self.pCap = pCap
        self.bCap = bCap
        self.br = br

class Circle():
    def __init__(self,x,y,r):
        self.x = x
        self.y = y
        self.r = r

class Ring():
    def __init__(self,x,y,outR,inR):
        self.x = x
        self.y = y
        self.outR = outR
        self.inR = inR

class TwoPoints():
    # This could define a rectangle, membrane filament or tether point or crosslink
    def __init__(self,x1,y1,x2,y2):
        self.pos = (x1, y1, x2, y2)

class ThreePoints():
    # This could define a triangle
    def __init__(self,x1,y1,x2,y2, x3, y3):
        self.pos = (x1, y1, x2, y2, x3, y3)

class MemWall():
    def __init__(self,y,length):
        self.y = y
        self.length = length

class RotRegion():
    def __init__(self,x1,y1,x2,y2,x3,y3,x4,y4):
        self.corners = (x1,y1,x2,y2,x3,y3,x4,y4)

parser = argparse.ArgumentParser(description="Plotting actin frames")
parser.add_argument("-I", type=str,help="Input filename")
parser.add_argument("-O",type=str,help="Output directory")
parser.add_argument("-TPF", type=str,help="Time per frame (inverse of fps)")
parser.add_argument("-FT", type=str,help="File type (only give png, eps, svg or pdf)")
parser.add_argument("-FL", type=str,help="Frame limits (xmin,ymin,xmax,ymax) (m)")

args = parser.parse_args()
filename = args.I
outdirectory = args.O
timeperframe = float(args.TPF) #0.2 # time between frames in seconds
filetype = args.FT
framelimitsString = args.FL
#print(timeperframe)
#print(filetype)
#print(framelimitsString)

frameLimits = [x.strip() for x in framelimitsString.split(',')]
frameLimits = list(map(float, frameLimits))

actins = []
circles = []
rects = []
tris = []
memWalls = []
memFilas = []
cortexs = []
tPs = []
cLs = []
nucs = []
brs = []
caps = []
aCaps = []
sevs = []

timeframe = -1
time = 0
framestep = 0
infile = open(filename)

for line in infile:
    timestep = float(line.split("Time between frames (s): ")[1])
    framestep = max(int(timeperframe / timestep),1)
    break

#next(infile) # Skip the next line

for line in infile:
    if line.startswith("#-"):
        # Reached end of the frame
        if (timeframe % framestep) == 0:
            # Plot this frame
            plotframescpp(outdirectory, time, actins, circles, rects, tris,
                          memWalls, memFilas, cortexs, tPs, cLs, nucs, brs, caps,
                          aCaps, sevs, filetype, frameLimits)
        # Clear data
        del actins [:]
        del circles [:]
        del rects [:]
        del tris [:]
        del memWalls [:]
        del memFilas [:]
        del cortexs [:]
        del tPs [:]
        del cLs [:]
        del nucs [:]
        del brs [:]
        del caps [:]
        del aCaps [:]
        del sevs [:]

    elif line.startswith("# ") or line.startswith("nActin"):
        pass

    elif line.startswith("Time (s)"):
        time = float(line.split("Time (s): ")[1])
        timeframe += 1

    elif line.startswith("!C,"):
        values = line.split(", ")
        circles.append(Circle(float(values[1]), float(values[2]), float(values[3])))

    elif line.startswith("!R,"):
        values = line.split(", ")
        rects.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))

    elif line.startswith("!TR,"):
        values = line.split(", ")
        tris.append(ThreePoints(float(values[1]), float(values[2]), float(values[3]), float(values[4]), float(values[5]), float(values[6])))

    elif line.startswith("!M,"):
        values = line.split(", ")
        memWalls.append(MemWall(float(values[1]), float(values[2])))

    elif line.startswith("!MF,"):
        values = line.split(", ")
        memFilas.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))

    elif line.startswith("!CX,"):
        values = line.split(", ")
        cortexs.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))

    elif line.startswith("!TP"):
        values = line.split(", ")
        tPs.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))

    elif line.startswith("!CL"):
        values = line.split(", ")
        cLs.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))

    elif line.startswith('!N,'):
        values = line.split(", ")
        #nucs.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))
        nucs.append(RotRegion(float(values[1]), float(values[2]), float(values[3]), float(values[4]), float(values[5]), float(values[6]), float(values[7]), float(values[8])))

    elif line.startswith('!NC'):
        values = line.split(", ")
        #nucs.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))
        nucs.append(Circle(float(values[1]), float(values[2]), float(values[3])))

    elif line.startswith('!NR'):
        values = line.split(", ")
        #nucs.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))
        nucs.append(Ring(float(values[1]), float(values[2]), float(values[4]), float(values[3])))

    elif line.startswith('!B,'):
        values = line.split(", ")
        brs.append(RotRegion(float(values[1]), float(values[2]), float(values[3]), float(values[4]), float(values[5]), float(values[6]), float(values[7]), float(values[8])))
        #brs.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))

    elif line.startswith('!BC'):
        values = line.split(", ")
        #nucs.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))
        brs.append(Circle(float(values[1]), float(values[2]), float(values[3])))

    elif line.startswith('!Cap,'):
        values = line.split(", ")
        #caps.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))
        caps.append(RotRegion(float(values[1]), float(values[2]), float(values[3]), float(values[4]), float(values[5]), float(values[6]), float(values[7]), float(values[8])))

    elif line.startswith('!CapC,'):
        values = line.split(", ")
        #caps.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))
        caps.append(Circle(float(values[1]), float(values[2]), float(values[3])))

    elif line.startswith('!ACap,'):
        values = line.split(", ")
        #caps.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))
        aCaps.append(RotRegion(float(values[1]), float(values[2]), float(values[3]), float(values[4]), float(values[5]), float(values[6]), float(values[7]), float(values[8])))

    elif line.startswith('!ACapC'):
        values = line.split(", ")
        #caps.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))
        aCaps.append(Circle(float(values[1]), float(values[2]), float(values[3])))

    elif line.startswith('!S,'):
        values = line.split(", ")
        #caps.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))
        sevs.append(RotRegion(float(values[1]), float(values[2]), float(values[3]), float(values[4]), float(values[5]), float(values[6]), float(values[7]), float(values[8])))


    else:
        values = line.split(", ")

        actins.append(Actin_sub(int(values[0]), int(values[1]), float(values[3]), float(values[4]),
                      float(values[5]), float(values[6]), int(values[7]),
                      int(values[8]), int(values[9])))

        '''
        actins.append(Actin_sub(int(values[1]), float(values[3]), float(values[4]),
                      float(values[5]), float(values[6]), int(values[7]),
                      int(values[8]), int(values[9])))
        '''

        # Old version results.dat format with angles

        #comment this out
        '''
        actins.append(Actin_sub(int(values[1]), float(values[3]), float(values[4]),
                      float(values[5]), float(values[6]), int(values[9]),
                      int(values[10]), int(values[11]), int(values[12]),
                      int(values[13])))
        '''
