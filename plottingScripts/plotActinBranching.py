'''
'''
How to run this directly from terminal given a 'results' file
 
 PENDING
python3 plottingScripts/plotActinBranching.py -I output/yyyy-mm-dd-runX/results.dat -O output/yyyy-mm-dd-runX/plotting/frames -TPF 1 -FT png -FL "-0.75E-5, -0.5E-5, 0.75E-5, 1E-5" -TB 1 -TC 0

python3 plottingScripts/plotActinBranching.py -I output/UseableFullOutputs/2025-08-19-FullPhago1um/results.dat -O output/UseableFullOutputs/2025-08-19-FullPhago1um/plotting/branching -TPF 1 -FT png -FL "-0.75E-5, -0.5E-5, 0.75E-5, 1E-5" -TB 1 -TC 0

'''

from __future__ import unicode_literals
from matplotlib import pyplot as plt
from os import makedirs, chdir
from scipy.spatial import ConvexHull
import numpy as np
import argparse
from matplotlib import colormaps
from matplotlib.colors import Normalize


plt.rcParams["font.size"] = 12
plt.rcParams["font.family"] = "DejaVu Sans"

def plotframe(directory, time, actins, circles, rects, tris, memWalls,
                  memFilas, cortexs, tPs, cLs, nucs, brs, caps, aCaps, sevs, Ftype,
                  Fls, toggleBranch, toggleCap, loggedBranch, cmap, norm, fig, ax):

    fac = 1E6 # distance factor, so we plot in um not m
    plt.xlim([Fls[0]*fac,Fls[2]*fac])
    plt.ylim([Fls[1]*fac,Fls[3]*fac])


    plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,2))
    
    use_col = cmap(norm(time))
    
    for i in range(len(circles)):
        circle = plt.Circle((circles[i].x*fac,circles[i].y*fac), circles[i].r*fac, facecolor="g", edgecolor="k", linewidth=0.5, fill=True)
        ax.add_artist(circle)

    for i in range(len(rects)):
        ax.plot([rects[i].pos[0]*fac,rects[i].pos[0]*fac],[rects[i].pos[1]*fac,rects[i].pos[3]*fac], color="g", linewidth=0.5)
        ax.plot([rects[i].pos[0]*fac,rects[i].pos[2]*fac],[rects[i].pos[3]*fac,rects[i].pos[3]*fac], color="g", linewidth=0.5)
        ax.plot([rects[i].pos[2]*fac,rects[i].pos[2]*fac],[rects[i].pos[3]*fac,rects[i].pos[1]*fac], color="g", linewidth=0.5)
        ax.plot([rects[i].pos[2]*fac,rects[i].pos[0]*fac],[rects[i].pos[1]*fac,rects[i].pos[1]*fac], color="g", linewidth=0.5)

    for i in range(len(tris)):
        ax.plot([tris[i].pos[0]*fac,tris[i].pos[2]*fac],[tris[i].pos[1]*fac,tris[i].pos[3]*fac], color="g", linewidth=0.5)
        ax.plot([tris[i].pos[2]*fac,tris[i].pos[4]*fac],[tris[i].pos[3]*fac,tris[i].pos[5]*fac], color="g", linewidth=0.5)
        ax.plot([tris[i].pos[4]*fac,tris[i].pos[0]*fac],[tris[i].pos[5]*fac,tris[i].pos[1]*fac], color="g", linewidth=0.5)



    for i in range(len(memWalls)):
        ax.plot([(-memWalls[i].length*fac)/2, (memWalls[i].length*fac)/2], [memWalls[i].y*fac,memWalls[i].y*fac],
                color = use_col, alpha = 0.5)

    for i in range(len(memFilas)):
        ax.plot([memFilas[i].pos[0]*fac, memFilas[i].pos[2]*fac], [memFilas[i].pos[1]*fac, memFilas[i].pos[3]*fac],
                color = use_col, lw=0.5, alpha = 0.5)

    prevActinID = -1
    for i in range(len(actins)):
        if toggleBranch:
            if actins[i].br == 1:
                numOldBranches = len(loggedBranch) ## Need to only plot new ID branches
                loggedBranch.add(actins[i].actinID)
                if numOldBranches != len(loggedBranch):
                    numOldBranches = len(loggedBranch)
                    # Plot a circle at the start to represent arp2/3
                    ax.plot(actins[i].pos[0]*fac,actins[i].pos[1]*fac,'o', markersize=2, color = use_col, alpha = 0.5)
        if toggleCap == 1:
            if actins[i].bCap == 1:
                # Plot a green circle at the end to represent cap
                ax.plot(actins[i].pos[2]*fac,actins[i].pos[3]*fac,'o', markersize=2, color = use_col, alpha = 0.5)


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
parser.add_argument("-TB", type=str,help="Toggle branch")
parser.add_argument("-TC", type=str,help="Toggle cap")

args = parser.parse_args()
filename = args.I
outdirectory = args.O
timeperframe = float(args.TPF) #0.2 # time between frames in seconds
filetype = args.FT
framelimitsString = args.FL
toggleCap = args.TC
toggleBranch = args.TB
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
loggedBranch = set()

timeframe = -1
time = 0
framestep = 0
infile = open(filename)

#next(infile) # Skip the next line
for line in infile:
    if line.startswith("Time (s)"):
        time = float(line.split("Time (s): ")[1])
    
    

## arrange colourmaps
cmap = colormaps['rainbow']  # new API, replaces cm.get_cmap
norm = Normalize(vmin=0, vmax=time)

fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.8)
ax.set_aspect('equal', adjustable='box')
    
plt.xlabel(u"x (\u03bcm)")
plt.ylabel(u"y (\u03bcm)")

time = 0

infile.close()
infile = open(filename)

for line in infile:
    timestep = float(line.split("Time between frames (s): ")[1])
    framestep = max(int(timeperframe / timestep),1)
    break
    
for line in infile:
    if line.startswith("#-"):
        # Reached end of the frame
        if (timeframe % framestep) == 0:
            # Plot this frame
            plotframe(outdirectory, time, actins, circles, rects, tris,
                          memWalls, memFilas, cortexs, tPs, cLs, nucs, brs, caps,
                          aCaps, sevs, filetype, frameLimits, toggleBranch, toggleCap, loggedBranch, cmap, norm, fig, ax)
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

makedirs(outdirectory, exist_ok=True) # check to see if directory exists,
#if it doesn't we create it here
#plt.legend(loc='lower right',prop={'size':6})
if (filetype == "png"):
    fig.savefig(outdirectory+"/Temporal_image.png",dpi=300)

elif (filetype == "svg"):
    fig.savefig(outdirectory+"/Temporal_image.svg",format='svg')

elif (filetype == "eps"):
    fig.savefig(outdirectory+"/Temporal_image.eps",format='eps')

elif (filetype == "pdf"):
    fig.savefig(outdirectory+"/Temporal_image.pdf",format='pdf')

else:
    print("Invalid filetype given, not plotting")


plt.close()
