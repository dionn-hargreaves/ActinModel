'''
How to run this directly from terminal given a 'results' file

python3 plottingScripts/plotActinframes.py -I output/yyyy-mm-dd-runX/results.dat -O output/yyyy-mm-dd-runX/plotting/frames -TPF 1 -FT png -BA o


'''

from __future__ import unicode_literals
from matplotlib import pyplot as plt
from os import makedirs, chdir
from scipy.spatial import ConvexHull
import numpy as np
import argparse

plt.rcParams["font.size"] = 12
plt.rcParams["font.family"] = "DejaVu Sans"

def plotframescpp(directory, time, actins, Ftype, angles, nucs, bins, reference, memFilas):

    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)


    ax.text(0.7,1.05,
        (u"Time: %.2f s"
        % (time)),
        bbox=dict(facecolor='none', edgecolor='black'),
        transform=ax.transAxes)

    prevActinID = -1
    if reference == 'o':
        for i in range(len(actins)):


            if actins[i].pCap == 1:
                if actins[i].br == 1:

                    #### actins[i].pos[0]*fac, actins[i].pos[1]*fac are the x,y coordinates of the branch first position
                    #### we want to plot the angle made relative to the radial vector (0,0 for now) - 0 is toward centre, pi is toward membrane
                    #### a.b/|a||b|
                    #### a = actins[i].pos[0], actins[i].pos[1]
                    #### b = actins[i].pos[2]-actins[i].pos[0], actins[i].pos[3]- actins[i].pos[1]
                    #### a.b = actins[i].pos[0]*(actins[i].pos[2]-actins[i].pos[0])+actins[i].pos[1]*(actins[i].pos[3]-actins[i].pos[1])
                    ##### |a||b| = (((actins[i].pos[0])**2 + (actins[i].pos[1])**2)**(1/2))*(((actins[i].pos[2]-actins[i].pos[0])**2+(actins[i].pos[3]- actins[i].pos[1])**2)**(1/2))
                    ##### store angle 
                    adotb = actins[i].pos[0]*(actins[i].pos[2]-actins[i].pos[0])+actins[i].pos[1]*(actins[i].pos[3]-actins[i].pos[1])
                    magab = (((actins[i].pos[0])**2 + (actins[i].pos[1])**2)**(1/2))*(((actins[i].pos[2]-actins[i].pos[0])**2+(actins[i].pos[3]- actins[i].pos[1])**2)**(1/2))
                        
                    theta = np.arccos(adotb/magab)
                    angles.append(TwoPoints(float(actins[i].pos[0]), float(actins[i].pos[1]), float(theta), float(time)))


    elif reference == 'm':
        Fila_x_all = [obj.pos[0] for obj in memFilas]
        Fila_x_all = np.array(Fila_x_all)
        Fila_y_all = [obj.pos[1] for obj in memFilas]
        Fila_y_all = np.array(Fila_y_all)

        
        for i in range(len(actins)):
            if actins[i].pCap == 1:
                if actins[i].br == 1:
                    
                    ### we first need to find the location of the nearest point on the membrane
                    minLength = 1.0
                    storeFila = 0
                    for j in range(len(memFilas)):
                        tempLength = (actins[i].pos[0]-Fila_x_all[j])**2 + (actins[i].pos[1]-Fila_y_all[j])**2
                        if minLength > tempLength:
                            minLength = tempLength
                            storeFila = j
                        

                    # we know the branch is closest to either segment j-1 or j (as they share the closest vertex)
                    # determine whether the midpoint of j-1 or j is closest 
                    if storeFila == 0:
                        tempLength_j1 =  (actins[i].pos[0]-((memFilas[-1].pos[2]+memFilas[-1].pos[0])/2))**2 + (actins[i].pos[1]-((memFilas[-1].pos[3]+memFilas[-1].pos[1])/2))**2
                        tempLength_j2 = (actins[i].pos[0]-((memFilas[0].pos[2]+memFilas[0].pos[0])/2))**2 + (actins[i].pos[1]-((memFilas[0].pos[3]+memFilas[0].pos[1])/2))**2
                        closestFila = min(tempLength_j1, tempLength_j2)

                        ## normal vector associated with this segment
                        ## s = (memFilas.pos[0]
                        if closestFila == tempLength_j1:
                            normalVec = np.array([memFilas[-1].pos[3]-memFilas[-1].pos[1], memFilas[-1].pos[0]-memFilas[-1].pos[2]])
                        else:
                            normalVec = np.array([memFilas[0].pos[3]-memFilas[0].pos[1], memFilas[0].pos[0]-memFilas[0].pos[2]])
                    else:
                        tempLength_j1 =  (actins[i].pos[0]-((memFilas[storeFila-1].pos[2]+memFilas[storeFila-1].pos[0])/2))**2 + (actins[i].pos[1]-((memFilas[storeFila-1].pos[3]+memFilas[storeFila-1].pos[1])/2))**2
                        tempLength_j2 = (actins[i].pos[0]-((memFilas[storeFila].pos[2]+memFilas[storeFila].pos[0])/2))**2 + (actins[i].pos[1]-((memFilas[storeFila].pos[3]+memFilas[storeFila].pos[1])/2))**2
                        closestFila = min(tempLength_j1, tempLength_j2)

                        ## normal vector associated with this segment
                        ## s = (memFilas.pos[0]
                        if closestFila == tempLength_j1:
                            normalVec = np.array([memFilas[storeFila-1].pos[3]-memFilas[storeFila-1].pos[1], memFilas[storeFila-1].pos[0]-memFilas[storeFila-1].pos[2]])
                        else:
                            normalVec = np.array([memFilas[storeFila].pos[3]-memFilas[storeFila].pos[1], memFilas[storeFila].pos[0]-memFilas[storeFila].pos[2]])

                    #### a.b/|a||b|
                    #### a = normalVec[0], normalVec[1]
                    #### b = actins[i].pos[2]-actins[i].pos[0], actins[i].pos[3]- actins[i].pos[1]
                    #### a.b = normalVec[0]*(actins[i].pos[2]-actins[i].pos[0])+normalVec[1]*(actins[i].pos[3]-actins[i].pos[1])
                    ##### |a||b| = (((normalVec[0])**2 + (normalVec[1])**2)**(1/2))*(((actins[i].pos[2]-actins[i].pos[0])**2+(actins[i].pos[3]- actins[i].pos[1])**2)**(1/2))
                    ##### store angle
                    adotb = normalVec[0]*(actins[i].pos[2]-actins[i].pos[0])+normalVec[1]*(actins[i].pos[3]-actins[i].pos[1])
                    magab = (((normalVec[0])**2 + (normalVec[1])**2)**(1/2))*(((actins[i].pos[2]-actins[i].pos[0])**2+(actins[i].pos[3]- actins[i].pos[1])**2)**(1/2))
                    theta = np.arccos(adotb/magab)
                    angles.append(TwoPoints(float(actins[i].pos[0]), float(actins[i].pos[1]), float(theta), float(time)))

    elif reference == 'a':
    
        for i in range(len(actins)):
            if actins[i].pCap == 1:
                if actins[i].br == 1:
                    
                    ### we first need to find the location of the nearest point on the membrane
                    minLength = 1.0
                    storeFila = 0
                    for j in range(len(brs)):
                        tempLength = (actins[i].pos[0]-nucs[j].corners[0])**2 + (actins[i].pos[1]-nucs[j].corners[1])**2
                        if minLength > tempLength:
                            minLength = tempLength
                            storeFila = j
                        

                    # we know the branch is closest to either segment j-1 or j (as they share the closest vertex)
                    # determine whether the midpoint of j-1 or j is closest 
                    if storeFila == 0: ## if at the first edge of the active region
                        distVec = np.array([(nucs[0].corners[0]+nucs[0].corners[6])*0.5-actins[i].pos[0], (nucs[0].corners[1]+nucs[0].corners[7])*0.5-actins[i].pos[1]])
                        
                        
                    else:
                        tempLength_j1 =  (actins[i].pos[0]-((nucs[storeFila-1].corners[6]+nucs[storeFila-1].corners[0])/2))**2 + (actins[i].pos[1]-((nucs[storeFila-1].corners[7]+nucs[storeFila-1].corners[1])/2))**2
                        tempLength_j2 = (actins[i].pos[0]-((nucs[storeFila].corners[6]+nucs[storeFila].corners[0])/2))**2 + (actins[i].pos[1]-((nucs[storeFila].corners[7]+nucs[storeFila].corners[1])/2))**2
                        closestFila = min(tempLength_j1, tempLength_j2)


                        if closestFila == tempLength_j1:
                            distVec = np.array([(nucs[storeFila-1].corners[0]+nucs[storeFila].corners[0])*0.5-actins[i].pos[0], (nucs[storeFila-1].corners[1]+nucs[storeFila].corners[1])*0.5-actins[i].pos[1]])
                        else:
                            distVec = np.array([(nucs[storeFila].corners[0]+nucs[storeFila].corners[6])*0.5-actins[i].pos[0], (nucs[storeFila].corners[1]+nucs[storeFila].corners[7])*0.5-actins[i].pos[1]])
                            
                    

                    #### a.b/|a||b|
                    #### a = distVec[0], distVec[1]
                    #### b = actins[i].pos[2]-actins[i].pos[0], actins[i].pos[3]- actins[i].pos[1]
                    #### a.b = distVec[0]*(actins[i].pos[2]-actins[i].pos[0])+distVec[1]*(actins[i].pos[3]-actins[i].pos[1])
                    ##### |a||b| = (((distVec[0])**2 + (distVec[1])**2)**(1/2))*(((actins[i].pos[2]-actins[i].pos[0])**2+(actins[i].pos[3]- actins[i].pos[1])**2)**(1/2))
                    ##### store angle
                    adotb = distVec[0]*(actins[i].pos[2]-actins[i].pos[0])+distVec[1]*(actins[i].pos[3]-actins[i].pos[1])
                    magab = (((distVec[0])**2 + (distVec[1])**2)**(1/2))*(((actins[i].pos[2]-actins[i].pos[0])**2+(actins[i].pos[3]- actins[i].pos[1])**2)**(1/2))
                    theta = np.arccos(adotb/magab)
                    angles.append(TwoPoints(float(actins[i].pos[0]), float(actins[i].pos[1]), float(theta), float(time)))
                    

                  
                    
    Angles = [obj.pos[2] for obj in angles]
    Angles = np.array(Angles)

    ax.hist(Angles, bins=bins)
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    

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

class TwoPoints():
    # This could define a rectangle, membrane filament or tether point or crosslink
    def __init__(self,x1,y1,x2,y2):
        self.pos = (x1, y1, x2, y2)

class RotRegion():
    def __init__(self,x1,y1,x2,y2,x3,y3,x4,y4):
        self.corners = (x1,y1,x2,y2,x3,y3,x4,y4)


parser = argparse.ArgumentParser(description="Plotting branch direction frames")
parser.add_argument("-I", type=str,help="Input filename")
parser.add_argument("-O",type=str,help="Output directory")
parser.add_argument("-TPF", type=str,help="Time per frame (inverse of fps)")
parser.add_argument("-FT", type=str,help="File type (only give png, eps, svg or pdf)")
parser.add_argument("-BA", type=str, help="Reference branch direction (either give o for origin, for m for normal to membrane)")


args = parser.parse_args()
filename = args.I
outdirectory = args.O
timeperframe = float(args.TPF) #0.2 # time between frames in seconds
filetype = args.FT
reference = args.BA

#print(timeperframe)
#print(filetype)
#print(framelimitsString)


angles = []
actins = []
memFilas = []
nucs= []

timeframe = -1
time = 0
framestep = 0
infile = open(filename)
bins = np.linspace(0,np.pi, 19)

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
            plotframescpp(outdirectory, time, actins, filetype, angles, nucs, bins, reference, memFilas)
        # Clear data
        del actins [:]
        del memFilas [:]
        del nucs[:]
        del angles[:]

        
    elif line.startswith("# ") or line.startswith("nActin"):
        pass

    elif line.startswith("Time (s)"):
        time = float(line.split("Time (s): ")[1])
        timeframe += 1

    elif line.startswith("!C,"):
        pass

    elif line.startswith("!R,"):
        pass

    elif line.startswith("!TR,"):
        pass

    elif line.startswith("!M,"):
        pass

    elif line.startswith("!MF,"):
        values = line.split(", ")
        memFilas.append(TwoPoints(float(values[1]), float(values[2]), float(values[3]), float(values[4])))
        
    elif line.startswith("!DT,"):
        pass

    elif line.startswith("!CX,"):
        pass

    elif line.startswith("!TP"):
        pass
    
    elif line.startswith("!CL"):
        pass

    elif line.startswith('!N,'):
        values = line.split(", ")
        nucs.append(RotRegion(float(values[1]), float(values[2]), float(values[3]), float(values[4]), float(values[5]), float(values[6]), float(values[7]), float(values[8])))


    elif line.startswith('!NC'):
        pass

    elif line.startswith('!NR'):
        pass

    elif line.startswith('!B,'):
        pass

    elif line.startswith('!BC'):
        pass

    elif line.startswith('!Cap,'):
        pass

    elif line.startswith('!CapC,'):
        pass

    elif line.startswith('!ACap,'):
        pass

    elif line.startswith('!ACapC'):
        pass

    elif line.startswith('!S,'):
        pass


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
