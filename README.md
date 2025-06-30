# ActinModel

Actin polymerisation model written by Sheffield University physics PhD student James Bradford under the supervision of Rhoda Hawkins.
Please read through all of this before doing anything.

## Getting the code

First you need to download the code and install any dependencies.

Open a terminal and create the directory where you want the code to be stored. The code has been setup to work with the directory structures in this readme, if you want to put things somewhere else or name them differently, you will need to change the code (source files: `createVid.cpp`, `pythonCaller.cpp`, `fileOutput.cpp`) and the makefile(s).

`mkdir /home/$USER/ActinModelling/`

Navigate to that directory you have just made

`cd /home/$USER/ActinModelling/`

Grab the latest version of the code off gitlab. If you do not have git you will need to install it with `sudo apt install git`

`git clone https://gitlab.com/JEBradford/actinmodel`

You will now have folder called `actinmodel` containing the code.

Two makefiles are provided. One for a local machine and one for Sheffield's high performance computing cluster ShARC.

## Dependencies

I will go through the dependencies presuming you are building on a local machine first. If you are building on ShARC you will need to download Eigen to your local machine and then upload it to your data directory. I will cover this later. 
 
### Eigen

**[Eigen](https://eigen.tuxfamily.org/)** is a linear algebra library available for free under the Mozilla Public License 2.0. Download the source code from their [website](https://eigen.tuxfamily.org/) or [gitlab](https://gitlab.com/libeigen/eigen/-/releases). Put Eigen in your home directory for a local install or in your data directory for ShARC. The easiest way is to download the zip and then extract it.

The makefile is setup to find eigen in a folder called `eigen-3.3.9` which is the latest version at the time of writing. If you have a different version please adjust the line in the makefile to point to your eigen folder, so changing `EIGEN := /home/$(USER)/eigen-3.3.9` to, for example, `EIGEN := /home/$(USER)/eigen-3.3.7`

### Boost libraries

**[Boost](https://www.boost.org/)** is a set of libraries for C++ and is licensed under the Boost Software License. Getting Boost is simple, I will assume you are using Ubuntu or another Debian based distribution. If you are not using Ubuntu/Debian you will have to look at how to get Boost using your package manager. Open a terminal

`sudo apt-get update`

`sudo apt-get install libboost-dev`

`sudo apt-get install libboost-program-options-dev`

### Anaconda Python 3

**[Anaconda](https://www.anaconda.com/products/individual) Python 3** gives you Python plus a whole suite of packages and software.

The way I install it is explained [here](https://docs.anaconda.com/anaconda/install/linux/). So go to the [website](https://www.anaconda.com/products/individual) and download the installer bash script. Make it executable by running `chmod +x yesAnaconda-3-2020.11-Linux-x86_64.sh` (the name of the script may have changed so don't just blindly copy that!). Then run the installer script

`./Anaconda-3-2020.11-Linux-x86_64.sh`

accept the terms and answer 'yes' when it asks 'Do you wish the installer to initialize Anaconda3
by running conda init?'
 
By default this installs it in the user's home directory (`home/$USER`), and the makefile is setup to find it there.

Close your terminal. Reopen your terminal and it should say '(base)' before `username@pcname:~$`.

The makefile is written for the latest version of Anaconda Python at the time of writing this which is Python 3.8. If you have a different version you will need to change the lines

`PYT := /home/$(USER)/anaconda3/include/python3.8/`

`LIB := -L $(CONPYT) -lpython3.8 -lboost_program_options -lutil`

Replacing python3.8 with your version. Sometimes this is python3.xm instead of python3.x. If you are not sure then you can check by looking in the relevant directories:

`cd /home/$USER/anaconda3/include` 

and look for a folder called `python 3.8` or `python3.7m` etc.

`cd /home/$USER/anaconda3/lib`

and look for a group of files called something like `libpython3.7m.a`, `libpython3.7m.so`, `libpython3.7m.so.1`.


### Geometric tools library

**[Geometric Tools](https://www.geometrictools.com/)** is a library covered under the Boost License. The newer versions appear to have changed the header files and does not contain the file `GTMathematics.h` which I #include. So I have provided the version I use which is version 3.8 from 2017.

**You do not need to download Geometric Tools.** You could remake `GTMathematics.h` (which is just a file that #includes all the mathematics header files) for the latest version but I see no reason to.

### FFmpeg

**[FFmpeg](https://ffmpeg.org/)** is a suite of libraries for handling video and audio. We need it to create videos from our frames. Installation is simple, again assuming you are on Ubuntu/Debian

`sudo apt-get install ffmpeg`

### Others

You will also need other stuff that you probably already have, if in doubt run the following

`sudo apt install make`

`sudo apt install g++`

## Installation on a local machine 

Once you have all the necessary dependencies open a terminal and navigate to the directory where the makefile is (`/home/$USER/ActinModelling/actinmodel`). For the first time building you need to make a directory called bin.

`mkdir bin`

Then run the makefile.

`make`

This will compile the code into the executable `model` which will be located in `bin`. Running `make debug` instead of `make` compiles with the debug option `-g` into the executable `model_debug`, this is useful for profiling.

To remove any compiled code run `make clean` or `make clean_debug` depending on how you compiled. This will remove the relevant executable and the build/build_debug directory with all the object files within it.

## Create output directories

Before running your first simulation you need to create output directories where output files, frames and videos will be saved.

We need an output directory for output files

`mkdir /home/$USER/ActinModelling/output`

a directory for the frames

`mkdir -p /home/$USER/ActinModelling/plotting/frames`

and a directory for the videos

`mkdir -p /home/$USER/ActinModelling/plotting/videos`

## Running a simulation

To run a simulation you need to supply it with a configuration file with parameter inputs. Examples of configuration files are included in the `configExamples` directory. For a list of all possible inputs run

`bin/model -h`

This gives a list of inputs together with any default values (in brackets) and a brief description.

To run with a configuration file run

`bin/model configFile.cfg`

So for running one of the examples run

`bin/model configExamples/SingleFila.cfg`

If a configuration file is not supplied, the simulation will run using the `default.cfg` as the input.

It is probably a good idea to make your own config directory to store the config files you write.

`mkdir /home/$USER/ActinModelling/actinmodel/config/`

## Output files

On completion of a simulation, output files are generated in the directory `/home/$USER/ActinModelling/output/YYYY-MM-DD-runN`, where N is the runnumber of the day (e.g. if this is the first simulation that has been executed that particular day N will be 1). The output files you will find in this directory are...
1. `results.dat`
2. `configInput.cfg`
3. `log.dat`

`results.dat` gives information on the state of the system at each frame. This includes the coordinates of the points making up the actin filaments, the state of actin filaments (capped, branched), the coordinates of the points making up the membrane, the positions of any other objects like exclusion zones, tethers or biochemical regions.

`configInput.cfg` is a copy of the config file used to run the particular simulation

`log.dat` gives details about the simulation including the random number seed.

Other output files can be produced such as the total lengths of all filaments in the system at each frame, information regarding phagocytosis. You can choose to turn these on through use of the config file.

## Visualisation

A python script `plotActinframes.py` is provided which can be set to run at the end of the simulation, printing out frames. The script takes several arguments. To see a description of the arguments run 

`plotActinframes.py -h`

You can also print frames manually using a results file as input. To do this run

`python plotActinframes.py -I results.dat -O pathToFrames/yyyy-mm-dd-runNumber -TPF 1 -FT png -FL "-5E-6, -5E-6, 5E-6, 5E-6"`

So the above will create png images from the results.dat file with a time per frame (inverse of framerate) of 1 second and the frame size limits will be plus minus five microns in x and y.

Frames can be stacked into a video using FFmpeg. A script `stackToVid.sh` is provided which does this, and this can be set to run at the end of a simulation. To manually create a video from frames as input

`stackToVid.sh fullPathToFrames/yyyy-mm-dd-runNumber fullPathToVideo 20 0`

Which creates a video in `fullPathToVideo` from the images in `fullPathToFrames/yyyy-mm-dd-runNumber`, with a framerate of 20. The last argument is to do with whether or not the simulation was a test run, but you can ignore this. **Make sure not to have forward slashes at the end of the input and output path, else it will not work.**

## Installation and running on ShARC

**For general help using the hpc services at Sheffield see [here](https://docs.hpc.shef.ac.uk/en/latest/hpc/index.html)**

Upload the actinmodel directory (which contains the code and GeometricTools) and Eigen to your data directory on ShARC. I tend to use sftp for this so something like this

`cd /home/$USER/ActinModelling/`

`sftp username@sharc.shef.ac.uk`

`cd /data/username`

`mkdir ActinModelling`

`put -r actinmodel`

`put -r eigen-3.3.9`

Boost, Anaconda Python and FFmpeg are provided as modules on ShARC, **so you should not try to upload them.**

SSH into ShARC

`ssh -X username@sharc.shef.ac.uk`

Start an interactive session

`qsh` or `qrshx`

Navigate to where the code is and rename `makefile` to something else, perhaps `makefile_desktop` and rename `makefile_ShARC` to `makefile`.

`mv makefile makefile_desktop`

`mv makefile_ShARC makefile`

`makefile_ShARC` works at the time of writing this but may need editing if ShARC update libraries/modules in the future.

Make sure you have an output directory, and plotting directorys if you are plotting frames/videos.

Make sure you create a `bin` directory as you did on your local machine.

**Make sure to run `make clean` to remove any compiled code from your desktop you may have uploaded.** 

You need to create an anaconda environment (see [here](https://docs.hpc.shef.ac.uk/en/latest/sharc/software/apps/python.html)). My conda environment is called `forActin`. To create one, in an interactive session first load python

`module load apps/python/conda`

then create your environment, we need Python 3.5 and some modules for plotting

`conda create -n forActin python=3.5 numpy matplotlib scipy`

A batch file (`batchForSharc.sh`) is provided in the gitlab repository. You can use this as a template for submitting a job. **Make sure to add your email and username for changing to the directory.** You can see that this batch file requests the necessary hardware from ShARC, loads the dependencies (anaconda and boost), activates the anaconda environment, deals with the python library path, then compiles the code and runs the simulation.

You can run this by

`qsub batchForSharc.sh`


### FFmpeg on ShARC

You can download the results file from ShARC and then produce frames and a video locally. However, if you want to do this on ShARC you can. Producing frames works fine as is but to use FFmpeg on ShARC you need to make some minor alterations.

In `src/createVid.cpp` change the line

`std::string outPath = "/home/"+boost::lexical_cast<std::string>(USER)+"/ActinModelling";`

to where your plotting directories are in ShARC, if you have followed this guide exactly then this will be `/data/$USER/ActinModelling` therefore it is 

`std::string outPath = "/data/"+boost::lexical_cast<std::string>(USER)+"/ActinModelling";`

In your batch file you need to load FFmpeg so add the lines 

```
#Load ffmpeg
module load apps/ffmpeg/4.3.2/gcc-8.2-cmake-3.17.1
```

You can see that this also loads the compiler, just like loading Boost.

## Troubleshooting

I may add more stuff here in the future.

### Python problems

For compiling and running the only thing that will probably give any problems is Python. If you have not followed the above and are using on a machine with Python/Anaconda already installed you may have to add it to your library PATH. **You should not need to do this if you have followed above and answered 'yes' to 'Do you wish the installer to initialize Anaconda3 by running conda init?' when installing Anaconda Python. Doing this may cause problems with FFMpeg**

To do this you need to put an `export` command in your .bashrc file to set it at login automatically. So first open a terminal, make sure you are in the `home/user` directory and open your .bashrc file using a text editor like gedit

`gedit .bashrc`

and copy and paste something like the following to the very bottom

```
#Add conda library to library path
export LD_LIBRARY_PATH=/home/$USER/anaconda3/lib
```

but adjusted to point to your Python lib directory.

After doing this you need to run 

`source .bashrc`

to activate it for this session.

### Adjusting the makefile

You may need to adjust the makefile to point to stuff located on your machine.


