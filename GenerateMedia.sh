 #!/bin/zsh


#//  GenerateMedia.sh
#//
#//
#//  Created by mbcx4dh9 on 25/06/2025.
#//  Script generates frame images and movies from data file
# Takes input of 1. path to current folder, e.g. /Users/user/Documents/Sheffield/ActinModelling/actinmodel/
 #               2. data file date and run in the form: yyyy-mm-dd_runX

. .venv/bin/activate
echo $1
echo $2
echo $3
python3 plottingScripts/plotActinframes.py -I $1/output/$2/results.dat -O $1/output/$2/plotting/frames -TPF 1 -FT png -FL "-0.9E-5, -0.8E-5, 0.9E-5, 1.0E-5"
python3 plottingScripts/plotBarbedOrientations-arrows.py -I $1/output/$2/results.dat -O $1/output/$2/plotting/barbframes -TPF 1 -FT png -FL "-0.9E-5, -0.8E-5, 0.9E-5, 0.8E-5" -BA m
python3 plottingScripts/plotBarbOrientations-RosePlot.py -I $1/output/$2/results.dat -O $1/output/$2/plotting/barbRose -TPF 1 -FT png -BA m
plottingScripts/stackToVid.sh $1/output/$2/plotting/frames $1/output/$2/plotting/movie 20 0
plottingScripts/stackToVid.sh $1/output/$2/plotting/barbframes $1/output/$2/plotting/barb 20 0
plottingScripts/stackToVid.sh $1/output/$2/plotting/barbRose $1/output/$2/plotting/RoseBarbmovie 20 0
