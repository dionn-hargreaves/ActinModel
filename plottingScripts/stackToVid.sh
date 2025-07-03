#!/bin/bash

# Bash script that creates a video using ffmpeg given an image directory
# Needs to read the file, create a copy, rename to something simply like %d.png
# then after doing that to all files, run ffmpeg
# then delete the tmp files

# 1st argument is the input directory (full path)
# 2nd argument is the output filename (full path)
# 3rd argument is the fps
# 4th argument is whether it is a test or not (1 = test, 0 = not test)

# An example running from the terminal is...
# plottingScripts/stackToVid.sh /...Full Path.../output/yyyy-mm-dd-runX/plotting/frames /...Full Path.../output/yyyy-mm-dd-runX/plotting/movie 20 0

mkdir $1/tmp
count=0
for file in `ls $1/*.png | sort -V`; do
	echo "$file"
	cp $file $1/tmp/$count.png
	#convert $file tmp/$count.jpg
	count=$((count + 1))
done

if [ "$4" == "0" ]
then
	ffmpeg -framerate $3 -pattern_type sequence -i ''"$direct"'/'"$1"'/tmp/%d.png' -vf scale=1200:-1 -c:v libx264 -crf 24 -preset veryslow -pix_fmt yuv420p $2.mp4
else
	ffmpeg -framerate $3 -pattern_type sequence -i ''"$direct"'/'"$1"'/tmp/%d.png' -vf scale=1200:-1 -c:v libx264 -crf 24 -preset veryslow -pix_fmt yuv420p $2.mp4
fi
rm -r $1/tmp
