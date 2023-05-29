### requirement

1. VTK 8
2. CMake 2.8.12.1 or newer
3. MacOS 10.9 or newer or Linux with kernel 4.12 or newer
4. ffmpeg 4 or newer

### Usage

for Linux, CmakeList.txt need to remove MACOS_BUNDLE in line add_executable

`cmake.`

`make`

`./project1F`

this program will generate 1000 images with different camera angles.

to encode these images to video

`ffmpeg -r 1/5 -i frame%03d.png -c:v libx264 -vf fps=24 -pix_fmt yuv420p out.mp4`

 

