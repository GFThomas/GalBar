#!/bin/sh

####################################################
## Install the GalBar code in the NEMO repository ##
## This assume that Nemo is already installed and ##
## activated                                      ##
## By G. THOMAS @ IAC.es March 2023               ##
####################################################



compiler="g++"

current_path=`pwd`

# Create the repo thomas
mkdir ${NEMO}/usr/thomas

# Copy the GalBar file in Nemo
cp GalBar.cc ${NEMO}/usr/thomas/.


# Compile it
echo "Compile GalBar"
cd ${NEMO}/usr/thomas
${compiler} -o ${NEMO}/obj/acc/GalBar.so GalBar.cc -Iinc/ -Iinc/utils/ -I${NEMO}/inc -I${NEMO}/nemo/lib -D_FILE_OFFSET_BITS=64  -mfpmath=sse -mpreferred-stack-boundary=4 --param inline-unit-growth=50 -ggdb3 -Wall -Wextra -Winit-self  -Wshadow -Wno-format-security -O2 -fPIC -funroll-loops -fforce-addr   -Woverloaded-virtual   -Llib/ -lfalcON -L../utils/lib -lWDutils -L${NEMO}/lib -lnemo -ldl -shared -w
mv acc/GalBar.so ${NEMO}/obj/acc/.

# Use mknemo to load the GalBar potential in GyrfalcON
echo "Do mknemo gyrfalcON"
mknemo gyrfalcON

# Back to the current path
cd ${current_path}
