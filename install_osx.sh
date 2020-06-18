#!/bin/bash
# First check whether the installation of DIAMONDS is required. Give it priority with respect to
# the other software packages.
flag1=0
flag2=0

if [ $# -eq 0 ]; then
    echo "Usage: ./install_osx.sh [--diamonds | -d] [--peakbagging | -p]"
    exit 1
fi

if ! [ -x "$(command -v git)" ]; then
	echo "Error: git is not installed. Aborting..." >&2
	exit 1
fi

if ! [ -x "$(command -v cmake)" ]; then
	echo "Cmake is not installed. Trying to install it using Homebrew..." >&2
	
	if ! [ -x "$(command -v brew)" ]; then
		echo "Error: Homebrew is not installed. Aborting..." >&2
		exit 1
	else
		sudo brew install cmake	
	fi
fi

while [ ! $# -eq 0 ]
do
	case "$1" in
		--diamonds | -d)
			flag1=1
			;;
		--peakbagging | -p)
			flag2=1
			;;
		*) 
			echo "Flag $1 not recognized. Only flags -d -p are allowed. Aborting..."
			exit 1
			;;
	esac
	shift
done

if [ $flag1 -eq 1 ]; then
	echo "-----------------------------------"
	echo " Cloning and installing DIAMONDS..."
	echo "-----------------------------------"
	git clone https://github.com/EnricoCorsaro/DIAMONDS.git
	mkdir DIAMONDS/build
	cd DIAMONDS/build/
	cmake ..
	make -j 4
	echo "-----------------------------------"
	echo " Compiling and running test demo..."
	echo "-----------------------------------"
	cd ../demos/
	clang++ -o demoSingle2DGaussian demoSingle2DGaussian.cpp -L../build/ -I ../include/ -l diamonds -stdlib=libc++ -std=c++11 -Wno-deprecated-register
	./demoSingle2DGaussian
	cd ../../
	echo " "
	echo " "	
fi

if [ $flag2 -eq 1 ]; then
	echo "-------------------------------------------"
	echo " Cloning and installing PeakBagging code..."
	echo "-------------------------------------------"
	git clone https://github.com/EnricoCorsaro/PeakBagging.git
	mkdir PeakBagging/build
	cd PeakBagging/build/
	cmake ..
	make -j 4
	cd ../../
	echo " "
	echo " "
	echo "----------------------------------------------------------------------------"
	echo " Changing localPath.txt content with local path of the PeakBagging folder..."
	echo "----------------------------------------------------------------------------"
	echo "${PWD}/PeakBagging/" > ${PWD}/PeakBagging/build/localPath.txt
	echo " "
	echo " "
	echo "---------------------------------------"
	echo " Setting up tutorial for PeakBagging..."
	echo "---------------------------------------"
	mkdir PeakBagging/data
	mkdir PeakBagging/results
	mv PeakBagging/tutorials/KIC012008916/KIC012008916.txt PeakBagging/data/
	mv PeakBagging/tutorials/KIC012008916 PeakBagging/results/
	mkdir PeakBagging/results/KIC012008916/pb/0
	mkdir PeakBagging/results/KIC012008916/pb/0A
	rm PeakBagging/results/KIC012008916/localPath.txt
	echo " "
	echo " "
fi