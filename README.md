# DimerDetection
Code used to detect possible dimers in our cryoEM data

# Prerequisites
- OS: Windows 10 
- Visual Studio (tested with VS2017)
- .NET Framework 2.7

# Installation
There are precompiled binaries in the bin/Release folder alongside neccessary .dll files for running the tool. If you want to compile the code yourself, open the .sln file in VisualStudio and build the solution

# Usage
There is a file with example Data in the Data folder. After cloning the git repository, simply extract the `.tar.gz` This is the .star file of all Warp picked particles for our described data  set, already merged with orientations from a 3D refinement of monomeric RdRp. You can run the following tools:

- `bin\Release\NearestNeighbor.exe`: calculates for all particles in a .star file the NN distances and relative orientations. Output is a tab-seperated file with all info (to be read by pandas,excel,...) and one .star file that contains rot and tilt angles between NN monomers
-  `.\bin\Release\PredictMonomers.exe`: based on monomer positions, predicts for all monomers that are not part of a dimer, the position on the micrograph where a second monomer would have to be
-  `.\bin\Release\PickDimers.exe`: Calculates all NNs and then writes coordinate files (compatible with relion extraction) of all potential dimers

A simple script `runTools.batch` is also included that contains example commands for all three tools
