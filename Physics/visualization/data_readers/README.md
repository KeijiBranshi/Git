# Data Parsing

Here is where I back up some of the code that I use for data/file parsing in my research. Right now I only have programs that help with analyzing and manipulating .xyz files (a common format found in scietific computing simulations). The two main files are:

* XYZFile.cs

        Reads and parses .xyz frames for animation in Unity. I use this in my data visualization projects to help render atoms in molecular dynamics simulations. 

* xyz_transform.rb

        I use this Ruby script to perform quick transformations on coordinates from a .xyz file. The script takes in a .xyz filepath and 6 floating point numbers as command-line arguments. The first 3 numbers are rotations in degrees along the x, y, and z axes, while the second 3 are translation amounts along those axes. The coordinates in the .xyz file are parsed into a matrix, then the user-supplied transformations are applied to that matrix. This currently only supports .xyz files with a single frame of particles.

The subdirectory, Ruby-lib, just stores some supplementary Ruby modules/functions that are useful for quick data manipulation in my research.  
