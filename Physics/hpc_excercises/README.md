#High Performance Computing Excercises
Here is a collection of parallel programs am I either working on or have finished. Most, if not all, of these programs are meant to run on USC's High Performance Computer (some are sequential, yet to be translated into parallel).  
The programs utilize a combination of MPI and OpenMP, and I'm currently working on one that employs CUDA. If you wish to run any them on your own computer then you will need those libraries.

##Samples Include:
* Parallel Quicksort Implementation:

        Uses a sequence of hypercube MPI calls to sort a list of numbers. Built to run on 8 cores.
* Wavelet Transformation: 

        Recursively performs Daubechies wavelet transform on an image, Lenna512x512.pgm, keeping only the smooth 
        components. Theoretically can run on any perfect square number of processors, however, due to compiler issues 
        with the USC High Performance Computer, I hardcoded it to run on 16.
* Molecular Dynamics:

        Currently a bare-bones C# implementation of simple molecular dynamics simulation in a Leonard-Jones potential. 
        The C# code was used in a Data Visualization research project (stored on a private repository due to the nature
        of the research). Translating this to C++ for use on the USC HPC. 
