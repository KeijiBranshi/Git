# Fixed Potential Barrier Quantum Dynamics Simulation using CUDA and MPI

This code simulates the time-evolution of a 1D quantum system, propogating towards a square potential barrier. Parallelizes the simulation using a combination of MPI and CUDA to perform the potential and kinetic propogations on a wave function (Psi).

The file, qmd_mpi.c, only uses MPI, while the file, qmd.cu, uses a combination of both. The wave function gets split into equally sized intervals along the x-axis. MPI handles the periodic boundary conditions, and CUDA actually performs the kinetic and potential propogation steps.

The graph, Kinetic-Potential-Total Energy, plots the time evolution of the system. Time is plotted on the x-axis, while the energy is plotted on the y-axis. The red line represents the potential energy, the blue line represents the kinetic energy, and the green line represents the total energy of the system.

Plot:
![Kinetic, Potential, Total Energy vs Time][plot]


[logo]: (https://github.com/KeijiBranshi/Git/blob/master/Physics/hpc_excercises/quantum_dynamics/Kinetic-Potential-Total%20Energy.png)
