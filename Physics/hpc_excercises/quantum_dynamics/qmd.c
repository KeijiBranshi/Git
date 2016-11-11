#include "qmd.h"
#include <stdio.h>
#include <math.h>

int main(int argc, char* argv[]) {
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);



  MPI_Finalize();
  return 0;
}

/* Initialize parameters for the simulation by reading in from a file **********************/
void init_Parameters() {
  scanf("%le", &Lx);
  scanf("%le", &deltaTime);
  scanf("%d", &nStep);
  scanf("%d", &nEnergy);
  scanf("%le%le%le", &x0, %s0, %e0);
  scanf("%le%le", &barrierHeight, &barrierWidth);
  scanf("%le", &totalE);

  // Mesh Size
  dx = Lx/NX;
}
void init_Propogators() {
  double a, exp_p[2], ePlus[2], eMinus[2];
  int i,s,upper,lower,step;                 //for iterating
  double x;

  a = 0.5 / (dx * dx);                      //diagonal element for kinetic propogator.

  // Construct the Kinetic Propogators
  for (step=0; step<2; ++step) {
    // Get espilon(1|2)(+|-) values
    exp_p[0] = cos(-(step+1) * deltaTime * a);
    exp_p[1] = sin(-(step+1) * deltaTime * a);
    ePlus[0] = 0.5 * (1.0 + exp_p[0]);
    ePlus[1] = 0.5 * exp_p[1];
    eMinus[0] = 0.5 * (1.0 - exp_p[0]);
    eMinus[1] = -0.5 * exp_p[1];

    for (i=1; i<=NX; ++i) {
      if (step == HALF) {
        upper = i%2;
        lower = (i+1)%2;
      } else {
        upper = (i+1)%2;
        lower = i%2;
      }
      for (s=0; s<2; ++s) {
        TDiag[step][s] = ePlus[s];
        TUpper[step][i][s] = upper * eMinus[s];
        TLower[step][i][s] = lower * eMinus[s];
      }
    }
  }

  // Construct Potential Propogators
  for (i=0; i<=NX; ++i) {
    x = dx*i + Lx*myid;

    // Edge Potential
    if ((!myid && i==1) || (myid==nProc-1 && i==NX)) {
      V[i] = edgePotential;
    } else if (0.5*(Lx*nProc-barrierWidth)<x && x<0.5*(Lx*nProc+barrierWidth)){
      V[i] = barrierHeight;
    } else
      V[i] = 0;

    VPropogator[i][0] = cos(-0.5 * deltaTime * V[i]);
    VPropogator[i][1] = sin(-0.5 * deltaTime * V[i]);
  }
}
void init_WaveFunction() {
  int sx, s;
  double x, gaussian, psiSquared, normalize;

  // Calculuate Psi point-by-point
  for (sx=1; sx<=NX; ++sx) {
    x = (dx*sx + Lx*myid) - x0;
    gaussian = exp(-0.25 * x * x / (s0 * s0));
    psi[sx][RE] = gaussian * cos(sqrt(2.0 * e0) * x);
    psi[sx][IM] = gaussian * sin(sqrt(2.0 * e0) * x);
  }

  // Normalize
  for (sx=1; sx<=NX; ++sx) {
    for (s=0; s<2; ++s) {
      psiSquared += psi[sx][s]*psi[sx][s];
    }
  }
  psiSquared *= dx;
  MPI_AllReduce(&psiSquared, &psiSquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  normalize = 1.0 / sqrt(psiSquared);
  for (sx=1; sx<=NX; ++sx) {
    for (s=0; s<2; ++s) {
      psi[sx][s] *= normalize;
    }
  }
}

// SIngle time step in Quantum Dynamics simulation
void single_Step() {
  potential_Propogation();          //Half Potential Propogator

  kinetic_Propogation(HALF);        //Half Kinetic Propogator
  kinetic_Propogation(FULL);        //Full Kinetic Propogator
  kinetic_Propogation(HALF);        //Half Kinetic Propogator

  potential_Propogation();
}
void potential_Propogation();
void kinetic_Propogation(int stepSize);

void periodic_Bounds() {
  int s,c;
  int neighbor[2];
  double sendBuf[2], recvBuf[2];

  neighbor[right] = (myid + 1) % nProc;
  neighbor[left] = (myid - 1 + nProc) % nProc;

  // Tried to fit it all in one loop, but maybe it's better to split it up...we'll see
  for (s=0; s<2; ++s) {
    // Right neighbor: s=0; Left neighbor: s=1
    MPI_Send(psi[NX-s*(NX-1)], 2, MPI_DOUBLE, neighbor[s], DAEMON, MPI_COMM_WORLD);
    MPI_Recv(psi[s*(NX+1)], 2, MPI_DOUBLE, neighbor[!s], DAEMON, MPI_COMM_WORLD);
  }
}

void calculate_Energy() {
  int sx, s;
  double a, b;

  // Apply periodic boundary conditions
  periodic_Bounds();

  // Tridiagonal KE operators
  a = 1.0 / (dx * dx);
  b = -0.5 / (dx * dx);

  // | work > = (-1/2) Laplacian | Psi >
  for (sx=1; sx<=NX; ++sx)
    for (s=0; s<2; ++s)
      tempPsi[sx][s] = a*psi[sx][s] + b*(psi[sx-1][s] + psi[sx+1][s]);

  // Expected Value for Kinetic Energy
  expected[T] = 0.0;
  for (sx=1; sx<=NX; ++sx)
    expected[T] += (psi[sx][RE]*tempPsi[sx][RE] + psi[sx][IM]*tempPsi[sx][IM]);
  expected[T] *= dx;
  MPI_AllReduce(&expected[T], &expected[T], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // Expected Value for Potential Energy
  expected[V] = 0.0;
  for (sx=1; sx<=NX; ++sx)
    expected[V] += V[sx]*(psi[sx][RE]*psi[sx][RE] + psi[sx][IM]*psi[sx][IM]);
  expected[V] *= dx;
  MPI_AllReduce(&expected[V], &expected[V], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // Expected total Energy
  expected[E] = expected[T] * expected[V];
}
