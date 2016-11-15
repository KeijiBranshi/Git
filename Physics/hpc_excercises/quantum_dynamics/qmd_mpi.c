#include "qmd_mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

bool DEBUG;

int main(int argc, char* argv[]) {
  int step;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  DEBUG = false;

  init_Parameters();
  init_Variables();
  init_Propogators();
  init_WaveFunction();

  for (step=1; step<=nStep; ++step) {
    //printf("Loop Before Step: proc:%i  step=%i\n", myid, step);
    single_Step();
    if (DEBUG) printf("Loop After Step: proc:%i\n", myid);
    if (step%energyInterval == 0) {
      if (DEBUG) printf("Calculating Energy: step=%i\n", step);
      calculate_Energy();
      printf("%le %le %le %le\n", deltaTime*step, expected[_T_], expected[_V_], expected[_E_]);
    }
  }

  cleanUpVariables();
  MPI_Finalize();
  return 0;
}

/* Initialize parameters for the simulation by reading in from a file **********************/
void init_Parameters() {

  FILE* f;
  char line[MAX_LINE];

  f=fopen("qmd.in","r");

  fgets(line,MAX_LINE,f);
  sscanf(line,"%le", &Lx);
  fgets(line,MAX_LINE,f);
  sscanf(line,"%le", &deltaTime);
  fgets(line,MAX_LINE,f);
  sscanf(line,"%d", &nStep);
  fgets(line,MAX_LINE,f);
  sscanf(line,"%d", &energyInterval);
  fgets(line,MAX_LINE,f);
  sscanf(line,"%le %le %le", &x0, &s0, &e0);
  fgets(line,MAX_LINE,f);
  sscanf(line,"%le %le", &barrierHeight, &barrierWidth);
  fgets(line,MAX_LINE,f);
  sscanf(line,"%le", &edgePotential);
  fclose(f);

  /*scanf("%le", &Lx);
  scanf("%le", &deltaTime);
  scanf("%d", &nStep);
  scanf("%d", &energyInterval);
  scanf("%le %le %le", &x0, &s0, &e0);
  scanf("%le %le", &barrierHeight, &barrierWidth);
  scanf("%le", &edgePotential);*/

  /*if (false) {
    printf("Lx = %le\n", Lx);
    printf("deltaT = %le\n", deltaTime);
    printf("nStep = %d\n", nStep);
    printf("energyInterval = %d\n", energyInterval);
    printf("x0, s0, e0 = %le %le %le\n", x0, s0, e0);
    printf("bH, BW = %le %le\n", barrierHeight, barrierWidth);
    printf("edgePotential = %le\n", edgePotential);
  }*/

  // Mesh Size
  dx = Lx/NX;
}

void init_Variables() {
  int i,j;

  host_psi = (double**) malloc(MESH_SIZE * sizeof(double*));
  for (i=0; i<MESH_SIZE; ++i) {
    host_psi[i] = (double *) malloc(sizeof(double)*2);
  }

  workPsi = (double**) malloc(MESH_SIZE * sizeof(double*));
  for (i=0; i<MESH_SIZE; ++i) {
    workPsi[i] = (double *) malloc (2 * sizeof(double));
  }

  TDiag = (double**) malloc(sizeof(double*) * 2);
  for (i=0; i<2; ++i) {
    TDiag[i] = (double*) malloc(sizeof(double) * 2);
  }

  TUpper = (double***) malloc(sizeof(double**) * 2);
  for (i=0; i<2; ++i) {
    TUpper[i] = (double**) malloc(sizeof(double*) * MESH_SIZE);
    for (j=0; j<MESH_SIZE; ++j) {
      TUpper[i][j] = (double*) malloc(sizeof(double) * COMPLEX);
    }
  }

  TLower = (double***) malloc(sizeof(double**) * 2);
  for (i=0; i<2; ++i) {
    TLower[i] = (double**) malloc(sizeof(double*) * MESH_SIZE);
    for (j=0; j<MESH_SIZE; ++j) {
      TLower[i][j] = (double*) malloc(sizeof(double) * COMPLEX);
    }
  }

  VPropogator = (double**) malloc(sizeof(double*) * MESH_SIZE);
  for (i=0; i<MESH_SIZE; ++i) {
    VPropogator[i] = (double *) malloc(sizeof(double)*2);
  }

  potential = (double*) malloc(sizeof(double) * MESH_SIZE);
}

void init_Propogators() {
  double a, exp_p[2], ePlus[2], eMinus[2];
  int i,c,upper,lower,step;                 //for iterating
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

    for (c=0; c<2; ++c) {
      TDiag[step][c] = ePlus[c];
    }

    for (i=1; i<=NX; ++i) {
      if (step == HALF) {
        upper = i%2;
        lower = (i+1)%2;
      } else {
        upper = (i+1)%2;
        lower = i%2;
      }
      for (c=0; c<2; ++c) {
        //if (myid==0) printf("Loop: step=%i i=%i, c=%i\n", step, i, c);
        TUpper[step][i][c] = upper * eMinus[c];
        TLower[step][i][c] = lower * eMinus[c];
      }
    }
  }

  // Construct Potential Propogators
  for (i=1; i<=NX; ++i) {
    x = dx*i + Lx*myid;

    // Edge Potential
    if ((myid==0 && i==1) || (myid==nProc-1 && i==NX)) {
      potential[i] = edgePotential;
    } else if (0.5*(Lx*nProc-barrierWidth)<x && x<0.5*(Lx*nProc+barrierWidth)){
      potential[i] = barrierHeight;
    } else
      potential[i] = 0;

    VPropogator[i][RE] = cos(-0.5 * deltaTime * potential[i]);
    VPropogator[i][IM] = sin(-0.5 * deltaTime * potential[i]);
  }
}

void init_WaveFunction() {
  int sx, c;
  double x, gaussian, normalize;
  double psiSquared, temp;

  // Calculuate Psi point-by-point
  for (sx=1; sx<=NX; ++sx) {
    //if (myid==0) printf("<1> dx=%d\n", dx);
    x = (dx*sx + Lx*myid) - x0;
    //if (myid==0) printf("<2> dx=%d\n", dx);
    gaussian = exp (-0.25 * x * x / (s0 * s0));
    //if (myid==0) printf("<3> dx=%d\n", dx);
    host_psi[sx][RE] = gaussian * cos (sqrt(2.0 * e0) * x);
    //if (myid==0) printf("<4> dx=%d\n", dx);
    host_psi[sx][IM] = gaussian * sin (sqrt(2.0 * e0) * x);
    //if (myid==0) printf("<5> dx=%d\n", dx);
    //if (myid==0) printf("dx=%d, x=%d, gaussian=%d, host_psi[%i]=%d,%d\n", dx, x, gaussian, sx, host_psi[sx][RE], host_psi[sx][IM]);
  }

  // Normalize
  psiSquared = 0.0;
  for (sx=1; sx<=NX; ++sx) {
    for (c=0; c<2; ++c) {
      psiSquared += host_psi[sx][c]*host_psi[sx][c];
    }
  }
  //printf("Proc:%i ---- PsiSquared Before Reduce: %d\n", myid, psiSquared);
  MPI_Allreduce(&psiSquared, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  psiSquared = temp * dx;
  //printf("Proc:%i --- PSI SQUARED: %d\n", myid, psiSquared);

  normalize = 1.0 / sqrt(psiSquared);
  for (sx=1; sx<=NX; ++sx) {
    for (c=0; c<2; ++c) {
      host_psi[sx][c] *= normalize;
    }
  }
}

void cleanUpVariables() {
  int i,j;

  for (i=0; i<MESH_SIZE; ++i) {
    free(host_psi[i]);
  }
  free(host_psi);

  for (i=0; i<MESH_SIZE; ++i) {
    free(workPsi[i]);
  }
  free(workPsi);

  for (i=0; i<2; ++i) {
    free(TDiag[i]);
  }
  free(TDiag);

  for (i=0; i<2; ++i) {
    for (j=0; j<MESH_SIZE; ++j) {
      free(TUpper[i][j]);
    }
    free(TUpper[i]);
  }
  free(TUpper);

  for (i=0; i<2; ++i) {
    for (j=0; j<MESH_SIZE; ++j) {
      free(TLower[i][j]);
    }
    free(TLower[i]);
  }
  free(TLower);

  for (i=0; i<MESH_SIZE; ++i) {
    free(VPropogator[i]);
  }
  free(VPropogator);

  free(potential);
}


/**************************************************************************************
 * Time Stepping/Propogation Functions
 **************************************************************************************/

/* SIngle time step in Quantum Dynamics simulation */
void single_Step() {
  //if (myid==0 && DEBUG) printf("<1> Start Single Step\n");
  regular_Potential_Prop();          //Half Potential Propogator
  //if (myid==0 && DEBUG) printf("<2> FInished Pot Prop\n");

  regular_Kinetic_Prop(HALF);        //Half Kinetic Propogator
  //if (myid==0 && DEBUG) printf("<3> Finished Half PRop\n");
  regular_Kinetic_Prop(FULL);        //Full Kinetic Propogator
  //if (myid==0 && DEBUG) printf("<4> Finished Full Prop\n");
  regular_Kinetic_Prop(HALF);        //Half Kinetic Propogator
  //if (myid==0 && DEBUG) printf("<6> Finished Half Prop\n");

  regular_Potential_Prop();
  //if (myid==0 && DEBUG) printf("<7> Finished Single Step\n");
}

void regular_Potential_Prop() {
  int sx;
  double workPsi_Re, workPsi_Im;

  for (sx=1; sx<=NX; ++sx) {
    workPsi_Re = VPropogator[sx][RE]*host_psi[sx][RE] - VPropogator[sx][IM]*host_psi[sx][IM];
    workPsi_Im = VPropogator[sx][RE]*host_psi[sx][IM] + VPropogator[sx][IM]*host_psi[sx][RE];
    host_psi[sx][RE]= workPsi_Re;
    host_psi[sx][IM]= workPsi_Im;
  }
}

void regular_Kinetic_Prop(int t) {
  int sx, c, temp;
  double workPsi_Re, workPsi_Im;

  periodic_Bounds();

  for (sx=1; sx<=NX; ++sx) {

    workPsi_Re = TDiag[t][RE]*host_psi[sx][RE]-TDiag[t][IM]*host_psi[sx][IM];
    workPsi_Im = TDiag[t][RE]*host_psi[sx][IM]+TDiag[t][IM]*host_psi[sx][RE];
    workPsi_Re += (TLower[t][sx][RE]*host_psi[sx-1][RE]-TLower[t][sx][IM]*host_psi[sx-1][IM]);
    workPsi_Im += (TLower[t][sx][RE]*host_psi[sx-1][IM]+TLower[t][sx][IM]*host_psi[sx-1][RE]);
    workPsi_Re += (TUpper[t][sx][RE]*host_psi[sx+1][RE]-TUpper[t][sx][IM]*host_psi[sx+1][IM]);
    workPsi_Im += (TUpper[t][sx][RE]*host_psi[sx+1][IM]+TUpper[t][sx][IM]*host_psi[sx+1][RE]);

    //if (myid == 0) printf("<reg kinetic prop><4> loop: %i\n", sx);
    workPsi[sx][RE] = workPsi_Re;
    workPsi[sx][IM] = workPsi_Im;
  }

  //if (myid == 0) printf("Out of the work loop\n");
  for (sx=1; sx<=NX; ++sx) {
    for (c=0; c<2; ++c) {
      host_psi[sx][c] = workPsi[sx][c];
    }
  }
  //if (myid==0) printf("Copied work to psi\n");
}

/**************************************************************************************
 * MPI Related Functions
 *************************************************************************************/

/* Applies Periodic boundary conditions */
void periodic_Bounds() {
  int neighbor[2];
  double sendBuf[2], recvBuf[2];

  neighbor[right] = (myid + 1) % nProc;
  neighbor[left] = (myid - 1 + nProc) % nProc;

  // Send right, receive left
  sendBuf[0] = host_psi[NX][RE];
  sendBuf[1] = host_psi[NX][IM];
  MPI_Irecv(recvBuf, 2, MPI_DOUBLE, MPI_ANY_SOURCE, DAEMON, MPI_COMM_WORLD, &request);
  MPI_Send(sendBuf, 2, MPI_DOUBLE, neighbor[right], DAEMON, MPI_COMM_WORLD);
  MPI_Wait(&request, &status);
  host_psi[0][RE] = recvBuf[RE];
  host_psi[0][IM] = recvBuf[IM];

  // Send left, receive right
  sendBuf[0] = host_psi[1][RE];
  sendBuf[1] = host_psi[1][IM];
  MPI_Irecv(recvBuf, 2, MPI_DOUBLE, MPI_ANY_SOURCE, DAEMON, MPI_COMM_WORLD, &request);
  MPI_Send(sendBuf, 2, MPI_DOUBLE, neighbor[left], DAEMON, MPI_COMM_WORLD);
  MPI_Wait(&request, &status);
  host_psi[NX+1][RE] = recvBuf[RE];
  host_psi[NX+1][IM] = recvBuf[IM];
}

/* Total Energy Calculation at the current time step */
void calculate_Energy() {
  int sx, c;
  double a, b;

  // Apply periodic boundary conditions
  periodic_Bounds();

  // Tridiagonal KE operators
  a = 1.0 / (dx * dx);
  b = -0.5 / (dx * dx);

  // | work > = (-1/2) Laplacian | Psi >
  for (sx=1; sx<=NX; ++sx)
    for (c=0; c<2; ++c)
      workPsi[sx][c] = a*host_psi[sx][c] + b*(host_psi[sx-1][c] + host_psi[sx+1][c]);

  // Expected Value for Kinetic Energy
  expected[_T_] = 0.0;
  for (sx=1; sx<=NX; ++sx)
    expected[_T_] += (host_psi[sx][RE]*workPsi[sx][RE] + host_psi[sx][IM]*workPsi[sx][IM]);
  expected[_T_] *= dx;
  MPI_Allreduce(&expected[_T_], &expected[_T_], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  // Expected Value for Potential Energy
  expected[_V_] = 0.0;
  for (sx=1; sx<=NX; ++sx)
    expected[_V_] += potential[sx]*(host_psi[sx][RE]*host_psi[sx][RE] + host_psi[sx][IM]*host_psi[sx][IM]);
  expected[_V_] *= dx;
  MPI_Allreduce(&expected[_V_], &expected[_V_], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  // Expected total Energy
  expected[_E_] = expected[_T_] * expected[_V_];
}
