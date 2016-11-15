#include "qmd.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int DEBUG;

int main(int argc, char* argv[]) {
  int step;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  cudaSetDevice(myid%2);
  usingCUDA = true;
  DEBUG = false;

  init_Parameters();
  init_Variables();
  init_Propogators();
  init_WaveFunction();
  if (usingCUDA) init_GPU();

  for (step=1; step<=nStep; ++step) {
    single_Step();
    if (step%energyInterval == 0) {
      calculate_Energy();
      if (myid==0) printf("%le %le %le %le\n", deltaTime*step, expected[_T_], expected[_V_], expected[_E_]);
      //if (myid == 0) printf("%le, %le\n", deltaTime*step, expected[_V_]);
    }
  }

  cleanUpVariables();
  if (usingCUDA) GPU_Finalize();
  MPI_Finalize();
  return 0;
}

/* Initialize GPU communication and variables **********************************************/
void init_GPU() {
  size_t size = 2*(MESH_SIZE)*sizeof(double);

  cudaMalloc((double**) &dev_psi, size);
  cudaMalloc((double**) &dev_workPsi, size);
  cudaMalloc((double**) &dev_TDiag_half, 2 * 2 * sizeof(double));
  cudaMalloc((double**) &dev_TDiag_full, 2 * 2 * sizeof(double));
  cudaMalloc((double**) &dev_TUpper_half, size);
  cudaMalloc((double**) &dev_TUpper_full, size);
  cudaMalloc((double**) &dev_TLower_half, size);
  cudaMalloc((double**) &dev_TLower_full, size);
  cudaMalloc((double**) &dev_VPropogator, size);

  hostToDevice(&TDiag[FULL], dev_TDiag_full, 1);
  hostToDevice(TUpper[FULL], dev_TUpper_full, MESH_SIZE);
  hostToDevice(TLower[FULL], dev_TLower_full, MESH_SIZE);

  hostToDevice(&TDiag[HALF], dev_TDiag_half, 1);
  hostToDevice(TUpper[HALF], dev_TUpper_half, MESH_SIZE);
  hostToDevice(TLower[HALF], dev_TLower_half, MESH_SIZE);
}

/* Cleans up cudaMalloced memory */
void GPU_Finalize() {
  cudaFree(&dev_psi);
  cudaFree(&dev_workPsi);
  cudaFree(&dev_TDiag_half);
  cudaFree(&dev_TDiag_full);
  cudaFree(&dev_TUpper_half);
  cudaFree(&dev_TUpper_full);
  cudaFree(&dev_TLower_half);
  cudaFree(&dev_TLower_full);
  cudaFree(&dev_VPropogator);
}

/* Initialize parameters for the simulation by reading in from a file **********************/
void init_Parameters() {
  scanf("%le", &Lx);
  scanf("%le", &deltaTime);
  scanf("%d", &nStep);
  scanf("%d", &energyInterval);
  scanf("%le %le %le", &x0, &s0, &e0);
  scanf("%le %le", &barrierHeight, &barrierWidth);
  scanf("%le", &edgePotential);

  if (DEBUG) {
    printf("Lx = %le\n", Lx);
    printf("deltaT = %le\n", deltaTime);
    printf("nStep = %d\n", nStep);
    printf("energyInterval = %d\n", energyInterval);
    printf("x0, s0, e0 = %le %le %le\n", x0, s0, e0);
    printf("bH, BW = %le %le\n", barrierHeight, barrierWidth);
    printf("edgePotential = %le\n", edgePotential);
  }

  MPI_Bcast(&Lx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&deltaTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nStep, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&energyInterval, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&x0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&s0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&e0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&barrierWidth, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&barrierHeight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&edgePotential, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // Mesh Size
  dx = (double) Lx / NX;
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
    x = (dx*i + Lx*myid);

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
  double DXSX, LXMYID;          //for debugging

  // Calculuate Psi point-by-point
  for (sx=1; sx<=NX; ++sx) {
    DXSX = (double) dx * sx;
    LXMYID = (double) Lx * myid;
    temp = DXSX + LXMYID;

    x = temp - x0;
    gaussian = exp((-0.25 * x * x) / (s0 * s0));
    host_psi[sx][RE] = gaussian * cos(sqrt(2.0 * e0) * x);
    host_psi[sx][IM] = gaussian * sin(sqrt(2.0 * e0) * x);
  }

  // Normalize
  temp = 0.0;
  for (sx=1; sx<=NX; ++sx) {
    for (c=0; c<2; ++c) {
      temp += (host_psi[sx][c]*host_psi[sx][c]);
    }
  }
  MPI_Allreduce(&temp, &psiSquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  psiSquared *= dx;
  //printf("<Proc: %i> PSI SQUARED: %le\n", myid, psiSquared);

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
  potential_Propogation();          //Half Potential Propogator

  kinetic_Propogation(HALF);        //Half Kinetic Propogator
  kinetic_Propogation(FULL);        //Full Kinetic Propogator
  kinetic_Propogation(HALF);        //Half Kinetic Propogator

  potential_Propogation();

  //getPsiSquared();
}

void getPsiSquared() {
  double temp, psiSquared;
  int sx, c;

  temp = 0.0;
  psiSquared = 0.0;
  for (sx=1; sx<=NX; ++sx) {
    for (c=0; c<2; ++c) {
      temp += (host_psi[sx][c]*host_psi[sx][c]);
    }
  }
  MPI_Allreduce(&temp, &psiSquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  psiSquared *= dx;

  if (myid==0) printf("Psi Squared: %le\n", psiSquared);
}

/* Half Potential Energy propogator Function */
void potential_Propogation() {
  if (usingCUDA) {
    hostToDevice(host_psi, dev_psi, MESH_SIZE);
    gpu_Potential_Prop <<<gridDim, blockDim>>> (dev_psi, dev_VPropogator);
    deviceToHost(dev_psi, host_psi, MESH_SIZE);
  }
  else {
    regular_Potential_Prop();
  }
}

__global__ void gpu_Potential_Prop(double *psi, double *vProp) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int sx, s_Re, s_Im;
  double workPsi_Re, workPsi_Im;
  sx= tid+1;
  s_Re = 2*sx;
  s_Im = 2*sx+1;

  workPsi_Re = vProp[s_Re]*psi[s_Re] - vProp[s_Im]*psi[s_Im];
  workPsi_Im = vProp[s_Re]*psi[s_Im] + vProp[s_Im]*psi[s_Re];
  psi[s_Re] = workPsi_Re;
  psi[s_Im] = workPsi_Im;
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

/* Half|Full Kinetic Energy Propogator Function */
void kinetic_Propogation(int stepSize) {
  //First apply periodic bounds
  periodic_Bounds();

  if (usingCUDA) {
    hostToDevice(host_psi, dev_psi, MESH_SIZE);

    if (stepSize == FULL) {
      gpu_Kinetic_Prop<<<gridDim, blockDim>>>(dev_psi, dev_workPsi, dev_TDiag_full, dev_TLower_full, dev_TUpper_full);
    }
    else {
      gpu_Kinetic_Prop<<<gridDim, blockDim>>>(dev_psi, dev_workPsi, dev_TDiag_half, dev_TLower_half, dev_TUpper_half);
    }

    gpu_WorkToPsi<<<gridDim, blockDim>>>(dev_workPsi, dev_psi);
    deviceToHost(dev_psi, host_psi, MESH_SIZE);
  }
  else {
    regular_Kinetic_Prop(stepSize);
  }
}

__global__ void gpu_Kinetic_Prop(double *psi, double *work, double *al, double *blx, double *bux) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int sx, s_Re, s_Im, l_Re, l_Im, u_Re, u_Im;
  double workPsi_Re, workPsi_Im;
  sx= tid+1;
  s_Re = 2*sx;
  s_Im = 2*sx+1;
  l_Re = 2*(sx-1);
  l_Im = 2*(sx-1)+1;
  u_Re = 2*(sx+1);
  u_Im = 2*(sx+1)+1;

  workPsi_Re = al[0]*psi[s_Re] - al[1]*psi[s_Im];
  workPsi_Im = al[0]*psi[s_Im] + al[1]*psi[s_Re];
  workPsi_Re += blx[s_Re]*psi[l_Re] - blx[s_Im]*psi[l_Im];
  workPsi_Im += blx[s_Re]*psi[l_Im] + blx[s_Im]*psi[l_Re];
  workPsi_Re += bux[s_Re]*psi[u_Re] - bux[s_Im]*psi[u_Im];
  workPsi_Im += bux[s_Re]*psi[u_Im] + bux[s_Im]*psi[u_Re];

  work[s_Re] = workPsi_Re;
  work[s_Im] = workPsi_Im;
}

void regular_Kinetic_Prop(int t) {
  int sx, c;
  double workPsi_Re, workPsi_Im;

  for (sx=1; sx<=NX; ++sx) {
    workPsi_Re = TDiag[t][RE]*host_psi[sx][RE] - TDiag[t][IM]*host_psi[sx][IM];
    workPsi_Im = TDiag[t][RE]*host_psi[sx][IM] + TDiag[t][IM]*host_psi[sx][RE];
    workPsi_Re += (TLower[t][sx][RE]*host_psi[sx-1][RE]) - (TLower[t][sx][IM]*host_psi[sx-1][IM]);
    workPsi_Im += (TLower[t][sx][RE]*host_psi[sx-1][IM]) + (TLower[t][sx][IM]*host_psi[sx-1][RE]);
    workPsi_Re += (TUpper[t][sx][RE]*host_psi[sx+1][RE]) - (TUpper[t][sx][IM]*host_psi[sx+1][IM]);
    workPsi_Im += (TUpper[t][sx][RE]*host_psi[sx+1][IM]) + (TUpper[t][sx][IM]*host_psi[sx+1][RE]);

    workPsi[sx][RE] = workPsi_Re;
    workPsi[sx][IM] = workPsi_Im;
  }

  for (sx=1; sx<=NX; ++sx) {
    for (c=0; c<2; ++c) {
      host_psi[sx][c] = workPsi[sx][c];
    }
  }
}

/*************************************************************************************
 * HOST-DEVICE CONVERSION Functions
 ************************************************************************************/

__global__ void gpu_WorkToPsi(double *work, double *psi) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int sx, s_Re, s_Im;

  sx= tid+1;
  s_Re = 2*sx;
  s_Im = 2*sx+1;
  psi[s_Re] = work[s_Re];
  psi[s_Im] = work[s_Im];
}

void hostToDevice(double** host, double* device, int size) {
  int i,j;
  double *hostBuf;
  hostBuf = (double*) malloc(sizeof(double) * size * 2);

  for (i=0; i<size; ++i) {
    for (j=0; j<2; ++j) {
      hostBuf[2*i + j] = host[i][j];
    }
  }

  cudaMemcpy((void*) device, hostBuf, 2*size*sizeof(double), cudaMemcpyHostToDevice);
  free(hostBuf);
}

void deviceToHost(double* device, double** host, int size) {
  int i,j;
  double *devBuf;
  devBuf = (double*) malloc(sizeof(double) * size * 2);

  cudaMemcpy((void*) devBuf, device, 2*size*sizeof(double), cudaMemcpyDeviceToHost);
  for (i=0; i<size; ++i) {
    for (j=0; j<2; ++j) {
      host[i][j] = devBuf[2*i + j];
    }
  }
  free(devBuf);
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
  MPI_Irecv(recvBuf, 2, MPI_DOUBLE, neighbor[left], DAEMON, MPI_COMM_WORLD, &request);
  MPI_Send(sendBuf, 2, MPI_DOUBLE, neighbor[right], DAEMON, MPI_COMM_WORLD);
  MPI_Wait(&request, &status);
  host_psi[0][RE] = recvBuf[RE];
  host_psi[0][IM] = recvBuf[IM];

  // Send left, receive right
  sendBuf[0] = host_psi[1][RE];
  sendBuf[1] = host_psi[1][IM];
  MPI_Irecv(recvBuf, 2, MPI_DOUBLE, neighbor[right], DAEMON+1, MPI_COMM_WORLD, &request);
  MPI_Send(sendBuf, 2, MPI_DOUBLE, neighbor[left], DAEMON+1, MPI_COMM_WORLD);
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
  expected[_E_] = expected[_T_] + expected[_V_];
}
