/************************************************************************************
 * 1D Quantum Dynamics using MPI and CUDA
 * Built off of/modified from from Dr. Nakano's qd1.c code
 ************************************************************************************/
#include <cuda.h>
#include <mpi.h>

#define NX 2496
#define MESH_SIZE 2498
#define PIE 3.141592653589793   /* extra 'E' is intentional */
#define NUM_BLOCKS 13
#define NUM_THREADS 192
#define DAEMON 666

#define RE 0
#define IM 1
#define COMPLEX 2
#define HALF 0
#define FULL 1

#define _T_ 0
#define _V_ 1
#define _E_ 2

#define right 0
#define left 1

#define MAX_LINE 256
/* MPI/GPU Parameters/Constants ******************************************************/
int nProc;
int myid;
MPI_Status status;
MPI_Request request;

bool usingCUDA;
const dim3 gridDim(NUM_BLOCKS, 1, 1);     //Grid dimensions
const dim3 blockDim(NUM_THREADS, 1, 1);   //Block dimensions

/* Input Parameters ******************************************************************/
//int NX;
//int MESH_SIZE;

double Lx;          /* Simulation box lengths in the x direction */
double deltaTime;   /* Time discretization unit */
int nStep;          /* Total # of simulation steps */
int energyInterval;        /* Interval to calculate energies */
double x0,s0,e0;     /* Center-of-mass, spread, & energy of initial wave function */
double barrierHeight, barrierWidth;
double edgePotential;

double dx;
double expected[3];

/* Wave Function and Propogators (CPU) *************************************************/
double **host_psi;                        /*[Mesh point][Re|Im]               */
double **workPsi;                         /*[Mesh point][Re|Im]               */
double **TDiag;                           /*[HALF|FULL][Re|Im]                */
double ***TUpper;                         /*[HALF|FULL][Mesh point][Re|Im]    */
double ***TLower;                         /*[HALF|FULL][Mesh point][Re|Im]    */
double **VPropogator;                     /*[Mesh point][Re|Im]               */
double *potential;

/* Wave Function and Propogators (GPU) *************************************************/
double *dev_psi;
double *dev_workPsi;

double *dev_TDiag_half;
double *dev_TDiag_full;
double *dev_TUpper_half;
double *dev_TUpper_full;
double *dev_TLower_half;
double *dev_TLower_full;

double *dev_VPropogator;

/* Simulation Functions **************************************************************/
void init_Parameters();
void init_Variables();
void init_Propogators();
void init_WaveFunction();
void init_GPU();
void GPU_Finalize();
void cleanUpVariables();

void single_Step();
void potential_Propogation();
void kinetic_Propogation(int);
void regular_Kinetic_Prop(int);
void regular_Potential_Prop();

void getPsiSquared();
void periodic_Bounds();
void calculate_Energy();

/* CUDA Related Functions ***********************************************************/
__global__ void gpu_Potential_Prop(double*,double*);
__global__ void gpu_Kinetic_Prop(double*,double*,double*,double*,double*);
__global__ void gpu_WorkToPsi(double*,double*);
void hostToDevice(double**,double*,int);
void deviceToHost(double*,double**,int);
