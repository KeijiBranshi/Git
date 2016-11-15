/************************************************************************************
 * 1D Quantum Dynamics using MPI and CUDA
 * Built off of/modified from from Dr. Nakano's qd1.c code
 ************************************************************************************/
#include <mpi.h>

#define NX 2496    /* Mesh points in x direction */
#define PIE 3.141592653589793   /* extra 'E' is intentional */
#define NUM_BLOCKS 13
#define NUM_THREADS 192
#define DAEMON 666

#define RE 0
#define IM 1
#define COMPLEX 2
#define HALF 0
#define FULL 1
#define MESH_SIZE 2498

#define _T_ 0
#define _V_ 1
#define _E_ 2

#define right 0
#define left 1

#define MAX_LINE 256
/* MPI Related Variables *************************************************************/
int nProc;
int myid;
MPI_Status status;
MPI_Request request;

/* Input Parameters ******************************************************************/
double Lx;          /* Simulation box lengths in the x & y directions */
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

/* Simulation Functions **************************************************************/
void init_Parameters();
void init_Variables();
void init_Propogators();
void init_WaveFunction();
void cleanUpVariables();

void single_Step();
void regular_Kinetic_Prop(int);
void regular_Potential_Prop();

void periodic_Bounds();
void calculate_Energy();
