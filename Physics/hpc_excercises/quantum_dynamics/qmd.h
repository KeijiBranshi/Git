/************************************************************************************
 * 1D Quantum Dynamics using MPI and CUDA
 * Built off of/modified from from Dr. Nakano's qd1.c code
 ************************************************************************************/

#define NX 2496    /* Mesh points in x direction */
#define PIE 3.141592653589793   /* extra 'E' is intentional */

#define RE 0
#define IE 1
#define HALF 0
#define FULL 1

#define T 0
#define V 1
#define E 2

#define right 0
#define left 1

/* MPI Related Variables *************************************************************/
int nProc;
int myid;
MPI_Status status;
MPI_Request request;

/* Simulation Functions **************************************************************/
void init_Parameters();
void init_Propogators();
void init_WaveFunction();

void single_Step();
void potential_Propogation();
void kinetic_Propogation(int);

void periodic_Bounds();
void calculate_Energy();

/* Input Parameters ******************************************************************/
double Lx;          /* Simulation box lengths in the x & y directions */
double deltaTime;   /* Time discretization unit */
int nStep;          /* Total # of simulation steps */
int nEnergy;        /* Interval to calculate energies */
double x0,s0,e0     /* Center-of-mass, spread, & energy of initial wave function */
double barrierHeight, barrierWidth;
double edgePotential;

/* Wave Function and Propogators *****************************************************/
double psi[NX+2][2];
double tempPsi[NX+2][2];
double TDiag[2][2];
double TUpper[2][NX+2][2];
double TLower[2][NX+2][2];
double VPropogator[NX+2][2];
double V[NX+2];

/* Other Variables *******************************************************************/
double dx;
double expected[3];
