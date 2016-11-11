#define NDIM 2    /* Spacial dimensions of the system */
#define NX 2496    /* Mesh points in x direction */
#define PIE 3.141592653589793   /* extra 'E' is intentional */
#define RE 0
#define IE 1
#define HALF 0
#define FULL 1

/* MPI Related Variables *************************************************************/
int nProc;
MPI_Status status;
MPI_Request request;

/* Simulation Functions **************************************************************/
void init_Parameters();
void init_Propogation();
void init_WaveFunction();

void single_Step();
void potential_Prop_Step();
void kinetic_Prop_Step(int,int);

void periodic_Bounds();
void calculate_Energy();

/* Input Parameters ******************************************************************/
double Lx, Ly;      /* Simulation box lengths in the x & y directions */
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
double totalT;
double totalV;
double totalE;
