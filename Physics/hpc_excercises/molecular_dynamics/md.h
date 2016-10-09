/* Molecular Dynamics */
#define NMAX 10000 //Max # of atoms to be handled by the program

/* Variables ****************************************************/
int nAtom;                      // # of atoms to actually be run
int stepCount;
double r[NMAX][3];              // Position Vectors
double v[NMAX][3];              // Velocity Vectors
double a[NMAX][3];              // Acceleration Vectors

/* Constants ****************************************************/
double temperature;
double deltaT;

/* Initialization Functions *************************************/
void initParam(void);
void initSystem(void);

/* Step Calculation Functions ***********************************/
void initAccel(void);               // computes a(r), based on r(0)
void halfStepVel(void);             // computes v(t + delta/2)
void singleStepPos(void);           // computes r(t + delta)
void computeAccel(void);            // computes a(t + delta)
void singleStepVel(void);           // computes v(t + delta)
