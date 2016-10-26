#include <stdio.h>
#include <stdbool.h>
#include "mpi.h"

#define DAEMON 666
#define ROOT 0
#define MAX_DIM 512
#define MAX_GREY 255
#define MAX_LINE 1024

#define SQRT_2 1.4142136
#define NORMALIZED true

typedef FILE File;

/* System/Parallelism Related Variables *************************************/
int numberOfProcs;         // Total number of processors
int procsPerDirection[2];  // # of processors in the x/y direction
int myid;                  // Rank/ID for this processor
MPI_Status status;
MPI_Request request;

/* Daubechies' Constants ****************************************************/
double C0 = 0.4829629131;
double C1 = 0.8365163037;
double C2 = 0.2241438680;
double C3 = -0.1294095226;

/* File IO ******************************************************************/
void readImageLocal(File*,char*,double **image,int*,int*,int*); //reads in IMAGE locally (only root should do this)
void readImageParallel(double **image,int*,int*);               //send IMAGE read at root to other processors
void writeImageLocal(File*,char*,double **image,int,int);       //writes IMAGE locally (only root should do this)
void writeImageParallel(double **image,int,int);                //send a processor's portion of the image to root

/* Decomposition Functions **************************************************/
void exchangeCachedColumn(double **image,int,int);  // sends first 2 columns of IMAGE to horizontal neighbor
void exchangeCachedRow(double **image,int,int);     // sends first 2 rows of IMAGE to vertical neighbor
void compressCols(double **image,int,int);          // performs the avelet transform on columns
void compressRows(double **image,int,int);          // performs the avelet transform on rows
double getMax(double **image,int,int);              // gets max pixel value in an IMAGE (for scaling/normalization)

/* Neighbor ID Retrieval Functions **********************************************/
int getRightNeighborID() {
  if ((myid+1)%procsPerDirection[0]) return (myid + 1);
  else return (myid - (myid%procsPerDirection[0]));
}
int getLeftNeighborID() {
  if (!(myid%procsPerDirection[0])) return (myid + procsPerDirection[0] - 1);
  else return (myid - 1);
}
int getBottomNeighborID() {
  int neighborID = myid + procsPerDirection[0];
  if (neighborID < numberOfProcs) return neighborID;
  else return (neighborID%numberOfProcs);
}
int getTopNeighborID() {
  int neighborID = myid - procsPerDirection[0];
  if (neighborID >= ROOT) return neighborID;
  else return numberOfProcs + neighborID;
}
