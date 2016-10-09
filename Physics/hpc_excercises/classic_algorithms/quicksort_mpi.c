/******************************************************************************
 * Hypercube Quicksort using MPI
 * Format/order of MPI calls taken from Dr. Aiichiro Nakano's parallel hybercube
 * mergesort algorithm
 ******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

#define N 1024 /* Max list size */
#define MAX 99 /* Max value of list element */
#define MAXP 32
#define MAXD 5

int nprocs, dim, myid; /* Cube size, dimension, & node ID */

void swap (int *list, int i, int j) {
  int temp = list[i];
  list[i] = list[j];
  list[j] = temp;
}

int avg(int *list, int n) {
  int i, sum=0;
  for (i=0; i<n; ++i) {
    sum += list[i];
  }
  return sum / n;
}

int partition(int *list, int left, int right, int* pval, int* pindex) {
  int i,j,pivot;

  if (left < right) {
    i = left-1; j = right+1;

    if (pval) {
      i = left-1;
      pivot = *pval;
    }
    else {
      i = left;
      pivot = list[left];
    }

    do {
      while (list[++i]<pivot && i<=right);
      while (list[--j]>pivot && j>left-1);
      if (i<j) {
        swap(list,i,j);
      }
    } while (i<j);

    if (!pval) swap(list,left,j);
  }
  else
    j = right;

  *pindex = j;
}

void sequential_quicksort(int *list, int l, int h) {
 int pivot;
 if (l < h) {
   partition(list, l, h, NULL, &pivot);
   sequential_quicksort(list, l, pivot-1);
   sequential_quicksort(list, pivot+1, h);
 }
}

void parallel_quicksort(int *list, int *n) {
  int L, k, pivot, pivot_idx;
  int mask=pow(2,dim)-1, bitvalue=pow(2,dim-1), partner=myid^bitvalue;	//mask=111, bitvalue=100
  int nelem=*n, nsend, nrecv;
  int list_buff[2*N], nless;

  int nprocs_cube, ranks[nprocs], c=0;
  MPI_Status status;
  MPI_Comm cube[MAXD][MAXP];                    //2D array of communicators
  MPI_Group cube_group[MAXD][MAXP];             //2D array of groups associated with above communicators
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs_cube);
  cube[dim][0] = MPI_COMM_WORLD;                //set initial Commmunicator
 
  //bitvalue = nprocs >> 1;
  //mask = nprocs-1;

  for (L=dim; L>0; --L) { 
    if ((myid & mask) == 0) {
      pivot = avg(list, nelem);
    }
    MPI_Bcast(&pivot, 1, MPI_INT, 0, cube[L][myid/nprocs_cube]);  //Broadcast to nodes in the same group
   
    //mask = mask ^ bitvalue;
    //bitvalue = bitvalue >> 1;
    partition(list, 0, nelem-1, &pivot, &pivot_idx);
    nless = pivot_idx+1;
 
    partner = myid ^ bitvalue;
    if ((myid & bitvalue) == 0) {
      nsend = nelem - nless;
      MPI_Send(&nsend, 1, MPI_INT, partner, 666, MPI_COMM_WORLD);//Send #  of elements to pass to partner
      MPI_Recv(&nrecv, 1, MPI_INT, partner, 666, MPI_COMM_WORLD, &status);//Recv # of elements from partner
      MPI_Send((list+nless), nsend, MPI_INT, partner, 666, MPI_COMM_WORLD);//send right sublist
      MPI_Recv(list_buff, nrecv, MPI_INT, partner, 666, MPI_COMM_WORLD, &status);//receive partner's left sublist

      /*Append left received list to my left list*/
      memcpy(&list[nless], list_buff, 4*4*4*4*nrecv);
    }
    else {
      nsend = nless;
      MPI_Send(&nsend, 1, MPI_INT, partner, 666, MPI_COMM_WORLD);
      MPI_Recv(&nrecv, 1, MPI_INT, partner, 666, MPI_COMM_WORLD, &status);
      MPI_Send(list, nsend, MPI_INT, partner, 666, MPI_COMM_WORLD); //send left sublist to partner
      MPI_Recv(list_buff, nrecv, MPI_INT, partner, 666, MPI_COMM_WORLD, &status); //receive right sublist of partner;

      /* Append the received list to my right sublist*/
      memcpy(&list_buff[nrecv], &list[nsend], 4*(nelem-nsend));
      memcpy(list, list_buff, 4*(nelem-nsend+nrecv));
    }
    nelem = nelem - nsend + nrecv;
    memset(list_buff, 0, N);
    mask = mask ^ bitvalue;
    bitvalue >> 1;
   
    nprocs_cube = pow(2,L);
    c = myid/nprocs_cube;

    /* Creates to child communicators from the current one */
    MPI_Comm_group(cube[L][c], &(cube_group[L][c]));

    nprocs_cube = nprocs_cube/2;
    for (k=0; k<nprocs_cube; ++k) ranks[k] = k;
    MPI_Group_incl(cube_group[L][c], nprocs_cube, ranks, &(cube_group[L-1][2*c  ]));
    MPI_Group_excl(cube_group[L][c], nprocs_cube, ranks, &(cube_group[L-1][2*c+1]));
    MPI_Comm_create(cube[L][c], cube_group[L-1][2*c  ], &(cube[L-1][2*c  ]));
    MPI_Comm_create(cube[L][c], cube_group[L-1][2*c+1], &(cube[L-1][2*c+1]));

    MPI_Group_free(&(cube_group[L  ][c    ]));
    MPI_Group_free(&(cube_group[L-1][2*c  ]));
    MPI_Group_free(&(cube_group[L-1][2*c+1]));
    /*******************************************************/ 
  }
  sequential_quicksort(list, 0, nelem-1);
  *n = nelem;
}

int main(int argc, char* argv[]) {
  int list[N], n=4;
  int i;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  dim = log (nprocs + 1e-10) / log (2.0);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  srand((unsigned) myid+1);
  for (i=0; i<n; ++i) list[i] = (rand()) % MAX;

  printf("Before Sort: Node %2d :", myid);
  for (i=0; i<n; ++i) printf("%3d ", list[i]);
  printf("\n");

  //sequential_quicksort(list, 0, n-1);
  parallel_quicksort(list, &n);

  printf("After Sort: Node %2d :", myid);
  for (i=0; i<n; ++i) printf("%3d ", list[i]);
  printf("\n");

  MPI_Finalize();
  return 0;
}
