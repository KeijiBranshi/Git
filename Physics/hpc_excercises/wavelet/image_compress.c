#include "image_compress.h"
#include "omp.h"
#include <math.h>
#include <stdlib.h>

void recursiveCompression(double **image, int *nRowPerProc, int *nColPerProc, int nRowTotal, int nColTotal, int level);

int main(int argc, char* argv[]) {
  double **image;
  int nRow,nCol,nGrey,r,c,i;
  int nRow_old, nCol_old, nRow_new, nCol_new;
  char line[MAX_LINE];
  int level;             // # of levels to recursively perform the wavelet transform
  File *f;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  image = malloc(MAX_DIM * sizeof(double));
  for (i=0; i<MAX_DIM; ++i) image[i] = malloc(MAX_DIM * sizeof(double));
  level = 3;
  procsPerDirection[0] = procsPerDirection[1] = 4;  //sqrt(numberOfProcs);    <== had to hardcode the number due to compiling issues on the HPC

  if (myid == ROOT) {
    readImageLocal(f, line, image, &nRow_old, &nCol_old, &nGrey);
    nRow = nRow_old;
    nCol = nCol_old;
  }

  recursiveCompression(image, &nRow, &nCol, nRow_old, nCol_old, level);

  if (myid == ROOT) {
    nRow_new = nRow_old / pow(2, level);
    nCol_new = nCol_old / pow(2, level);
    writeImageLocal(f, line, image, nRow_new, nCol_new);
  }

  for (i=0; i<MAX_DIM; ++i) {
    free(image[i]);
  }
  free(image);

  MPI_Finalize();

}

/*
 * Recursive function that calls all the necessary image compression and file sending/receiving
 * functions. Determines how many levels to perform the wavelet transform.
 */
void recursiveCompression(double **image, int *nRowPerProc, int *nColPerProc, int nRowTotal, int nColTotal, int level) {
  if (level == 0) return;
  if (myid == 0) {
    *nRowPerProc = (nRowTotal) / procsPerDirection[0];
    *nColPerProc = (nColTotal) / procsPerDirection[1];
  }
  readImageParallel(image, nRowPerProc, nColPerProc);
  exchangeCachedColumn(image, *nRowPerProc, *nColPerProc);
  exchangeCachedRow(image, *nRowPerProc, *nColPerProc);
  compressCols(image, *nRowPerProc, *nColPerProc);
  compressRows(image, *nRowPerProc, *nColPerProc);
  writeImageParallel(image, *nRowPerProc, *nColPerProc);

  recursiveCompression(image, nColPerProc, nColPerProc, nRowTotal/2, nColTotal/2, --level);
}

/*******************************************************************************
 * File/Image IO
 *******************************************************************************/

/*
 * Only the root processor call this function. This reads in the local file.
 */
void readImageLocal(File *f, char *line, double **image, int *nRow, int *nCol, int *nGrey) {
  int r,c;

  f=fopen("Lenna512x512.pgm","r");

  fgets(line,MAX_LINE,f);
  fgets(line,MAX_LINE,f);
  fgets(line,MAX_LINE,f);
  sscanf(line,"%d %d",nCol,nRow);
  fgets(line,MAX_LINE,f);
  sscanf(line,"%d",nGrey);

  for (r=0; r<(*nRow); ++r) {
    for (c=0; c<(*nCol); ++c) {
      image[r][c] = (double)fgetc(f);
    }
  }

  fclose(f);
}

/*
 * Called by each processor. The root will send the IMAGE out to all other processors.
 * All other processors will read in their portion of the IMAGE
 */
void readImageParallel(double **image, int *nRow, int *nCol) {
  int id, idRow, idCol, r=0, c=0, i=0;
  int totalPixels;
  double *buffer;

  // Broadcast row and column numbers to all pocessors
  MPI_Bcast(nRow, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(nCol, 1, MPI_INT, 0, MPI_COMM_WORLD);

  totalPixels = (*nRow) * (*nCol);
  buffer = malloc(totalPixels * sizeof(double));

  // The root will send each processor it's portion of the IMAGE to handle
  if (myid == ROOT) {
    for (id=1; id<numberOfProcs; ++id) {
      idRow = id / procsPerDirection[1];
      idCol = id % procsPerDirection[1];

      i=0;
      for (r=0; r<(*nRow); ++r) {
        for (c=0; c<(*nCol); ++c) {
          buffer[i++] = image[(*nRow)*idRow + r][(*nCol)*idCol + c];
        }
      }

      MPI_Send(buffer, totalPixels, MPI_DOUBLE, id, DAEMON, MPI_COMM_WORLD);
    }
  }
  else {
    // Each processor other than the root reads in the image
    MPI_Recv(buffer, totalPixels, MPI_DOUBLE, ROOT, DAEMON, MPI_COMM_WORLD, &status);

    i=0;
    for (r=0; r<(*nRow); ++r) {
      for (c=0; c<(*nCol); ++c) {
        image[r][c] = buffer[i++];
      }
    }
  }

  free(buffer);
}

/*
 * Writes the final IMAGE to a file. Only the root calls this function
 */
void writeImageLocal(File *f, char *line, double **image, int nRow, int nCol) {
  int r,c,work;
  double max;

  max = getMax(image, nRow, nCol);  //gets max value for normalization

  f=fopen("Lenna64x64.pgm","w");

  fprintf(f,"P5\n");
  fprintf(f,"# Simple image test\n");
  fprintf(f,"%d %d\n",nCol,nRow);
  fprintf(f,"%d\n",MAX_GREY);

  for (r=0; r < nRow; r++) {
    for (c=0; c < nCol; c++) {
      work = (int) ( (image[r][c] / max) * MAX_GREY);
      fputc((char)work,f);
    }
  }

  fclose(f);
}

/*
 * Called by all processes. Everyone send their compressed IMAGE data to the root.
 */
void writeImageParallel(double **image, int nRow, int nCol) {
  double *buffer;
  int r,c,id,i=0;
  int idRow, idCol;
  buffer = malloc ((nRow/2) * (nCol/2) * sizeof(double));

  if (myid > ROOT) {
    //Store image in buffer in row-major order
    i = 0;
    for (r=0; r<nRow/2; ++r) {
      for (c=0; c<nCol/2; ++c) {
        buffer[i++] = image[r][c];
      }
    }

    //Send my compressed IMAGE in BUFFER to the to root node
    MPI_Send(buffer, (nCol/2)*(nRow/2), MPI_DOUBLE, ROOT, DAEMON, MPI_COMM_WORLD);
  }
  else {
    //piece together the image from other processors
    for (id=1; id<numberOfProcs; ++id) {
      idRow = id/procsPerDirection[0];
      idCol = id%procsPerDirection[0];
      //printf("Attempting to receive from %i\n", id);
      MPI_Recv(buffer, (nCol/2)*(nRow/2), MPI_DOUBLE, id, DAEMON, MPI_COMM_WORLD, &status);
      //printf("Receiving image from %i\n", id);

      i = 0;
      for (r=0; r<nRow/2; ++r) {
        for (c=0; c<nCol/2; ++c) {
          image[idRow*(nRow/2)+r][idCol*(nCol/2)+c] = buffer[i++];
        }
      }
    }
  }

  free(buffer);
}

/*******************************************************************************
 * Decomposition Functions
 *******************************************************************************/

/*
 * Sends the current processes' first two columns in IMAGE to its left
 * neighbor process. Then it receives two columns from its right neighbor process.
 */
void exchangeCachedColumn(double **image, int nRow, int nCol) {
  double *sendBuf, *recvBuf;
  int r,c,i=0;

  sendBuf = malloc( (nRow * 2) * sizeof(double));
  recvBuf = malloc( (nRow * 2) * sizeof(double));

  //Copy my first 2 columns into the SENDBUF
  for (r=0; r<nRow; ++r) {
    for (c=0; c<2; ++c) {
      sendBuf[i++] = image[r][c];
    }
  }

  //Send my first 2 columns to my left niehgbor processor, wait for columns from my right neighbor
  MPI_Irecv(recvBuf, (nRow * 2), MPI_DOUBLE, MPI_ANY_SOURCE, DAEMON, MPI_COMM_WORLD, &request);
  MPI_Send(sendBuf, (nRow * 2), MPI_DOUBLE, getLeftNeighborID(), DAEMON, MPI_COMM_WORLD);
  MPI_Wait(&request, &status);

  //Append message to my image
  i=0;
  for (r=0; r<nRow; ++r) {
    for (c=0; c<2; ++c) {
      image[r][nCol+c] = recvBuf[i++];
    }
  }

  free(sendBuf);
  free(recvBuf);
}

/*
 * Sends the first two rows of the current process to its top neighbor process.
 * Then it receives two rows from its lower neighbor process.
 */
void exchangeCachedRow(double **image, int nRow, int nCol) {
  double *sendBuf, *recvBuf;
  int r,c,i=0;

  sendBuf = malloc( (nCol * 2) * sizeof(double));
  recvBuf = malloc( (nCol * 2) * sizeof(double));

  //copy my first 2 rows into SENDBUF
  for (c=0; c<nCol; ++c) {
    for (r=0; r<2; ++r) {
      sendBuf[i++] = image[r][c];
    }
  }

  //Send my first 2 rows to my top neighbor processor, wait for rows from my bottom neighbor
  MPI_Irecv(recvBuf, (nCol * 2), MPI_DOUBLE, MPI_ANY_SOURCE, DAEMON, MPI_COMM_WORLD, &request);
  MPI_Send(sendBuf, (nCol * 2), MPI_DOUBLE, getTopNeighborID(), DAEMON, MPI_COMM_WORLD);
  MPI_Wait(&request, &status);

  //Append message to my image
  i=0;
  for (c=0; c<nCol; ++c) {
    for (r=0; r<2; ++r) {
      image[nRow+r][c] = recvBuf[i++];
    }
  }

  free(sendBuf);
  free(recvBuf);
}

/*
 * Compresses NROW rows of IMAGE down to NROW/2 rows.
 * Compresses according to Daubechies method using OpenMP.
 * Only keeps the smooth portion of the IMAGE (detail is not computed).
 */
void compressRows(double **image, int nRow, int nCol) {
  #pragma omp parallel num_threads(4)
  {
    int r,c;

    #pragma omp for private(r)
    for (c=0; c<nCol; ++c) {
      for (r=0; r<nRow/2; ++r) {
        image[r][c] = C0 * image[2*r][c] + C1 * image[2*r+1][c]
                          + C2 * image[2*r+2][c] + C3 * image[2*r+3][c];
      }
    }
  }
}

/*
 * Compresses NCOL columns of IMAGE down to NCOL/2 columns.
 * Compresses according to Daubechies method using OpenMP.
 * Only keeps the smooth portion of the IMAGE (detail is not computed).
 */
void compressCols(double **image, int nRow, int nCol) {
  #pragma omp parallel num_threads(4)
  {
    int r,c;

    #pragma omp for private(c)
    for (r=0; r<nRow; ++r) {
      for (c=0; c<nCol/2; ++c) {
        image[r][c] = C0 * image[r][2*c] + C1 * image[r][2*c+1]
                          + C2 * image[r][2*c+2] + C3 * image[r][2*c+3];
      }
    }
  }
}

/*
 * Gets the maximum pixel value in the IMAGE from row 0:nRow and column 0:nCol.
 */
double getMax(double **image, int nRow, int nCol) {
  int r, c;
  double max = 0.0;

  for (r=0; r<nRow; ++r) {
    for (c=0; c<nCol; ++c) {
      if (image[r][c] > max) max = image[r][c];
    }
  }

  return max;
}
