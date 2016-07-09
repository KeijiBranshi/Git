#include "matrix.h"
#include <iostream>

int main() {
  Matrix* mat = new Matrix(2,2);

  for (int i=0; i<row; ++i){
    for (int j=0; j<col; ++i){
      if (i==j) mat.set(i,j,1);
      else mat.set(i,j,0);
    }
  }

  mat.print();
}
