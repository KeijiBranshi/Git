#include "matrix.h"
#include <iostream>

int main() {
  std::cout << "Starting main.cpp" << std::endl;

  double array[3][3] =
  {
    {5, 8, 1},
    {3, 9, 2},
    {0, 6, 7}
  };
  Matrix mat(3,3);
  for (int i=0; i<3; ++i){
    for (int j=0; j<3; ++j){
      mat[i][j] = array[i][j];
    }
  }

  Matrix sum = mat + mat;

  mat.print();
  sum.print();
}
