#include "matrix.h"
#include <iostream>

int main() {
  std::cout << "Starting main.cpp" << std::endl;

  double array1[3][3] =
  {
    {5, 8, 1},
    {3, 9, 2},
    {0, 6, 7}
  };
  double array2[3][4] =
  {
    {5, 8, 1, 7},
    {3, 9, 2, 1},
    {0, 6, 7, 1}
  };

  Matrix mat1(3,3);
  Matrix mat2(3,4);
  for (int i=0; i<3; ++i){
    for (int j=0; j<3; ++j){
      mat1[i][j] = array1[i][j];
      mat2[i][j] = array2[i][j];
      if (j==2) mat2[i][j+1] = array2[i][j+1];
    }
  }

  Matrix sum = mat1 + mat2;

  mat1.print();
  std::cout << " + \n";
  mat2.print();
  std::cout << " != \n";
  sum.print();
}
