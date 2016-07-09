/* If there's anything I learned in undergrad physics, it's that matrices are god. Here's my attempt at understanding them better*/
#ifndef MATRIX_H
#define MATRIX_H

class Matrix {
private:
  /*Variable*/
  int nRows;
  int nCols;
  double** mat;

  /*Functions*/
  void clear();

public:
  /*Functions*/
  Matrix(int,int);
  Matrix(int,int,double[][])
  ~Matrix();

  bool setAll(double**);
  bool setRow(int,double*);
  bool setColumn(int,double*);
  bool set(int,int,double);

  Matrix operator+(const Matrix&);
  Matrix operator-(const Matrix&);
  Matrix operator*(const Matrix&);
  double* operator[](const Matrix&);
  Matrix& operator=(const Matrix&);
  bool operator==(const Matrix&);

  Matrix getZero();
  Matrix getIdentity();
  double* eigenvalues();
  double** eigenvectors();

  string print();
};

#endif
