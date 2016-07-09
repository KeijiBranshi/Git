/* If there's anything I learned in undergrad physics, it's that matrices are god. Here's my attempt at understanding them better*/
#ifndef MATRIX_H
#define MATRIX_H
#include <string>

class Matrix {
private:
  /*Variable*/
  int nRow;
  int nCol;
  double** mat;

  /*Functions*/
  void clear();
  void reset();

public:
  /*Functions*/
  Matrix(int,int);
  template<typename T> Matrix(int,int,T&);
  Matrix(const Matrix&);
  ~Matrix();

  //bool setAll(double**);
  //bool setRow(int,double*);
  //bool setColumn(int,double*);
  bool set(int,int,double);

  double get(int,int) const;
  int rows() const;
  int columns() const;

  Matrix operator+(const Matrix&);
  //Matrix operator-(const Matrix&);
  //Matrix operator*(const Matrix&);
  Matrix& operator=(const Matrix&);
  //bool operator==(const Matrix&);
  double* operator[](int);

  Matrix getZero();
  Matrix getIdentity();
  Matrix getZero(int,int);
  Matrix getIdentity(int,int);
  //double* eigenvalues();
  //double** eigenvectors();

  std::string print() const;
  std::string details() const;
};

#endif
