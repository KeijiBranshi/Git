#include "matrix.h"
#include <iostream>
#include <string>
#include <sstream>

Matrix::Matrix(int row, int col)
{
  this->nRow = row;
  this->nCol = col;
  this->mat = new double*[row];
  for (int i=0; i<row; ++i){
    this->mat[i] = new double[col];
    for (int j=0; j<col; ++j){
      this->mat[i][j] = 0.0;          //initialize all elements to 0
    }
  }
}

Matrix::~Matrix()
{
  this->clear();
}

void Matrix::clear()
{
  for (int i=0; i<this->nRow; ++i){
    delete this->mat[i];
  }
  delete this->mat;
}

bool Matrix::set(int row, int col, double val)
{
  this->mat[row][col] = val;
}

Matrix& Matrix::operator=(const Matrix& other)
{
  this->clear();
  for (int i=0; i<this->nRow; ++i){
    for (int j=0; j<this->nCol; ++i){
      this->set(i,j,other->[i][j]);
    }
  }
}

Matrix Matrix::getZero(int row, int col)
{
  Matrix zero = Matrix(row,col);
  for (int i=0; i<row; ++i){
    for (int j=0; j<col; ++i){
      this->zero.set(i,j,0);
    }
  }

  return zero;
}
Matrix Matrix::getIdentity(int row, int col)
{
  Matrix identity = Matrix(row,col);
  for (int i=0; i<row; ++i){
    for (int j=0; j<col; ++i){
      if (i==j) this->identity.set(i,j,1);
      else this->identity.set(i,j,0);
    }
  }

  return identity;
}

string Matrix::print()
{
  std::ostringstream oss;
  for (int i=0; i<this->nRow; ++i) {
    for (int j=0; j<this->nCol; ++j) {
      if (j=0) oss << "|";              //left bracket
      oss << this->mat[i][j];           //element i,j of matrix
      if (j=this->nCol-1) oss << "|";   //right bracket
      else oss << " " ;                 //space within matrix
    }
  }

  std::cout << oss.str();
  return oss.str();
}
