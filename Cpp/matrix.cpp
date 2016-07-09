#include "matrix.h"
#include <iostream>
#include <sstream>

/*Constructors*/
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

template<typename T>
Matrix::Matrix(int row, int col, T& original)
{
  this->nRow = row;
  this->nCol = col;
  this->mat = new double*[this->nRow];
  for (int i=0; i<this->nRow; ++i){
    this->mat[i] = new double[this->nCol];
    for (int j=0; j<this->nCol; ++j){
      this->mat[i][j] = original[i][j];
    }
  }
}

Matrix::Matrix(const Matrix& original)
{
  this->nRow = original.rows();
  this->nCol = original.columns();
  this->mat = new double*[this->nRow];
  for (int i=0; i<this->nRow; ++i){
    this->mat[i] = new double[this->nCol];
    for (int j=0; j<this->nCol; ++j){
      this->mat[i][j] = original.get(i,j);
    }
  }
}

Matrix::~Matrix()
{
  this->clear();
}

/*Private Functions*/
void Matrix::clear()
{
  for (int i=0; i<this->nRow; ++i){
    if (this->mat[i] != NULL)
      delete this->mat[i];
  }
  if (this->mat != NULL)
    delete this->mat;
}

void Matrix::reset()
{
  this->clear();

  /*re-instantiate all values to 0*/
  this->mat = new double*[this->nRow];
  for (int i=0; i<this->nRow; ++i){
    this->mat[i] = new double[this->nCol];
    for (int j=0; j<this->nCol; ++j){
      this->mat[i][j] = 0.0;
    }
  }
}

/*Public Functions*/
bool Matrix::set(int row, int col, double val)
{
  if (row<this->nRow && col<this->nCol && row>=0 && col>=0)
  {
    this->mat[row][col] = val;
    return true;
  }
  else
    return false;
}

double Matrix::get(int row, int col) const
{
  return this->mat[row][col];
}

int Matrix::rows() const
{
  return this->nRow;
}

int Matrix::columns() const
{
  return this->nCol;
}

Matrix& Matrix::operator=(const Matrix& other)
{
  this->reset();
  for (int i=0; i<this->nRow; ++i){
    for (int j=0; j<this->nCol; ++i){
      this->set(i,j,other.get(i,j));
    }
  }

  return *this;
}

Matrix Matrix::operator+(const Matrix& other)
{
  if (this->nRow != other.rows() && this->nCol != other.columns()){
    std::cout << "Invalid Matrix Addition:\n" << "Attempting to add a " <<
    this->details() << "with a " << other.details();

    return *this;
  }
  else {
    Matrix sum(this->nRow, this->nCol);
    for (int i=0; i<this->nRow; ++i){
      for (int j=0; j<this->nCol; ++j){
        sum[i][j] = this->get(i,j) + other.get(i,j);
      }
    }

    return sum;
  }

}

double* Matrix::operator[](int index)
{
  return this->mat[index];
}

Matrix Matrix::getZero()
{
  Matrix zero(this->nRow,this->nCol);

  for (int i=0; i<this->nRow; ++i){
    for (int j=0; j<this->nCol; ++i){
      zero.set(i,j,0);
    }
  }

  return zero;
}
Matrix Matrix::getIdentity()
{
  Matrix identity(this->nRow,this->nCol);

  for (int i=0; i<this->nRow; ++i){
    for (int j=0; j<this->nCol; ++j){
      if (i==j) identity.set(i,j,1);
      else identity.set(i,j,0);
    }
  }

  return identity;
}

Matrix Matrix::getZero(int row, int col)
{
  Matrix zero(row,col);
  for (int i=0; i<row; ++i){
    for (int j=0; j<col; ++i){
      zero.set(i,j,0);
    }
  }
  return zero;
}
Matrix Matrix::getIdentity(int row, int col)
{
  Matrix identity(row,col);
  for (int i=0; i<row; ++i){
    for (int j=0; j<col; ++i){
      if (i==j) identity.set(i,j,1);
      else identity.set(i,j,0);
    }
  }
  return identity;
}

std::string Matrix::print() const
{
  std::ostringstream oss;
  for (int i=0; i<this->nRow; ++i) {
    for (int j=0; j<this->nCol; ++j) {
      if (j==0) oss << "|";                 //left bracket
      oss << this->mat[i][j];               //element i,j of matrix
      if (j==this->nCol-1) oss << "|\n";      //right bracket
      else oss << " " ;                     //space within matrix
    }
  }

  std::cout << oss.str();
  return oss.str();
}

std::string Matrix::details() const
{
  std::ostringstream oss;
  oss << "[" << this->nRow << "]";
  oss << " x ";
  oss << "[" << this->nCol << "]";
  oss << " Matrix ";

  return oss.str();
}
