#include "TransitionMatrix.h"
#include <fstream>
#include <mutex>


TransitionMatrix::TransitionMatrix(const std::string& input_file, Heme* hemes, int nHemes)
: hemes_(hemes), nHemes_(nHemes) {}

void TransitionMatrix::initialize_data(const std::string& input_file)
{
  std::ifstream fin(input_file);

  transfer_ = new double*[nHemes_];
  injection_ = new double[nHemes_];
  ejection_ = new double[nHemes_];
  currentState_ = new bool[nHemes_];
  animationState_ = new bool[nHemes_];

  int row, col;
  double rate;
  for(int i = 0; i < nHemes_; ++i) {
    transfer_[i] = new double[nHemes_];

    for(int j = 0; j < nHemes_; ++j) {
      if(i == j)
        transfer_[i][j] = 0;
      else
      {
        fin >> row;
        fin >> col;
        fin >> rate;
        transfer_[row-1][col-1] = rate;
      }
    }

    injection_[i] = 0;
    ejection_[i] = 0;
    currentState_[i] = false;
    animationState_[i] = false;

    //Manually change this to desired values
    injection_[0] = 150000;
    ejection_[nHemes_-1] = 150000;
  }



}
