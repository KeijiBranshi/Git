#ifndef TRANSITIONMATRIX_H
#define TRANSITIONMATRIX_H

class TransitionMatrix {
private:
  double time_ = 0;

  int nHemes_;
  Heme* hemes_;

  double** transfer_;
  double* injection_;
  double* ejection_;

  bool* animationState_;
  bool* currentState_;

  std::queue<Event> eventBuffer_;

public:
  TransitionMatrix(const std::string& input_file, Heme* hemes, int nHeme);

  /*Simulation Methods*/
  void initialize_data(const std::string& input_file);
  void run();
  void runIteration();

  /*Accessors, Modifiers, and animation stuff*/
  double getTransferRate(int from_heme, int to_heme);
  void setCurrentState(int heme, bool val);
  void set AnimationState(int heme, bool val);
  public Even dequeueEvent();
  void clearBuffer();

  double getTimeElapsed();
  int getNumHemes();
  Heme getHemeObject(int heme);
};

#endif
