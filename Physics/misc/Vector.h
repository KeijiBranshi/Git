#ifndef VECTOR_H
#define VECTOR_H

#include <vector>

class Vector  {
private:
  std::vector<double> values;
public:
  Vector();
  Vector(double[], int);
  Vector(std::vector<double>);
  ~Vector();

  Vector operator+(Vector);
  Vector operator*(Vector);
  Vector operator==(Vector);
}

#endif
