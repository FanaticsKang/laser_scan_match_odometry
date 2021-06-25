#pragma once
#include "csm/laser_data.h"

class LidarData : public laser_data {
 public:
  LidarData(const laser_data& base) : laser_data(base){};
  int CountEqual(const int* v, int n, int value);
  bool ValidLidar();
  void InvalidIfOutside(const double min_reading, const double max_reading);
  bool ValidRay(const int i);
  void ComputeCartesian();
};
