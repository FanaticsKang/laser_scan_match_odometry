#pragma once
#include "csm/laser_data.h"

class LidarData : public laser_data {
 public:
  LidarData(const laser_data& base) : laser_data(base){};
  int CountEqual(const int* v, int n, int value) {
    int num = 0;
    for (int i = 0; i < n; ++i) {
      if (value == v[i]) {
        ++num;
      }
    }
    return num;
  };
};
