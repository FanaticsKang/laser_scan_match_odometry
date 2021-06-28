#pragma once
#include "csm/laser_data.h"

class LidarData : public laser_data {
 public:
  LidarData(const laser_data& base) : laser_data(base){};
  int CountEqual(const int* v, int n, int value);
  bool ValidLidar();
  void InvalidIfOutside(const double min_reading, const double max_reading);
  bool ValidRay(const int i) const;
  void ComputeCartesian();
  void ComputeWorldCoords(const double* pose);
  void ComputeWorldCoords(Eigen::Vector3d* const pose);

  void SetNullCorrespondence(const int i);

  // 该激光点可能的起始点, start cell 为中心值，from 和 to 在start cell左右.
  void PossibleInterval(const double* p_i_w, double max_angular_correction_deg,
                        double max_linear_correction, int* from, int* to,
                        int* start_cell);

  void SetCorrespondence(int i, int j1, int j2);
  int NumValidCorrespondences();

  unsigned int CorrHash();

  int MinMax(int from, int to, int x) {
    return (std::max)((std::min)(x, to), from);
  }

  int NextValid(int i, int dir) {
    int j = i + dir;
    // 在范围内，且有效
    while ((j < this->nrays) && (j >= 0) && !this->ValidRay(j)) {
      j += dir;
    }
    return this->ValidRay(j) ? j : -1;
  }

  int NextValidUp(int i) { return NextValid(i, +1); }

  int NextValidDown(int i) { return NextValid(i, -1); }
};
