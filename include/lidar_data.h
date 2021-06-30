#pragma once
// #include "csm/laser_data.h"
#include "correspondence.h"
#include "point2d.h"

#include <vector>

class LidarData {
 public:
  LidarData();
  void Initialization(const size_t point_size);
  int CountEqual(const int* v, int n, int value);
  bool ValidLidar();
  void InvalidIfOutside(const double min_reading, const double max_reading);
  bool ValidRay(const int i) const;
  void ComputeCartesian();
  void ComputeWorldCoords(const Eigen::Vector3d& pose);

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

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  int nrays;
  double min_theta;
  double max_theta;

  std::vector<double> theta;
  std::vector<int> valid;
  std::vector<double> readings;
  std::vector<int> cluster;
  std::vector<double> alpha;
  std::vector<double> cov_alpha;
  std::vector<int> alpha_valid;
  std::vector<double> readings_sigma;
  std::vector<double> true_alpha;

  std::vector<Correspondence> corr;

  double true_pose[3];
  double odometry[3];
  double estimate[3];

  /** Cartesian representation */
  std::vector<Point2d> points;
  /** Cartesian representation, in "world" (laser_ref) coordinates.
      Computed using ld_compute_world_coords() */
  std::vector<Point2d> points_w;
};
