#pragma once

#include "eigen3/Eigen/Core"

class IcpResult {
 public:
  IcpResult() = default;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /** 1 if the result is valid */
  int valid_;

  /** Scan matching result (x,y,theta) */
  Eigen::Vector3d x_;

  /** Number of iterations done */
  int iterations_;
  /** Number of valid correspondence in the end */
  int nvalid_;
  /** Total correspondence error */
  double error_;
};