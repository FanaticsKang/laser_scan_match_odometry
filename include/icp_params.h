#pragma once

#include "csm/algos.h"
#include "icp_result.h"
#include "lidar_data.h"
#include "math_utils.h"

struct GpcCorr {
  double p[2];
  double q[2];

  double C[2][2];

  int valid;
};

class IcpParams : public sm_params {
 public:
  IcpParams(const sm_params& base);
  void PLIcp(IcpResult* const result);
  int IcpLoop(double* const q0, double* const x_new, double* const total_error,
              int* const valid, int* const iterations);

  void FindCorrespondences();
  
  // 删除多余的点
  // 有两个laser_sens的点和一个laser_ref的点（j1对应）的距离都很近，则删除离得远的
  void KillOutliersDouble();
  double KillOutliersTrim();

  void SwapDouble(double* a, double* b);
  void QuickSort(std::vector<double>& array, int begin, int end);

  int ComputeNextEstimate(const double x_old[3], double x_new[3]);

  int TerminationCriterion(const double* delta);
  int Compatible(const int i, const int j);
};