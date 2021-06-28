#pragma once

#include "csm/algos.h"
#include "icp_result.h"
#include "lidar_data.h"
#include "math_utils.h"

struct GpcCorr {
  // current point
  double p[2];
  // reference point
  double q[2];

  double C[2][2];

  int valid;
};

class IcpParams : public sm_params {
 public:
  IcpParams(const sm_params& base);
  void PLIcp(IcpResult* const result);
  // 循环计算ICP, q0为初始位置(x,y,theta)；
  int IcpLoop(double* const q0, double* const x_new, double* const total_error,
              int* const valid, int* const iterations);

  // laser_ref 和 laser_sens 之间的数据关联
  void FindCorrespondences();

  // 删除多余的点
  // 有两个laser_sens的点和一个laser_ref的点（j1对应）的距离都很近，则删除离得远的
  void KillOutliersDouble();

  // 设置阈值，删除相应的点，如误差小于全部的90%，则有效, 否则删除。
  double KillOutliersTrim();

  void SwapDouble(double* a, double* b);
  void QuickSort(std::vector<double>& array, int begin, int end);

  int ComputeNextEstimate(const double x_old[3], double x_new[3]);

  int TerminationCriterion(const double* delta);
  int Compatible(const int i, const int j);
  int GpcSolve(int K, const std::vector<GpcCorr>& c, const double* x0,
               const double* cov_x0, double* x_out);
  double GpcTotalError(const std::vector<GpcCorr>& co, int n, const double* x);

  // PLICP解析解，见PL-ICP原文附录
  double GpcError(const struct GpcCorr* co, const double* x);

  int PolyGreatestRealRoot(int n, const double* a, double* root);
};