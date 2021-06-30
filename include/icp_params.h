#pragma once

// #include "csm/algos.h"
#include "icp_result.h"
#include "lidar_data.h"
#include "math_utils.h"

struct GpcCorr {
  // current point
  double p[2];
  // reference point
  double q[2];
  Eigen::Vector3d line;

  double C[2][2];

  int valid;
};

class IcpParams {
 public:
  IcpParams();
  void PLIcp(IcpResult* const result);
  // 循环计算ICP, q0为初始位置(x,y,theta)；
  int IcpLoop(const Eigen::Vector3d& x_old, Eigen::Vector3d* const x_new,
              double* const total_error, int* const valid,
              int* const iterations);

  // laser_ref 和 laser_sens 之间的数据关联
  void FindCorrespondences();

  // 删除多余的点
  // 有两个laser_sens的点和一个laser_ref的点（j1对应）的距离都很近，则删除离得远的
  void KillOutliersDouble();

  // 设置阈值，删除相应的点，如误差小于全部的90%，则有效, 否则删除。
  double KillOutliersTrim();

  void SwapDouble(double* a, double* b);
  void QuickSort(std::vector<double>& array, int begin, int end);

  int ComputeNextEstimate(const Eigen::Vector3d x_old,
                          Eigen::Vector3d* const x_new);

  int TerminationCriterion(const Eigen::Vector3d& delta);
  int Compatible(const int i, const int j);
  bool GpcSolve(const int total_size, const std::vector<GpcCorr>& c,
                Eigen::Vector3d* const x_new);

  bool SolveOptimization(const std::vector<GpcCorr>& c,
                         const Eigen::Vector3d& x_old,
                         Eigen::Vector3d* const x_new);

  double GpcTotalError(const std::vector<GpcCorr>& co, const int n,
                       const Eigen::Vector3d& x);

  // PLICP解析解，见PL-ICP原文附录
  double GpcError(const GpcCorr& co, const Eigen::Vector3d& x);

  int PolyGreatestRealRoot(int n, const double* a, double* root);

 public:
  /** First scan ("ref"erence scan) */
  LidarData* laser_ref;
  /** Second scan ("sens"or scan) */
  LidarData* laser_sens;

  /** Where to start */
  double first_guess[3];

  /** Maximum angular displacement between scans (deg)*/
  double max_angular_correction_deg;
  /** Maximum translation between scans (m) */
  double max_linear_correction;

  /** When to stop */
  int max_iterations;
  /** A threshold for stopping. */
  double epsilon_xy;
  /** A threshold for stopping. */
  double epsilon_theta;

  /** Maximum distance for a correspondence to be valid */
  double max_correspondence_dist;
  /** Use smart tricks for finding correspondences. Only influences speed; not
   * convergence. */
  int use_corr_tricks;

  /** Restart if error under threshold (0 or 1)*/
  int restart;
  /** Threshold for restarting */
  double restart_threshold_mean_error;
  /** Displacement for restarting */
  double restart_dt;
  /** Displacement for restarting */
  double restart_dtheta;

  /* Functions concerning discarding correspondences.
     THESE ARE MAGIC NUMBERS -- and they need to be tuned. */

  /** Percentage of correspondences to consider: if 0.9,
      always discard the top 10% of correspondences with more error */
  double outliers_maxPerc;

  /** Parameters describing a simple adaptive algorithm for discarding.
      1) Order the errors.
          2) Choose the percentile according to outliers_adaptive_order.
             (if it is 0.7, get the 70% percentile)
          3) Define an adaptive threshold multiplying outliers_adaptive_mult
             with the value of the error at the chosen percentile.
          4) Discard correspondences over the threshold.

          This is useful to be conservative; yet remove the biggest errors.
  */
  double outliers_adaptive_order; /* 0.7 */
  double outliers_adaptive_mult;  /* 2 */

  /** Do not allow two different correspondences to share a point */
  int outliers_remove_doubles;

  /* Functions that compute and use point orientation for defining matches. */
  /** For now, a very simple max-distance clustering algorithm is used */
  double clustering_threshold;
  /** Number of neighbour rays used to estimate the orientation.*/
  int orientation_neighbourhood;
  /** Discard correspondences based on the angles */
  int do_alpha_test;
  double do_alpha_test_thresholdDeg;

  /** I believe this trick is documented in one of the papers by Guttman (but I
     can't find the reference). Or perhaps I was told by him directly.

          If you already have a guess of the solution, you can compute the polar
     angle of the points of one scan in the new position. If the polar angle is
     not a monotone function of the readings index, it means that the surface is
     not visible in the next position. If it is not visible, then we don't use
     it for matching.

          This is confusing without a picture! To understand what's going on,
     make a drawing in which a surface is not visible in one of the poses.

          Implemented in the function visibilityTest().
  */
  int do_visibility_test;

  /** If 1, use PlICP; if 0, use vanilla ICP. */
  int use_point_to_line_distance;

  /** If 1, the field "true_alpha" is used to compute the incidence
      beta, and the factor (1/cos^2(beta)) used to weight the impact
      of each correspondence. This works fabolously if doing localization,
      that is the first scan has no noise.
          If "true_alpha" is not available, it uses "alpha".
  */
  int use_ml_weights;

  /* If 1, the field "readings_sigma" is used to weight the correspondence by
   * 1/sigma^2 */
  int use_sigma_weights;

  /** Use the method in http://purl.org/censi/2006/icpcov to compute
      the matching covariance. */
  int do_compute_covariance;

  /** Checks that find_correspondences_tricks give the right answer */
  int debug_verify_tricks;

  /** Pose of sensor with respect to robot: used for computing
      the first estimate given the odometry. */
  double laser[3];

  /** Noise in the scan */
  double sigma;

  /** mark as invalid ( = don't use ) rays outside of this interval */
  double min_reading, max_reading;

};