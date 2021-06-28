#include <iostream>

#include <glog/logging.h>
#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/Polynomials>

#include "icp_params.h"
#include "math_utils.h"

IcpParams::IcpParams(const sm_params& base) : sm_params(base) {}

double Norm2D(const double p[2]) { return sqrt(p[0] * p[0] + p[1] * p[1]); }

void ominus_d(const double x[3], double res[3]) {
  double c = cos(x[2]);
  double s = sin(x[2]);
  res[0] = -c * x[0] - s * x[1];
  res[1] = s * x[0] - c * x[1];
  res[2] = -x[2];
}

/** safe if res == x1 */
void oplus_d(const double x1[3], const double x2[3], double res[3]) {
  double c = cos(x1[2]);
  double s = sin(x1[2]);
  double x = x1[0] + c * x2[0] - s * x2[1];
  double y = x1[1] + s * x2[0] + c * x2[1];
  double theta = x1[2] + x2[2];
  res[0] = x;
  res[1] = y;
  res[2] = theta;
}
void pose_diff_d(const double pose2[3], const double pose1[3], double res[3]) {
  double temp[3];
  ominus_d(pose1, temp);
  oplus_d(temp, pose2, res);

  while (res[2] > +M_PI) res[2] -= 2 * M_PI;
  while (res[2] < -M_PI) res[2] += 2 * M_PI;
}

int poly_greatest_real_root(int n, const double* a, double* root) {
  Eigen::VectorXd poly_coeffs(n);
  for (int i = 0; i < n; i++) {
    poly_coeffs(i) = a[i];
  }
  if (n != 5) {
    std::cerr << "ERROR: WRONG DEGREE POLYNOMIAL TO SOLVE." << std::endl;
    return 0;
  }
  Eigen::PolynomialSolver<double, 4> psolve(poly_coeffs);

  // std::cerr<<"dims "<<psolve.roots().rows()<<"
  // "<<psolve.roots().cols()<<std::endl;
  Eigen::Matrix<std::complex<double>, 4, 1, 0, 4, 1> eigen_roots =
      psolve.roots();

  int assigned = 0;
  double lambda = 0;
  for (unsigned int i = 0; i < eigen_roots.size(); i++) {
    if (eigen_roots(i).imag() == 0) {
      if (!assigned || eigen_roots(i).real() > lambda) {
        assigned = 1;
        lambda = eigen_roots(i).real();
      }
    }
  }

  if (!assigned) {
    fprintf(
        stderr,
        "poly_greatest_real_root: Could not find real root for polynomial.\n");
    fprintf(stderr, "polynomial coefficients : ");
    for (int i = 0; i < n; i++) fprintf(stderr, " %lf ", a[i]);
    fprintf(stderr, "\nRoots:\n");

    for (int i = 0; i < n - 1; i++)
      fprintf(stderr, "root z%d = %+.18f + %+.18f i \n", i,
              eigen_roots(i).real(), eigen_roots(i).imag());

    return 0;
  }

  *root = lambda;
  return 1;
}

int gpc_solve(int K, const std::vector<GpcCorr>& c, const double* x0,
              const double* cov_x0, double* x_out) {
  (void)x0;
  (void)cov_x0;
  Eigen::Matrix4d bigM;
  Eigen::Matrix<double, 4, 1> g;
  Eigen::Matrix<double, 2, 4> bigM_k;
  Eigen::Matrix<double, 4, 2> bigM_k_t;
  Eigen::Matrix<double, 2, 2> C_k;
  Eigen::Matrix<double, 2, 1> q_k;
  Eigen::Matrix<double, 4, 1> temp41;
  Eigen::Matrix<double, 2, 2> temp22;
  Eigen::Matrix<double, 2, 2> temp22b;
  Eigen::Matrix<double, 4, 2> temp42;
  Eigen::Matrix<double, 4, 4> temp44;
  Eigen::Matrix<double, 2, 1> temp21;
  Eigen::Matrix<double, 2, 2> temp22c;
  Eigen::Matrix<double, 1, 2> temp12;

  bigM.setZero();
  g.setZero();
  temp42.setZero();

  double d_bigM[4][4] = {
      {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
  double d_g[4] = {0, 0, 0, 0};
  int k;
  for (k = 0; k < K; k++) {
    if (!c[k].valid) continue;

    double C00 = c[k].C[0][0];
    double C01 = c[k].C[0][1];
    double C10 = c[k].C[1][0];
    double C11 = c[k].C[1][1];

    if (C01 != C10) {
      fprintf(stderr, "k=%d; I expect C to be a symmetric matrix.\n", k);
      return 0;
    }

    double qx = c[k].q[0];
    double qy = c[k].q[1];
    double px = c[k].p[0];
    double py = c[k].p[1];

    /* [ C00,  c01,  px C00 + c01 py , c01 px - py C00 ] */
    d_bigM[0][0] += C00;
    d_bigM[0][1] += C01;
    d_bigM[0][2] += +px * C00 + py * C01;
    d_bigM[0][3] += -py * C00 + px * C01;

    /*  [ C10 ,  C11 , py C11 + px C10 , px C11 - py C10 ] */
    d_bigM[1][0] += C10;
    d_bigM[1][1] += C11;
    d_bigM[1][2] += +px * C10 + py * C11;
    d_bigM[1][3] += +px * C11 - py * C10;

    /*Col 1 = [ py C10 + px C00 ]
     Col 2 = [ py C11 + c01 px ]
     Col 3 = [ py (py C11 + px C10) + px (px C00 + c01 py) ]
     Col 4 = [ py (px C11 - py C10) + px (c01 px - py C00) ]
    */
    d_bigM[2][0] += px * C00 + py * C10;
    d_bigM[2][1] += px * C01 + py * C11;
    d_bigM[2][2] +=
        (px * px) * (+C00) + (px * py) * (+C10 + C01) + (py * py) * (+C11);
    d_bigM[2][3] +=
        (px * px) * (+C01) + (px * py) * (-C00 + C11) + (py * py) * (-C10);

    /*Col 1 = [ px C10 - py C00 ]
      Col 2 = [ px C11 - c01 py ]
     Col 3 = [ px (py C11 + px C10) - py (px C00 + c01 py) ]
     Col 4 = [ px (px C11 - py C10) - py (c01 px - py C00) ]*/
    d_bigM[3][0] += -py * C00 + px * C10;
    d_bigM[3][1] += -py * C01 + px * C11;
    d_bigM[3][2] +=
        (px * px) * (+C10) + (px * py) * (-C00 + C11) + (py * py) * (-C01);
    d_bigM[3][3] +=
        (px * px) * (+C11) + (px * py) * (-C10 - C01) + (py * py) * (+C00);

    d_g[0] += C00 * qx + C10 * qy;
    d_g[1] += C01 * qx + C11 * qy;
    d_g[2] += qx * (C00 * px + C01 * py) + qy * (C10 * px + C11 * py);
    d_g[3] += qx * (C00 * (-py) + C01 * px) + qy * (C10 * (-py) + C11 * px);
  }

  {
    unsigned int a, b;
    for (a = 0; a < 4; a++) {
      g(a, 0) = -2 * d_g[a];
    }
    for (a = 0; a < 4; a++) {
      for (b = 0; b < 4; b++) {
        bigM(a, b) = 2 * d_bigM[a][b];
      }
    }
  }

  Eigen::Matrix<double, 2, 2> mA = bigM.block<2, 2>(0, 0);
  Eigen::Matrix<double, 2, 2> mB = bigM.block<2, 2>(0, 2);
  Eigen::Matrix<double, 2, 2> mD = bigM.block<2, 2>(2, 2);

  Eigen::Matrix<double, 2, 2> mS;
  Eigen::Matrix<double, 2, 2> mSa;

  temp22b = mA.inverse();
  /* temp22c = inv(A) * mB           */
  temp22c = temp22b * mB;
  /* temp22 = mB'               */
  temp22 = mB.transpose();
  temp22b = temp22 * temp22c;
  temp22b *= -1.0;
  mS = mD + temp22b;

  /* mSa = mS.inv * mS.det; */
  mSa = mS.inverse();
  mSa *= mS.determinant();

  Eigen::Matrix<double, 2, 1> g1;
  Eigen::Matrix<double, 2, 1> g2;
  Eigen::Matrix<double, 1, 2> g1t;
  Eigen::Matrix<double, 1, 2> g2t;
  Eigen::Matrix<double, 2, 2> mAi;
  Eigen::Matrix<double, 2, 2> mBt;

  g1(0, 0) = g(0, 0);
  g1(1, 0) = g(1, 0);
  g2(0, 0) = g(2, 0);
  g2(1, 0) = g(3, 0);
  g1t = g1.transpose();
  g2t = g2.transpose();
  mBt = mB.transpose();
  mAi = mA.inverse();

  Eigen::Matrix<double, 1, 2> m1t;
  Eigen::Matrix<double, 2, 1> m1;
  Eigen::Matrix<double, 1, 2> m2t;
  Eigen::Matrix<double, 2, 1> m2;
  Eigen::Matrix<double, 1, 2> m3t;
  Eigen::Matrix<double, 2, 1> m3;

  /* m1t = g1t*mAi*mB */
  temp12 = g1t * mAi;
  m1t = temp12 * mB;

  m1 = m1t.transpose();
  /*     m2t = m1t*mSa    */
  m2t = m1t * mSa;
  m2 = m2t.transpose();
  /* m3t = g2t*mSa     */
  m3t = g2t * mSa;
  m3 = m3t.transpose();

  double p[3] = {m2t.dot(m2) - 2 * m2t.dot(m3) + m3t.dot(m3),
                 4 * m2t.dot(m1) - 8 * m2t.dot(g2) + 4 * g2t.dot(m3),
                 4 * m1t.dot(m1) - 8 * m1t.dot(g2) + 4 * g2t.dot(g2)};

  double l[3] = {mS.determinant(), 2 * mS(0, 0) + 2 * mS(1, 1), 4};

  /* q = p - l^2        */
  double q[5] = {p[0] - (l[0] * l[0]), p[1] - (2 * l[1] * l[0]),
                 p[2] - (l[1] * l[1] + 2 * l[0] * l[2]), -(2 * l[2] * l[1]),
                 -(l[2] * l[2])};

  double lambda;
  if (!poly_greatest_real_root(5, q, &lambda)) return 0;

  Eigen::Matrix<double, 4, 4> W;
  W.setZero();
  W(2, 2) = 1.0;
  W(3, 3) = 1.0;
  Eigen::Matrix<double, 4, 1> x;
  W *= (2 * lambda);

  bigM += W;
  temp44 = bigM.inverse();
  x = temp44 * g;
  x *= -1.0;

  x_out[0] = x(0, 0);
  x_out[1] = x(1, 0);
  x_out[2] = atan2(x(3, 0), x(2, 0));

  return 1;
}

double gpc_error(const struct GpcCorr* co, const double* x) {
  double c = cos(x[2]);
  double s = sin(x[2]);
  double e[2];
  e[0] = c * (co->p[0]) - s * (co->p[1]) + x[0] - co->q[0];
  e[1] = s * (co->p[0]) + c * (co->p[1]) + x[1] - co->q[1];
  double this_error = e[0] * e[0] * co->C[0][0] +
                      2 * e[0] * e[1] * co->C[0][1] + e[1] * e[1] * co->C[1][1];

  if (0) /* due to limited numerical precision, error might be negative */
    if (this_error < 0) {
      fprintf(
          stderr,
          "Something fishy: error = %lf e = [%lf %lf]  C = [%lf,%lf;%lf,%lf]\n",
          this_error, e[0], e[1], co->C[0][0], co->C[0][1], co->C[1][0],
          co->C[1][1]);
    }
  return this_error;
}

double gpc_total_error(const std::vector<GpcCorr>& co, int n, const double* x) {
  int i;
  double error = 0;
  for (i = 0; i < n; i++) {
    if (!co[i].valid) continue;
    error += gpc_error(&(co.at(i)), x);
  }
  if (0) /* due to limited numerical precision, error might be negative */
    if (error < 0) {
      fprintf(stderr, "Something fishy!\n");
    }
  return error;
}

int IcpParams::Compatible(int i, int j) {
  // 目前始终为 true
  if (!this->do_alpha_test) {
    return 1;
  }

  double theta0 = 0; /* FIXME */
  if ((this->laser_sens->alpha_valid[i] == 0) ||
      (this->laser_ref->alpha_valid[j] == 0))
    return 1;

  double alpha_i = this->laser_sens->alpha[i];
  double alpha_j = this->laser_ref->alpha[j];
  double tolerance = MathUtils::Deg2Rad(this->do_alpha_test_thresholdDeg);

  /** FIXME remove alpha test */
  double theta = MathUtils::AngleDiff(alpha_j, alpha_i);
  if (fabs(MathUtils::AngleDiff(theta, theta0)) >
      tolerance + MathUtils::Deg2Rad(this->max_angular_correction_deg)) {
    return 0;
  } else {
    return 1;
  }
}

void IcpParams::FindCorrespondences() {
  // const LDP laser_ref  = this->laser_ref;
  // const LDP laser_sens = this->laser_sens;
  LidarData* laser_ref = reinterpret_cast<LidarData*>(this->laser_ref);
  LidarData* laser_sens = reinterpret_cast<LidarData*>(this->laser_sens);

  for (int i = 0; i < laser_sens->nrays; ++i) {
    // Check
    if (!laser_sens->ValidRay(i)) {
      laser_sens->SetNullCorrespondence(i);
      continue;
    }

    double* p_i_w = laser_sens->points_w[i].p;

    int j1 = -1;
    double best_dist = 10000;

    int from, to, start_cell;
    laser_ref->PossibleInterval(p_i_w, this->max_angular_correction_deg,
                                this->max_linear_correction, &from, &to,
                                &start_cell);

    // 选择最邻近的点1(j1)
    for (int j = from; j <= to; j++) {
      if (!laser_ref->ValidRay(j)) {
        continue;
      }

      double dist = MathUtils::DistanceSquared(p_i_w, laser_ref->points[j].p);
      if (dist > pow(this->max_correspondence_dist, 2)) {
        continue;
      }

      if ((-1 == j1) || (dist < best_dist)) {
        // 目前Compatible始终为true, 没哟使用alpha_test
        if (Compatible(i, j)) {
          j1 = j;
          best_dist = dist;
        }
      }
    }

    // j1 == -1, 则没有满足项;
    // j1不能是第一个和最后一个, 不处理极端情况
    if (j1 == -1 || j1 == 0 || (j1 == (laser_ref->nrays - 1))) {
      laser_sens->SetNullCorrespondence(i);
      continue;
    }

    const int j2up = laser_ref->NextValidUp(j1);
    const int j2down = laser_ref->NextValidDown(j1);
    if ((j2up == -1) && (j2down == -1)) {
      laser_sens->SetNullCorrespondence(i);
      continue;
    }

    int j2;
    if (j2up == -1) {
      j2 = j2down;
    } else if (j2down == -1) {
      j2 = j2up;
    } else {
      // 取距离近的
      const double dist_up =
          MathUtils::DistanceSquared(p_i_w, laser_ref->points[j2up].p);
      const double dist_down =
          MathUtils::DistanceSquared(p_i_w, laser_ref->points[j2down].p);
      j2 = dist_up < dist_down ? j2up : j2down;
    }

    laser_sens->SetCorrespondence(i, j1, j2);
    laser_sens->corr[i].dist2_j1 = best_dist;
    // use_point_to_line_distance default TRUE.
    laser_sens->corr[i].type = this->use_point_to_line_distance
                                   ? correspondence::corr_pl
                                   : correspondence::corr_pp;
  }
}

void IcpParams::KillOutliersDouble() {
  double threshold = 3; /* TODO: add as configurable */

  LDP laser_ref = this->laser_ref;
  LDP laser_sens = this->laser_sens;

  std::vector<double> dist2_i(laser_sens->nrays, 0.0);
  std::vector<double> dist2_j(laser_ref->nrays, 1000000);

  for (int i = 0; i < laser_sens->nrays; i++) {
    if (!laser_sens->corr[i].valid) {
      continue;
    }
    const int j1 = laser_sens->corr[i].j1;
    dist2_i[i] = laser_sens->corr[i].dist2_j1;
    dist2_j[j1] = (std::min)(dist2_j[j1], dist2_i[i]);
  }

  int nkilled = 0;
  for (int i = 0; i < laser_sens->nrays; i++) {
    if (!laser_sens->corr[i].valid) {
      continue;
    }
    const int j1 = laser_sens->corr[i].j1;
    // 有两个laser_sens的点和一个laser_ref的点（j1对应）的距离都很近，则删除离得远的
    if (dist2_i[i] > (threshold * threshold) * dist2_j[j1]) {
      laser_sens->corr[i].valid = 0;
      nkilled++;
    }
  }
  LOG(INFO) << "kill_outliers_double: killed " << nkilled << " correspondences";
}

void IcpParams::SwapDouble(double* a, double* b) {
  double t = *a;
  *a = *b;
  *b = t;
}

/** Code taken from Wikipedia */
// TODO 替换为标准库
void IcpParams::QuickSort(std::vector<double>& array, int begin, int end) {
  if (end > begin) {
    double pivot = array[begin];
    int l = begin + 1;
    int r = end + 1;
    while (l < r) {
      if (array[l] < pivot) {
        l++;
      } else {
        r--;
        SwapDouble(&(array[l]), &(array[r]));
      }
    }
    l--;
    SwapDouble(&(array[begin]), &(array[l]));
    if (l > begin) {
      QuickSort(array, begin, l);
    }
    if (end > r) {
      QuickSort(array, r, end);
    }
  }
}

double IcpParams::KillOutliersTrim() {
  LDP laser_ref = this->laser_ref;
  LDP laser_sens = this->laser_sens;

  /* dist2, indexed by k, contains the error for the k-th correspondence */
  int k = 0;

  // dist2 是顺序所有有效间距，dist会跳过没有的
  std::vector<double> dist2(laser_sens->nrays, 0.0);
  std::vector<double> dist(laser_sens->nrays, 0.0);

  /* for each point in laser_sens */
  for (int i = 0; i < laser_sens->nrays; i++) {
    /* which has a valid correspondence */
    if (!laser_sens->corr[i].valid) {
      dist[i] = std::numeric_limits<double>::quiet_NaN();
      continue;
    }
    double* p_i_w = laser_sens->points_w[i].p;

    int j1 = laser_sens->corr[i].j1;
    int j2 = laser_sens->corr[i].j2;
    /* Compute the distance to the corresponding segment */
    dist[i] = MathUtils::DistToSegment(laser_ref->points[j1].p,
                                       laser_ref->points[j2].p, p_i_w);
    dist2[k] = dist[i];
    k++;
  }

  /* two errors limits are defined: */
  /* In any case, we don't want more than outliers_maxPerc% */
  int order = (int)floor(k * (this->outliers_maxPerc));
  // Looks like MinMax
  order = (std::max)(0, (std::min)(order, k - 1));

  /* The dists for the correspondence are sorted
     in ascending order */
  QuickSort(dist2, 0, k - 1);
  // outliers_maxPerc default 0.9, 这里取90%的误差作为误差阈值
  double error_limit1 = dist2[order];

  /* Then we take a order statics (o*K) */
  /* And we say that the error must be less than alpha*dist(o*K) */
  int order2 = (int)floor(k * this->outliers_adaptive_order);
  order2 = (std::max)(0, (std::min)(order2, k - 1));
  double error_limit2 = this->outliers_adaptive_mult * dist2[order2];

  double error_limit = (std::min)(error_limit1, error_limit2);

  LOG(INFO) << "icp_outliers: maxPerc " << this->outliers_maxPerc
            << " error_limit: fix " << error_limit1 << " adaptive "
            << error_limit2;

  double total_error = 0;
  int nvalid = 0;
  for (int i = 0; i < laser_sens->nrays; i++) {
    if (!laser_sens->corr[i].valid) {
      continue;
    }
    if (dist[i] > error_limit) {
      laser_sens->corr[i].valid = 0;
      laser_sens->corr[i].j1 = -1;
      laser_sens->corr[i].j2 = -1;
    } else {
      nvalid++;
      total_error += dist[i];
    }
  }

  LOG(INFO) << "icp_outliers: valid " << nvalid << " / " << k
            << "(limit: " << error_limit
            << ") mean error = " << total_error / nvalid;
  return total_error;
}

int IcpParams::TerminationCriterion(const double* delta) {
  double a = Norm2D(delta);
  double b = fabs(delta[2]);
  return (a < this->epsilon_xy) && (b < this->epsilon_theta);
}

int IcpParams::IcpLoop(double* const q0, double* const x_new,
                       double* const total_error, int* const valid,
                       int* const iterations) {
  if (MathUtils::AnyNan(q0, 3)) {
    // LOG(ERROR) << "icp_loop: Initial pose contains nan: " <<
    // friendly_pose(q0);
    return 0;
  }

  // 新一帧的位置
  LidarData* laser_sens = reinterpret_cast<LidarData*>(this->laser_sens);
  double x_old[3], delta[3], delta_old[3] = {0, 0, 0};

  // x_old <- q0, q0是迭代的起点
  MathUtils::Copy(q0, 3, x_old);

  std::vector<unsigned int> hashes(this->max_iterations, 0);

  // sm_debug("icp: starting at  q0 =  %s  \n", friendly_pose(x_old));

  int all_is_okay = 1;

  int iteration;
  for (iteration = 0; iteration < this->max_iterations; iteration++) {
    /** Compute laser_sens's points in laser_ref's coordinates
        by roto-translating by x_old */
    laser_sens->ComputeWorldCoords(x_old);

    /** Find correspondences (the naif or smart way) */
    // TODO
    // if (this->use_corr_tricks)
    //   find_correspondences_tricks(this);
    // else
    FindCorrespondences();

    /** If debug_verify_tricks, make sure that find_correspondences_tricks()
        and find_correspondences() return the same answer */
    // if (this->debug_verify_tricks) debug_correspondences(this);

    /* If not many correspondences, bail out */
    const int num_corr = laser_sens->NumValidCorrespondences();
    double fail_perc = 0.05;
    if (num_corr < fail_perc * laser_sens->nrays) { /* TODO: arbitrary */
      LOG(ERROR) << "	: before trimming, only " << num_corr
                 << " correspondences.";
      all_is_okay = 0;
      break;
    }

    /* Kill some correspondences (using dubious algorithm) */
    if (this->outliers_remove_doubles) {
      KillOutliersDouble();
    }

    const int num_corr2 = laser_sens->NumValidCorrespondences();

    double error = this->KillOutliersTrim();

    const int num_corr_after = laser_sens->NumValidCorrespondences();

    LOG(INFO) << "Laser_sens has " << laser_sens->nrays
              << " rays valid. Corr found: " << num_corr
              << ", after double cut: " << num_corr2
              << ", after adaptive cut: " << num_corr_after;

    *total_error = error;
    *valid = num_corr_after;

    LOG(INFO) << "Icp_loop: total error: " << *total_error << "  valid "
              << *valid << "   mean = " << *total_error / *valid;

    /* If not many correspondences, bail out */
    if (num_corr_after < fail_perc * laser_sens->nrays) {
      LOG(INFO) << "  icp_loop: failed: after trimming, only " << num_corr_after
                << "correspondences.";
      all_is_okay = 0;
      break;
    }

    /* Compute next estimate based on the correspondences */
    if (!ComputeNextEstimate(x_old, x_new)) {
      // sm_error("  icp_loop: Cannot compute next estimate.\n");
      all_is_okay = 0;
      break;
    }
    pose_diff_d(x_new, x_old, delta);

    // {
    //   sm_debug(
    //       "  icp_loop: killing. laser_sens has %d/%d rays valid,  %d corr "
    //       "found -> %d after double cut -> %d after adaptive cut \n",
    //       count_equal(laser_sens->valid, laser_sens->nrays, 1),
    //       laser_sens->nrays, num_corr, num_corr2, num_corr_after);
    // }
    /** Checks for oscillations */
    hashes[iteration] = laser_sens->CorrHash();

    // {
    //   sm_debug("  icp_loop: it. %d  hash=%d nvalid=%d mean error = %f, x_new=
    //            % s\n
    //            ", iteration, hashes[iteration], *valid, *total_error /
    //            *valid, friendly_pose(x_new));
    // }

    /** PLICP terminates in a finite number of steps! */
    if (this->use_point_to_line_distance) {
      int loop_detected = 0; /* TODO: make function */
      int a;
      for (a = iteration - 1; a >= 0; a--) {
        if (hashes[a] == hashes[iteration]) {
          // sm_debug("icpc: oscillation detected (cycle length = %d)\n",
          //          iteration - a);
          loop_detected = 1;
          break;
        }
      }
      if (loop_detected) {
        // egsl_pop_named("icp_loop iteration");
        break;
      }
    }

    /* This termination criterium is useless when using
       the point-to-line-distance; however, we put it here because
       one can choose to use the point-to-point distance. */
    if (TerminationCriterion(delta)) {
      // egsl_pop_named("icp_loop iteration");
      break;
    }

    MathUtils::Copy(x_new, 3, x_old);
    MathUtils::Copy(delta, 3, delta_old);

    // egsl_pop_named("icp_loop iteration");
  }

  *iterations = iteration + 1;

  return all_is_okay;
}

void IcpParams::PLIcp(IcpResult* const result) {
  result->valid = 0;

  LidarData* laser_ref = reinterpret_cast<LidarData*>(this->laser_ref);
  LidarData* laser_sens = reinterpret_cast<LidarData*>(this->laser_sens);

  if (!laser_ref->ValidLidar() || !laser_sens->ValidLidar()) {
    return;
  }
  laser_ref->InvalidIfOutside(this->min_reading, this->max_reading);
  laser_sens->InvalidIfOutside(this->min_reading, this->max_reading);

  // TODO
  // if(this->use_corr_tricks ||this->debug_verify_tricks)
  // 	ld_create_jump_tables(laser_ref);

  // Origin laser is range and bearing, here compute as Cartesian.
  laser_ref->ComputeCartesian();
  laser_sens->ComputeCartesian();

  // TODO
  // if (this->do_alpha_test) {
  //   ld_simple_clustering(laser_ref,this->clustering_threshold);
  //   ld_compute_orientation(laser_ref,this->orientation_neighbourhood,
  // this->sigma);
  //   ld_simple_clustering(laser_sens,this->clustering_threshold);
  //   ld_compute_orientation(laser_sens,this->orientation_neighbourhood,
  // this->sigma);
  // }

  Eigen::Vector3d x_new;
  // 	gsl_vector * x_old = vector_from_array(3,this->first_guess);
  Eigen::Vector3d x_old(this->first_guess[0], this->first_guess[1],
                        this->first_guess[2]);

  double error;
  int iterations;
  int nvalid;

  if (!IcpLoop(x_old.data(), x_new.data(), &error, &nvalid, &iterations)) {
    // sm_error("icp: ICP failed for some reason. \n");
    result->valid = 0;
    result->iterations = iterations;
    result->nvalid = 0;
  } else {
  }
}

int IcpParams::ComputeNextEstimate(const double x_old[3], double x_new[3]) {
  LDP laser_ref = this->laser_ref;
  LDP laser_sens = this->laser_sens;

  // struct GpcCorr c[laser_sens->nrays];
  struct GpcCorr dummy;
  std::vector<GpcCorr> c(laser_sens->nrays, dummy);

  int k = 0;
  for (int i = 0; i < laser_sens->nrays; i++) {
    if (!laser_sens->valid[i]) {
      continue;
    }

    if (!laser_sens->corr[i].valid) {
      continue;
    }

    int j1 = laser_sens->corr[i].j1;
    int j2 = laser_sens->corr[i].j2;

    c[k].valid = 1;

    if (laser_sens->corr[i].type == correspondence::corr_pl) {
      c[k].p[0] = laser_sens->points[i].p[0];
      c[k].p[1] = laser_sens->points[i].p[1];
      c[k].q[0] = laser_ref->points[j1].p[0];
      c[k].q[1] = laser_ref->points[j1].p[1];

      /** TODO: here we could use the estimated alpha */
      double diff[2];
      diff[0] = laser_ref->points[j1].p[0] - laser_ref->points[j2].p[0];
      diff[1] = laser_ref->points[j1].p[1] - laser_ref->points[j2].p[1];
      double one_on_norm = 1 / sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
      double normal[2];
      normal[0] = +diff[1] * one_on_norm;
      normal[1] = -diff[0] * one_on_norm;

      double cos_alpha = normal[0];
      double sin_alpha = normal[1];

      c[k].C[0][0] = cos_alpha * cos_alpha;
      c[k].C[1][0] = c[k].C[0][1] = cos_alpha * sin_alpha;
      c[k].C[1][1] = sin_alpha * sin_alpha;

    } else {
      c[k].p[0] = laser_sens->points[i].p[0];
      c[k].p[1] = laser_sens->points[i].p[1];

      MathUtils::ProjectionOnSegment(laser_ref->points[j1].p,
                                     laser_ref->points[j2].p,
                                     laser_sens->points_w[i].p, c[k].q);

      /* Identity matrix */
      c[k].C[0][0] = 1;
      c[k].C[1][0] = 0;
      c[k].C[0][1] = 0;
      c[k].C[1][1] = 1;
    }

    double factor = 1;

    /* Scale the correspondence weight by a factor concerning the
       information in this reading. */
    if (this->use_ml_weights) {
      int have_alpha = 0;
      double alpha = 0;
      if (!isnan(laser_ref->true_alpha[j1])) {
        alpha = laser_ref->true_alpha[j1];
        have_alpha = 1;
      } else if (laser_ref->alpha_valid[j1]) {
        alpha = laser_ref->alpha[j1];
        ;
        have_alpha = 1;
      } else
        have_alpha = 0;

      if (have_alpha) {
        double pose_theta = x_old[2];
        /** Incidence of the ray
                Note that alpha is relative to the first scan (not the world)
                and that pose_theta is the angle of the second scan with
                respect to the first, hence it's ok. */
        double beta = alpha - (pose_theta + laser_sens->theta[i]);
        factor = 1 / pow(cos(beta), 2);
      } else {
        static int warned_before = 0;
        if (!warned_before) {
          LOG(WARNING)
              << "Param use_ml_weights was active, but not valid alpha[] or "
                 "true_alpha[]. Perhaps, if this is a single ray not having "
                 "alpha, you should mark it as inactive.";
          warned_before = 1;
        }
      }
    }

    /* Weight the points by the sigma in laser_sens */
    if (this->use_sigma_weights) {
      if (!isnan(laser_sens->readings_sigma[i])) {
        factor *= 1 / pow(laser_sens->readings_sigma[i], 2);
      } else {
        static int warned_before = 0;
        if (!warned_before) {
          LOG(WARNING) << "Param use_sigma_weights was active, but the field "
                          "readings_sigma[] was not filled in.";
        }
      }
    }

    c[k].C[0][0] *= factor;
    c[k].C[1][0] *= factor;
    c[k].C[0][1] *= factor;
    c[k].C[1][1] *= factor;

    k++;
  }

  /* TODO: use prior for odometry */
  double std = 0.11;
  const double inv_cov_x0[9] = {
      1 / (std * std), 0, 0, 0, 1 / (std * std), 0, 0, 0, 0};

  int ok = gpc_solve(k, c, 0, inv_cov_x0, x_new);
  if (!ok) {
    LOG(ERROR) << "gpc_solve_valid failed";
    return 0;
  }

  double old_error = gpc_total_error(c, k, x_old);
  double new_error = gpc_total_error(c, k, x_new);

  // sm_debug("\tcompute_next_estimate: old error: %f  x_old= %s \n",
  // old_error,
  //          friendly_pose(x_old));
  // sm_debug("\tcompute_next_estimate: new error: %f  x_new= %s \n",
  // new_error,
  //          friendly_pose(x_new));
  // sm_debug("\tcompute_next_estimate: new error - old_error: %g \n",
  //          new_error - old_error);

  double epsilon = 0.000001;
  if (new_error > old_error + epsilon) {
    // sm_error(
    //     "\tcompute_next_estimate: something's fishy here! Old error: %lf
    //     new " "error: %lf  x_old %lf %lf %lf x_new %lf %lf %lf\n",
    //     old_error, new_error, x_old[0], x_old[1], x_old[2], x_new[0],
    //     x_new[1], x_new[2]);
  }

  return 1;
}