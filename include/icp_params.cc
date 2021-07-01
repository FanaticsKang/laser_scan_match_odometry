#include <iostream>

#include <glog/logging.h>
#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/Polynomials>
#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "toc.h"

#include "edge_point_to_line.h"
#include "icp_params.h"
#include "math_utils.h"

IcpParams::IcpParams() {}

int IcpParams::PolyGreatestRealRoot(int n, const double* a, double* root) {
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

bool IcpParams::SolveOptimization(
    const std::vector<GpcCorrespondence>& connections,
    const Eigen::Vector3d& x_old, Eigen::Vector3d* const x_new) {
  // TODO some error in my g2o lib
  static g2o::SparseOptimizer optimizer;
  optimizer.edges().clear();
  optimizer.vertices().clear();

  typedef g2o::BlockSolver<g2o::BlockSolverTraits<3, 1> > BlockSolver_3_1;
  // 注意
  BlockSolver_3_1::LinearSolverType* linearSolver;

  linearSolver = new g2o::LinearSolverDense<BlockSolver_3_1::PoseMatrixType>();

  BlockSolver_3_1* solver_ptr = new BlockSolver_3_1(linearSolver);

  g2o::OptimizationAlgorithmGaussNewton* solver =
      new g2o::OptimizationAlgorithmGaussNewton(solver_ptr);
  optimizer.setAlgorithm(solver);

  g2o::SE2 se2(x_old);
  g2o::VertexSE2* vertex_se2 = new g2o::VertexSE2();
  vertex_se2->setEstimate(se2);
  vertex_se2->setId(0);
  vertex_se2->setFixed(false);
  optimizer.addVertex(vertex_se2);

  // const double th_huber = sqrt(3.841);

  for (auto& one_connection : connections) {
    Eigen::Vector2d point(one_connection.p[0], one_connection.p[1]);
    EdgePointToLine* edge = new EdgePointToLine(point, one_connection.line);
    edge->setVertex(
        0, static_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(0)));
    edge->setMeasurement(0);
    edge->setInformation(Eigen::Matrix<double, 1, 1>::Identity());

    // g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
    // edge->setRobustKernel(rk);
    // rk->setDelta(th_huber);

    optimizer.addEdge(edge);
  }
  // optimizer.setVerbose(true);
  optimizer.initializeOptimization();

  Alpha::TicToc tic;
  optimizer.optimize(3);
  std::cout << "Optimization: " << tic.TocMicroseconds()
            << " ms, edge size: " << optimizer.edges().size() << std::endl;

  g2o::VertexSE2* result = static_cast<g2o::VertexSE2*>(optimizer.vertex(0));
  *x_new = result->estimate().toVector();
  return true;
}

bool IcpParams::GpcSolve(const int total_size,
                         const std::vector<GpcCorrespondence>& connections,
                         Eigen::Vector3d* const x_out) {
  Eigen::Matrix4d d_bigM = Eigen::Matrix4d::Zero();
  Eigen::Vector4d d_g = Eigen::Vector4d::Zero();

  for (int k = 0; k < total_size; ++k) {
    if (!connections[k].valid) {
      continue;
    }

    const double C00 = connections[k].C[0][0];
    const double C01 = connections[k].C[0][1];
    const double C10 = connections[k].C[1][0];
    const double C11 = connections[k].C[1][1];

    if (C01 != C10) {
      std::cout << "I expect C to be a symmetric matrix, k = " << k
                << ",  C =\n"
                << connections[k].C[0][0] << std::endl;
      return false;
    }

    // 取出了参考点（q）, 当前点（p）和一个方向矩阵（C）
    const double qx = connections[k].q[0];
    const double qy = connections[k].q[1];
    const double px = connections[k].p[0];
    const double py = connections[k].p[1];

    //  [ C00,  c01,  px C00 + c01 py , c01 px - py C00 ]
    d_bigM(0, 0) += C00;
    d_bigM(0, 1) += C01;
    d_bigM(0, 2) += +px * C00 + py * C01;
    d_bigM(0, 3) += -py * C00 + px * C01;

    // [ C10 ,  C11 , py C11 + px C10 , px C11 - py C10 ]
    d_bigM(1, 0) += C10;
    d_bigM(1, 1) += C11;
    d_bigM(1, 2) += +px * C10 + py * C11;
    d_bigM(1, 3) += +px * C11 - py * C10;

    //  Col 1 = [ py C10 + px C00 ]
    //  Col 2 = [ py C11 + c01 px ]
    //  Col 3 = [ py (py C11 + px C10) + px (px C00 + c01 py) ]
    //  Col 4 = [ py (px C11 - py C10) + px (c01 px - py C00) ]
    d_bigM(2, 0) += px * C00 + py * C10;
    d_bigM(2, 1) += px * C01 + py * C11;
    d_bigM(2, 2) +=
        (px * px) * (+C00) + (px * py) * (+C10 + C01) + (py * py) * (+C11);
    d_bigM(2, 3) +=
        (px * px) * (+C01) + (px * py) * (-C00 + C11) + (py * py) * (-C10);

    // Col 1 = [ px C10 - py C00 ]
    // Col 2 = [ px C11 - c01 py ]
    // Col 3 = [ px (py C11 + px C10) - py (px C00 + c01 py) ]
    // Col 4 = [ px (px C11 - py C10) - py (c01 px - py C00) ]
    d_bigM(3, 0) += -py * C00 + px * C10;
    d_bigM(3, 1) += -py * C01 + px * C11;
    d_bigM(3, 2) +=
        (px * px) * (+C10) + (px * py) * (-C00 + C11) + (py * py) * (-C01);
    d_bigM(3, 3) +=
        (px * px) * (+C11) + (px * py) * (-C10 - C01) + (py * py) * (+C00);

    d_g[0] += C00 * qx + C10 * qy;
    d_g[1] += C01 * qx + C11 * qy;
    d_g[2] += qx * (C00 * px + C01 * py) + qy * (C10 * px + C11 * py);
    d_g[3] += qx * (C00 * (-py) + C01 * px) + qy * (C10 * (-py) + C11 * px);
  }

  Eigen::Matrix4d bigM = 2 * d_bigM;
  Eigen::Vector4d g = -2 * d_g;

  Eigen::Matrix2d mA = bigM.block<2, 2>(0, 0);
  Eigen::Matrix2d mB = bigM.block<2, 2>(0, 2);
  Eigen::Matrix2d mD = bigM.block<2, 2>(2, 2);

  Eigen::Matrix2d mS;
  Eigen::Matrix2d mSa;

  mS = mD - mB.transpose() * mA.inverse() * mB;

  /* mSa = mS.inv * mS.det; */
  mSa = mS.inverse();
  mSa *= mS.determinant();

  Eigen::Matrix<double, 2, 1> g1 = g.block<2, 1>(0, 0);
  Eigen::Matrix<double, 2, 1> g2 = g.block<2, 1>(2, 0);
  Eigen::Matrix<double, 1, 2> g1t = g1.transpose();
  Eigen::Matrix<double, 1, 2> g2t = g2.transpose();
  Eigen::Matrix2d mAi = mA.inverse();

  Eigen::Matrix<double, 1, 2> m1t = g1t * mAi * mB;
  Eigen::Matrix<double, 2, 1> m1 = m1t.transpose();

  Eigen::Matrix<double, 1, 2> m2t = m1t * mSa;
  Eigen::Matrix<double, 2, 1> m2 = m2t.transpose();

  Eigen::Matrix<double, 1, 2> m3t = g2t * mSa;
  Eigen::Matrix<double, 2, 1> m3 = m3t.transpose();

  double p[3] = {m2t.dot(m2) - 2 * m2t.dot(m3) + m3t.dot(m3),
                 4 * m2t.dot(m1) - 8 * m2t.dot(g2) + 4 * g2t.dot(m3),
                 4 * m1t.dot(m1) - 8 * m1t.dot(g2) + 4 * g2t.dot(g2)};

  double l[3] = {mS.determinant(), 2 * mS(0, 0) + 2 * mS(1, 1), 4};

  /* q = p - l^2        */
  double q[5] = {p[0] - (l[0] * l[0]), p[1] - (2 * l[1] * l[0]),
                 p[2] - (l[1] * l[1] + 2 * l[0] * l[2]), -(2 * l[2] * l[1]),
                 -(l[2] * l[2])};

  double lambda;
  if (!PolyGreatestRealRoot(5, q, &lambda)) return 0;

  Eigen::Matrix<double, 4, 4> W;
  W.setZero();
  W(2, 2) = 1.0;
  W(3, 3) = 1.0;
  Eigen::Matrix<double, 4, 1> x;
  W *= (2 * lambda);

  bigM += W;
  x = bigM.inverse() * g;
  x *= -1.0;

  (*x_out) << x(0, 0), x(1, 0), atan2(x(3, 0), x(2, 0));

  return true;
}

double IcpParams::GpcError(const GpcCorrespondence& co,
                           const Eigen::Vector3d& x) {
  double connections = cos(x[2]);
  double s = sin(x[2]);
  double e[2];
  e[0] = connections * (co.p[0]) - s * (co.p[1]) + x[0] - co.q[0];
  e[1] = s * (co.p[0]) + connections * (co.p[1]) + x[1] - co.q[1];
  double this_error = e[0] * e[0] * co.C[0][0] + 2 * e[0] * e[1] * co.C[0][1] +
                      e[1] * e[1] * co.C[1][1];

  return this_error;
}

double IcpParams::GpcTotalError(const std::vector<GpcCorrespondence>& co,
                                const int n, const Eigen::Vector3d& x) {
  double error = 0;
  for (int i = 0; i < n; i++) {
    if (!co[i].valid) {
      continue;
    }
    error += GpcError(co[i], x);
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
  // const LidarData* laser_ref  = this->laser_ref;
  // const LidarData* laser_sens = this->laser_sens;
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
                                   ? Correspondence::corr_pl
                                   : Correspondence::corr_pp;
  }
}

void IcpParams::KillOutliersDouble() {
  double threshold = 3; /* TODO: add as configurable */

  LidarData* laser_ref = this->laser_ref;
  LidarData* laser_sens = this->laser_sens;

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
  LidarData* laser_ref = this->laser_ref;
  LidarData* laser_sens = this->laser_sens;

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

int IcpParams::TerminationCriterion(const Eigen::Vector3d& delta) {
  double a = delta.squaredNorm();
  double b = fabs(delta[2]);
  return (a < this->epsilon_xy) && (b < this->epsilon_theta);
}

int IcpParams::IcpLoop(const Eigen::Vector3d& q0, Eigen::Vector3d* const x_new,
                       double* const total_error, int* const valid,
                       int* const iterations) {
  if (q0.hasNaN()) {
    LOG(ERROR) << "icp_loop: Initial pose contains nan: " << q0;
    return 0;
  }

  // 新一帧的位置
  LidarData* laser_sens = reinterpret_cast<LidarData*>(this->laser_sens);

  Eigen::Vector3d x_old, delta, delta_old;
  x_old = q0;
  delta_old.setZero();

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
    MathUtils::PoseDiff(x_new->data(), x_old.data(), delta.data());

    /** Checks for oscillations */
    hashes[iteration] = laser_sens->CorrHash();

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

    x_old = *x_new;
    delta_old = delta;

    // egsl_pop_named("icp_loop iteration");
  }

  *iterations = iteration + 1;

  return all_is_okay;
}

void IcpParams::PLIcp(IcpResult* const result) {
  result->valid_ = 0;

  LidarData* laser_ref = reinterpret_cast<LidarData*>(this->laser_ref);
  LidarData* laser_sens = reinterpret_cast<LidarData*>(this->laser_sens);

  if (!laser_ref->ValidLidar() || !laser_sens->ValidLidar()) {
    return;
  }
  laser_ref->InvalidIfOutside(this->min_reading, this->max_reading);
  laser_sens->InvalidIfOutside(this->min_reading, this->max_reading);

  // Origin laser is range and bearing, here compute as Cartesian.
  laser_ref->ComputeCartesian();
  laser_sens->ComputeCartesian();

  Eigen::Vector3d x_new;
  Eigen::Vector3d x_old(this->first_guess[0], this->first_guess[1],
                        this->first_guess[2]);

  double error;
  int iterations;
  int nvalid;

  if (!IcpLoop(x_old, &x_new, &error, &nvalid, &iterations)) {
    // sm_error("icp: ICP failed for some reason. \n");
    result->valid_ = 0;
    result->iterations_ = iterations;
    result->nvalid_ = 0;
  } else {
    /* It was succesfull */

    int restarted = 0;
    double best_error = error;
    Eigen::Vector3d best_x = x_new;

    if (this->restart &&
        (error / nvalid) > (this->restart_threshold_mean_error)) {
      LOG(WARNING) << "Restarting: " << (error / nvalid) << " > "
                   << this->restart_threshold_mean_error;
      restarted = 1;
      double dt = this->restart_dt;
      double dth = this->restart_dtheta;
      // sm_debug("icp_loop: dt = %f dtheta= %f deg\n", dt, rad2deg(dth));

      double perturb[6][3] = {{dt, 0, 0},  {-dt, 0, 0}, {0, dt, 0},
                              {0, -dt, 0}, {0, 0, dth}, {0, 0, -dth}};

      for (int a = 0; a < 6; a++) {
        // sm_debug("-- Restarting with perturbation #%d\n", a);
        IcpParams my_params = *this;
        Eigen::Vector3d start;
        start[0] = x_new[0] + perturb[a][0];
        start[1] = x_new[1] + perturb[a][1];
        start[2] = x_new[2] + perturb[a][2];
        Eigen::Vector3d x_a;
        double my_error;
        int my_valid;
        int my_iterations;
        if (!my_params.IcpLoop(start, &x_a, &my_error, &my_valid,
                               &my_iterations)) {
          // sm_error("Error during restart #%d/%d. \n", a, 6);
          break;
        }
        iterations += my_iterations;

        if (my_error < best_error) {
          // sm_debug("--Perturbation #%d resulted in error %f < %f\n", a,
          //          my_error, best_error);
          best_x = x_a;
          best_error = my_error;
        }
      }
    }

    /* At last, we did it. */
    result->valid_ = 1;
    for (int i = 0; i < 3; ++i) {
      result->x_[i] = best_x[i];
    }
    // sm_debug("icp: final x =  %s  \n", gsl_friendly_pose(best_x));

    if (restarted) {  // recompute correspondences in case of restarts
      laser_sens->ComputeWorldCoords(result->x_);
      // if (this->use_corr_tricks)
      //   find_correspondences_tricks(this);
      // else
      FindCorrespondences();
    }

    if (this->do_compute_covariance) {
      // TODO
      // val cov0_x, dx_dy1, dx_dy2;
      // compute_covariance_exact(laser_ref, laser_sens, best_x, &cov0_x,
      // &dx_dy1,
      //                          &dx_dy2);

      // val cov_x = sc(square(this->sigma), cov0_x);
      // /*			egsl_v2da(cov_x, result->cov_x); */

      // result->cov_x_m = egsl_v2gslm(cov_x);
      // result->dx_dy1_m = egsl_v2gslm(dx_dy1);
      // result->dx_dy2_m = egsl_v2gslm(dx_dy2);
    }

    result->error_ = best_error;
    result->iterations_ = iterations;
    result->nvalid_ = nvalid;
  }
}

int IcpParams::ComputeNextEstimate(const Eigen::Vector3d x_old,
                                   Eigen::Vector3d* const x_new) {
  LidarData* laser_ref = this->laser_ref;
  LidarData* laser_sens = this->laser_sens;

  // struct GpcCorrespondence connections[laser_sens->nrays];
  GpcCorrespondence dummy;
  std::vector<GpcCorrespondence> connections(laser_sens->nrays, dummy);

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

    connections[k].valid = 1;

    if (laser_sens->corr[i].type == Correspondence::corr_pl) {
      connections[k].p[0] = laser_sens->points[i].p[0];
      connections[k].p[1] = laser_sens->points[i].p[1];
      connections[k].q[0] = laser_ref->points[j1].p[0];
      connections[k].q[1] = laser_ref->points[j1].p[1];

      Eigen::Vector3d point_1(laser_ref->points[j1].p[0],
                              laser_ref->points[j1].p[1], 1);
      Eigen::Vector3d point_2(laser_ref->points[j2].p[0],
                              laser_ref->points[j2].p[1], 1);
      connections[k].line = point_1.cross(point_2);
      connections[k].line.normalize();

      double diff[2];
      diff[0] = laser_ref->points[j1].p[0] - laser_ref->points[j2].p[0];
      diff[1] = laser_ref->points[j1].p[1] - laser_ref->points[j2].p[1];
      double one_on_norm = 1 / sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
      // 误差归一化
      double normal[2];
      normal[0] = +diff[1] * one_on_norm;
      normal[1] = -diff[0] * one_on_norm;

      double cos_alpha = normal[0];
      double sin_alpha = normal[1];

      // 2D方向余旋矩阵
      connections[k].C[0][0] = cos_alpha * cos_alpha;
      connections[k].C[1][0] = connections[k].C[0][1] = cos_alpha * sin_alpha;
      connections[k].C[1][1] = sin_alpha * sin_alpha;

    } else {
      connections[k].p[0] = laser_sens->points[i].p[0];
      connections[k].p[1] = laser_sens->points[i].p[1];

      MathUtils::ProjectionOnSegment(
          laser_ref->points[j1].p, laser_ref->points[j2].p,
          laser_sens->points_w[i].p, connections[k].q);

      /* Identity matrix */
      connections[k].C[0][0] = 1;
      connections[k].C[1][0] = 0;
      connections[k].C[0][1] = 0;
      connections[k].C[1][1] = 1;
    }

    k++;
  }

  // Core
  Alpha::TicToc tic;
  SolveOptimization(connections, x_old, x_new);
  std::cout << "SolveOptimization: " << tic.TocMicroseconds() << " ms"
            << std::endl;
  tic.Tic();
  int ok = GpcSolve(k, connections, x_new);
  std::cout << " GpcSolve: " << tic.TocMicroseconds() << " ms" << std::endl;

  if (!ok) {
    LOG(ERROR) << "gpc_solve_valid failed";
    return 0;
  }

  // double old_error = GpcTotalError(connections, k, x_old);
  // double new_error = GpcTotalError(connections, k, *x_new);

  return 1;
}