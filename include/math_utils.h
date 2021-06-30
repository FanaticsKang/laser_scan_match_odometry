#pragma once
#include <math.h>
#include <iostream>

#include <eigen3/Eigen/Core>

class MathUtils {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  static double Deg2Rad(double deg) { return deg * (M_PI / 180); }
  static double Rad2Deg(double rad) { return rad * (180 / M_PI); }
  static double AngleDiff(const double a, const double b) {
    double t = a - b;
    while (t < -M_PI) {
      t += 2 * M_PI;
    }
    while (t > M_PI) {
      t -= 2 * M_PI;
    }
    return t;
  }
  static double Norm2D(const double p[2]) {
    return sqrt(p[0] * p[0] + p[1] * p[1]);
  }

  // TODO These function will be rebuild for Eigen
  static bool AnyNan(const double* d, const int n) {
    for (int i = 0; i < n; ++i) {
      if (isnan(d[i])) {
        return true;
      }
    }
    return false;
  }

  // copy the arrary
  static void Copy(const double* from, int n, double* const to) {
    for (int i = 0; i < n; ++i) {
      to[i] = from[i];
    }
  }

  static double DistanceSquared(const double a[2], const double b[2]) {
    double x = a[0] - b[0];
    double y = a[1] - b[1];
    return x * x + y * y;
  }

  static double ProjectionOnLine(const double a[2], const double b[2],
                                 const double p[2], double res[2]) {
    double t0 = a[0] - b[0];
    double t1 = a[1] - b[1];
    // 1/线段长
    double one_on_r = 1 / sqrt(t0 * t0 + t1 * t1);

    /* normal */
    double nx = t1 * one_on_r;
    double ny = -t0 * one_on_r;
    double c = nx, s = ny;
    double rho = c * a[0] + s * a[1];

    res[0] = c * rho + s * s * p[0] - c * s * p[1];
    res[1] = s * rho - c * s * p[0] + c * c * p[1];

    // Eigen::Vector2d point_a(a[0], a[1]);
    // Eigen::Vector2d point_b(b[0], b[1]);
    // Eigen::Vector2d point_p(p[0], p[1]);

    // Eigen::Vector2d ab = point_b - point_a;
    // Eigen::Vector2d ap = point_p - point_a;

    // const double ap_ab = ap.dot(ab);
    // const double r2 = ab.squaredNorm();
    // Eigen::Vector2d tmp = ab * ap_ab / r2;
    // Eigen::Vector2d test = point_a + tmp;

    // std::cout << "res[0]: " << res[0] << ", res[1]: " << res[1] << std::endl;
    // std::cout << "test: " << test.transpose() << std::endl;

    return fabs(rho - (c * p[0] + s * p[1]));
  }

  static double DistToSegment(const double a[2], const double b[2],
                              const double x[2]) {
    double proj[2];
    double distance = ProjectionOnLine(a, b, x, proj);
    if ((proj[0] - a[0]) * (proj[0] - b[0]) +
            (proj[1] - a[1]) * (proj[1] - b[1]) <
        0) {
      /* the projection is inside the segment */
      return distance;
    } else {
      return sqrt((std::min)(DistanceSquared(a, x), DistanceSquared(b, x)));
    }
  }

  static void ProjectionOnSegment(const double a[2], const double b[2],
                                  const double x[2], double proj[2]) {
    ProjectionOnLine(a, b, x, proj);
    if ((proj[0] - a[0]) * (proj[0] - b[0]) +
            (proj[1] - a[1]) * (proj[1] - b[1]) <
        0) {
      /* the projection is inside the segment */
    } else if (DistanceSquared(a, x) < DistanceSquared(b, x)) {
      Copy(a, 2, proj);
    } else {
      Copy(b, 2, proj);
    }
  }

  static void Ominus(const double x[3], double res[3]) {
    double c = cos(x[2]);
    double s = sin(x[2]);
    res[0] = -c * x[0] - s * x[1];
    res[1] = s * x[0] - c * x[1];
    res[2] = -x[2];
  }

  /** safe if res == x1 */
  static void Oplus(const double x1[3], const double x2[3], double res[3]) {
    double c = cos(x1[2]);
    double s = sin(x1[2]);
    double x = x1[0] + c * x2[0] - s * x2[1];
    double y = x1[1] + s * x2[0] + c * x2[1];
    double theta = x1[2] + x2[2];
    res[0] = x;
    res[1] = y;
    res[2] = theta;
  }

  static void PoseDiff(const double pose2[3], const double pose1[3],
                       double res[3]) {
    double temp[3];
    Ominus(pose1, temp);
    Oplus(temp, pose2, res);

    while (res[2] > +M_PI) res[2] -= 2 * M_PI;
    while (res[2] < -M_PI) res[2] += 2 * M_PI;
  }
};