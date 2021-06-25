#pragma once
#include "math.h"

class MathUtils {
 public:
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
  static double norm_d(const double p[2]) {
    return sqrt(p[0] * p[0] + p[1] * p[1]);
  }
};