#pragma once
#include "math.h"

class MathUtils {
 public:
  static double deg2rad(double deg) { return deg * (M_PI / 180); }
  static double rad2deg(double rad) { return rad * (180 / M_PI); }
};