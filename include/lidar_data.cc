#include <glog/logging.h>
#include <limits>

#include "lidar_data.h"

#include "math_utils.h"

int LidarData::CountEqual(const int* v, int n, int value) {
  int num = 0;
  for (int i = 0; i < n; ++i) {
    if (value == v[i]) {
      ++num;
    }
  }
  return num;
};

bool LidarData::ValidRay(const int i) {
  return (i >= 0) && (i < this->nrays) && (this->valid[i]);
}

void LidarData::InvalidIfOutside(const double min_reading,
                                 const double max_reading) {
  for (int i = 0; i < this->nrays; ++i) {
    if (!ValidRay(i)) {
      continue;
    }
    const double r = this->readings[i];
    if (r <= min_reading || r > max_reading) {
      this->valid[i] = 0;
    }
  }
}

void LidarData::ComputeCartesian() {
  for (int i = 0; i < this->nrays; i++) {
    // if(!ld_valid_ray(this,i)) continue;
    double x = cos(this->theta[i]) * this->readings[i];
    double y = sin(this->theta[i]) * this->readings[i];

    this->points[i].p[0] = x;
    this->points[i].p[1] = y;
    // TODO
    this->points[i].rho = std::numeric_limits<double>::quiet_NaN();
    this->points[i].phi = std::numeric_limits<double>::quiet_NaN();
  }
}

bool LidarData::ValidLidar() {
  const int min_nrays = 10;
  const int max_nrays = 10000;

  if (this->nrays < min_nrays || this->nrays > max_nrays) {
    LOG(ERROR) << "Invalid number of rays: " << this->nrays;
    return false;
  }

  if (isnan(this->min_theta) || isnan(this->max_theta)) {
    LOG(ERROR) << "Invalid min / max theta: min_theta = " << this->min_theta
               << " max_theta = " << this->max_theta;
    return false;
  }

  const double min_fov = MathUtils::deg2rad(20.0);
  const double max_fov = 2.01 * M_PI;
  const double fov = this->max_theta - this->min_theta;

  if (fov < min_fov || fov > max_fov) {
    LOG(ERROR) << "Strange FOV: " << fov << " rad : " << MathUtils::rad2deg(fov)
               << " deg";
    return false;
  }

  if (fabs(this->min_theta - this->theta[0]) > 1e-8) {
    LOG(ERROR) << "Min_theta " << this->min_theta << " should be theta[0] "
               << this->theta[0];
    return false;
  }

  if (fabs(this->max_theta - this->theta[this->nrays - 1]) > 1e-8) {
    LOG(ERROR) << "Max_theta " << this->max_theta << " should be theta[end] "
               << this->theta[this->nrays - 1];
    return false;
  }

  /* Check that there are valid rays */
  const double min_reading = 0;
  const double max_reading = 100;

  for (int i = 0; i < this->nrays; i++) {
    double th = this->theta[i];
    if (this->valid[i]) {
      double r = this->readings[i];
      if (isnan(r) || isnan(th)) {
        LOG(ERROR) << "Ray[" << i << "]: r = " << r << "  theta = " << th
                   << " but valid is " << this->valid[i];
        return false;
      }
      if (!(min_reading < r && r < max_reading)) {
        LOG(ERROR) << "Ray[" << i << "]: r = " << r << " is not in interval("
                   << min_reading << ", " << max_reading << ")";
        return false;
      }
    } else {
      /* ray not valid, but checking theta anyway */
      if (isnan(th)) {
        LOG(ERROR) << "Ray[" << i << "]: valid = " << this->valid[0]
                   << " but theta = " << th;
        return false;
      }

      if (this->cluster[i] != -1) {
        LOG(ERROR) << "Invalid ray " << i << " has cluster "
                   << this->cluster[i];
        return false;
      }
    }
    if (this->cluster[i] < -1) {
      LOG(ERROR) << "Ray[" << i << "]: Invalid cluster value "
                 << this->cluster[i];
      return false;
    }

    if (!isnan(this->readings_sigma[i]) && this->readings_sigma[i] < 0) {
      LOG(ERROR) << "Ray[" << i << "]: has invalid readings_sigma "
                 << this->readings_sigma[i];
      return false;
    }
  }
  /* Checks that there is at least 10% valid rays */
  const int num_valid = this->CountEqual(this->valid, this->nrays, 1);
  const int num_invalid = this->CountEqual(this->valid, this->nrays, 0);

  if (num_valid < this->nrays * 0.10) {
    LOG(ERROR) << "Valid: " << num_valid << "/" << this->nrays
               << " invalid: " << num_invalid;
    return false;
  }

  return true;
}