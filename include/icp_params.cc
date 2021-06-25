#include <iostream>

#include <glog/logging.h>

#include "icp_params.h"
#include "math_utils.h"

IcpParams::IcpParams(const sm_params& base) : sm_params(base) {}

int IcpParams::ValidLidar(LidarData* const lidar) {
  if (lidar == nullptr) {
    LOG(ERROR) << "NULL pointer given as lidar.";
    return 0;
  }

  const int min_nrays = 10;
  const int max_nrays = 10000;

  if (lidar->nrays < min_nrays || lidar->nrays > max_nrays) {
    LOG(ERROR) << "Invalid number of rays: " << lidar->nrays;
    return 0;
  }

  if (isnan(lidar->min_theta) || isnan(lidar->max_theta)) {
    LOG(ERROR) << "Invalid min / max theta: min_theta = " << lidar->min_theta
               << " max_theta = " << lidar->max_theta;
    return 0;
  }

  const double min_fov = MathUtils::deg2rad(20.0);
  const double max_fov = 2.01 * M_PI;
  const double fov = lidar->max_theta - lidar->min_theta;

  if (fov < min_fov || fov > max_fov) {
    LOG(ERROR) << "Strange FOV: " << fov << " rad = " << MathUtils::rad2deg(fov)
               << " deg";
    return 0;
  }

  if (fabs(lidar->min_theta - lidar->theta[0]) > 1e-8) {
    LOG(ERROR) << "Min_theta " << lidar->min_theta << " should be theta[0] "
               << lidar->theta[0];
    return 0;
  }

  if (fabs(lidar->max_theta - lidar->theta[lidar->nrays - 1]) > 1e-8) {
    LOG(ERROR) << "Max_theta " << lidar->max_theta << " should be theta[end] "
               << lidar->theta[lidar->nrays - 1];
    return 0;
  }

  /* Check that there are valid rays */
  const double min_reading = 0;
  const double max_reading = 100;

  for (int i = 0; i < lidar->nrays; i++) {
    double th = lidar->theta[i];
    if (lidar->valid[i]) {
      double r = lidar->readings[i];
      if (isnan(r) || isnan(th)) {
        LOG(ERROR) << "Ray[" << i << "]: r = " << r << "  theta = " << th
                   << " but valid is " << lidar->valid[i];
        return 0;
      }
      if (!(min_reading < r && r < max_reading)) {
        LOG(ERROR) << "Ray[" << i << "]: r = " << r << " is not in interval("
                   << min_reading << ", " << max_reading << ")";
        return 0;
      }
    } else {
      /* ray not valid, but checking theta anyway */
      if (isnan(th)) {
        LOG(ERROR) << "Ray[" << i << "]: valid = " << lidar->valid[0]
                   << " but theta = " << th;
        return 0;
      }

      if (lidar->cluster[i] != -1) {
        LOG(ERROR) << "Invalid ray " << i << " has cluster "
                   << lidar->cluster[i];
        return 0;
      }
    }
    if (lidar->cluster[i] < -1) {
      LOG(ERROR) << "Ray[" << i << "]: Invalid cluster value "
                 << lidar->cluster[i];
      return 0;
    }

    if (!isnan(lidar->readings_sigma[i]) && lidar->readings_sigma[i] < 0) {
      LOG(ERROR) << "Ray[" << i << "]: has invalid readings_sigma "
                 << lidar->readings_sigma[i];
      return 0;
    }
  }
  /* Checks that there is at least 10% valid rays */
  const int num_valid = lidar->CountEqual(lidar->valid, lidar->nrays, 1);
  const int num_invalid = lidar->CountEqual(lidar->valid, lidar->nrays, 0);

  if (num_valid < lidar->nrays * 0.10) {
    LOG(ERROR) << "Valid: " << num_valid << "/" << lidar->nrays
               << " invalid: " << num_invalid;
    return 0;
  }

  return 1;
}

void IcpParams::PLIcp(IcpResult* const result) {
  result->valid = 0;

  LidarData laser_ref(*(this->laser_ref));
  LidarData laser_sens(*(this->laser_sens));

  if (!ValidLidar(&laser_ref) || !ValidLidar(&laser_sens)) {
    return;
  }

  /** Mark as invalid the rays outside of (min_reading, max_reading] */
  //   ld_invalid_if_outside(laser_ref, params->min_reading,
  //   params->max_reading); ld_invalid_if_outside(laser_sens,
  //   params->min_reading, params->max_reading);
}