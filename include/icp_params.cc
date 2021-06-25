#include <iostream>

#include <glog/logging.h>

#include "icp_params.h"
#include "math_utils.h"

IcpParams::IcpParams(const sm_params& base) : sm_params(base) {}

void IcpParams::PLIcp(IcpResult* const result) {
  result->valid = 0;

  LidarData laser_ref(*(this->laser_ref));
  LidarData laser_sens(*(this->laser_sens));

  if (!laser_ref.ValidLidar() || !laser_sens.ValidLidar()) {
    return;
  }

  /** Mark as invalid the rays outside of (min_reading, max_reading] */
  //   ld_invalid_if_outside(laser_ref, params->min_reading,
  //   params->max_reading); ld_invalid_if_outside(laser_sens,
  //   params->min_reading, params->max_reading);
}