#pragma once

#include "csm/algos.h"
#include "icp_result.h"
#include "math_utils.h"
#include "lidar_data.h"

class IcpParams : public sm_params {
 public:
  IcpParams(const sm_params& base);
  void PLIcp(IcpResult* const result);
  int ValidLidar(LidarData* const ld);
};