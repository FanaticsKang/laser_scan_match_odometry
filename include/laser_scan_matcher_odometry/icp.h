#pragma once

#include "csm/algos.h"
class IcpParams : public sm_params {
 public:
  IcpParams(const sm_params& base) : sm_params(base){};
};

class IcpResult : public sm_result {
 public:
  IcpResult(const sm_result& base) : sm_result(base){};
};