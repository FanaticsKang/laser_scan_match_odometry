#pragma once

#include "csm/algos.h"

class IcpResult : public sm_result {
 public:
  IcpResult(const sm_result& base) : sm_result(base){};
};