#include <iostream>

#include <glog/logging.h>
#include <eigen3/Eigen/Core>

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
  laser_ref.InvalidIfOutside(this->min_reading, this->max_reading);
  laser_sens.InvalidIfOutside(this->min_reading, this->max_reading);

  // TODO
  // if(params->use_corr_tricks || params->debug_verify_tricks)
  // 	ld_create_jump_tables(laser_ref);

  laser_ref.ComputeCartesian();
  laser_sens.ComputeCartesian();

  // TODO
  // if (params->do_alpha_test) {
  //   ld_simple_clustering(laser_ref, params->clustering_threshold);
  //   ld_compute_orientation(laser_ref, params->orientation_neighbourhood,
  //                          params->sigma);
  //   ld_simple_clustering(laser_sens, params->clustering_threshold);
  //   ld_compute_orientation(laser_sens, params->orientation_neighbourhood,
  //                          params->sigma);
  // }

  Eigen::Vector3d x_new;
  // 	gsl_vector * x_old = vector_from_array(3, params->first_guess);
  Eigen::Vector3d x_old(this->first_guess[0], this->first_guess[1],
                        this->first_guess[2]);
}