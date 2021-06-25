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

bool LidarData::ValidRay(const int i) const {
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

unsigned int LidarData::CorrHash() {
  unsigned int hash = 0;
  unsigned int i = 0;

  for (unsigned int i = 0; i < (unsigned)this->nrays; i++) {
    int str =
        this->corr[i].valid ? (this->corr[i].j1 + 1000 * this->corr[i].j2) : -1;

    hash ^= ((i & 1) == 0) ? ((hash << 7) ^ (str) ^ (hash >> 3))
                           : (~((hash << 11) ^ (str) ^ (hash >> 5)));
  }

  return (hash & 0x7FFFFFFF);
}

void LidarData::PossibleInterval(const double* p_i_w,
                                 double max_angular_correction_deg,
                                 double max_linear_correction, int* from,
                                 int* to, int* start_cell) {
  const double angle_res = (this->max_theta - this->min_theta) / this->nrays;

  /* Delta for the angle */
  // max_angular_correction_deg 是两帧之间最大的角度差
  // max_linear_correction 是两帧之间最大的位移差
  const double delta =
      fabs(MathUtils::Deg2Rad(max_angular_correction_deg)) +
      fabs(atan(max_linear_correction / MathUtils::Norm2D(p_i_w)));

  /* Dimension of the cell range */
  const int range = (int)ceil(delta / angle_res);

  /* To be turned into an interval of cells */
  double start_theta = atan2(p_i_w[1], p_i_w[0]);

  /* Make sure that start_theta is in the interval [min_theta,max_theta].
     For example, -1 is not in [0, 2pi] */
  if (start_theta < min_theta) {
    start_theta += 2 * M_PI;
  }
  if (start_theta > max_theta) {
    start_theta -= 2 * M_PI;
  }

  // *start_cell = (start_theta - this->min_theta) / angle_res;
  *start_cell = (int)((start_theta - this->min_theta) /
                      (this->max_theta - this->min_theta) * this->nrays);

  *from = MinMax(0, this->nrays - 1, *start_cell - range);
  *to = MinMax(0, this->nrays - 1, *start_cell + range);

  // LOG(WARNING) << "from: " << *from << " to: " << *to << " delta: " << delta
  //              << " start_theta: " << start_theta
  //              << " min/max theta: " << min_theta << " / " << max_theta
  //              << " range: " << range << " start cell: " << *start_cell
  //              << " max angular deg: " << max_angular_correction_deg
  //              << " max_linear： " << max_linear_correction
  //              << " norm_p: " << MathUtils::Norm2D(p_i_w);
}

int LidarData::NumValidCorrespondences() {
  int num = 0;
  for (int i = 0; i < this->nrays; i++) {
    if (this->corr[i].valid) {
      num++;
    }
  }
  return num;
}
void LidarData::SetNullCorrespondence(const int i) {
  this->corr[i].valid = 0;
  this->corr[i].j1 = -1;
  this->corr[i].j2 = -1;
  this->corr[i].dist2_j1 =
      std::numeric_limits<double>::quiet_NaN();  // GSL_NAN;
}

void LidarData::SetCorrespondence(int i, int j1, int j2) {
  this->corr[i].valid = 1;
  this->corr[i].j1 = j1;
  this->corr[i].j2 = j2;
}

// TODO 改为Eigen
void LidarData::ComputeWorldCoords(const double* pose) {
  const double pose_x = pose[0];
  const double pose_y = pose[1];
  const double pose_theta = pose[2];
  const double cos_theta = cos(pose_theta);
  const double sin_theta = sin(pose_theta);
  const int nrays = this->nrays;

  point2d* points = this->points;
  point2d* points_w = this->points_w;

  for (int i = 0; i < nrays; i++) {
    if (!ValidRay(i)) {
      continue;
    }
    const double& x0 = points[i].p[0];
    const double& y0 = points[i].p[1];

    if (isnan(x0) || isnan(y0)) {
      LOG(ERROR)
          << "ld_compute_world_coords(): I expected that cartesian coords were "
             "already computed: ray"
          << i << x0 << y0;
    }

    points_w[i].p[0] = cos_theta * x0 - sin_theta * y0 + pose_x;
    points_w[i].p[1] = sin_theta * x0 + cos_theta * y0 + pose_y;
    /* polar coordinates */
  }

  for (int i = 0; i < nrays; i++) {
    double x = points_w[i].p[0];
    double y = points_w[i].p[1];
    points_w[i].rho = sqrt(x * x + y * y);
    points_w[i].phi = atan2(y, x);
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

  const double min_fov = MathUtils::Deg2Rad(20.0);
  const double max_fov = 2.01 * M_PI;
  const double fov = this->max_theta - this->min_theta;

  if (fov < min_fov || fov > max_fov) {
    LOG(ERROR) << "Strange FOV: " << fov << " rad : " << MathUtils::Rad2Deg(fov)
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