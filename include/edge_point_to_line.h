#pragma once

#include "g2o/core/base_unary_edge.h"
#include "g2o/types/slam2d/edge_se2.h"

class EdgePointToLine : public g2o::BaseUnaryEdge<1, double, g2o::VertexSE2> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  EdgePointToLine(const Eigen::Vector2d& point, const Eigen::Vector3d& line)
      : point_r_(point), line_w_(line) {}

  bool read(std::istream& is){};
  bool write(std::ostream& os) const {};

  void computeError() {
    const g2o::VertexSE2* v0 = static_cast<const g2o::VertexSE2*>(_vertices[0]);
    Eigen::Isometry2d T_wr = v0->estimate().toIsometry();

    Eigen::Vector2d point_w = T_wr * point_r_;
    // 齐次坐标系
    Eigen::Vector3d _point_w(point_w.x(), point_w.y(), 1.0);
    _error = line_w_.transpose() * _point_w;
  }

  virtual void linearizeOplus() {
    Eigen::Matrix3d dp_dt;
    dp_dt.setZero();
    dp_dt(0, 0) = 1;
    dp_dt(1, 1) = 1;
    const double p_x = point_r_[0];
    const double p_y = point_r_[1];
    g2o::VertexSE2* vi = static_cast<g2o::VertexSE2*>(_vertices[0]);
    g2o::SE2 se2 = vi->estimate();
    Eigen::Matrix2d rotation = se2.rotation().toRotationMatrix();
    // const double theta = se2.rotation().smallestPositiveAngle();
    // std::cout << "theta: " << theta << ", cos: " << cos(theta)
    //           << ", sin: " << sin(theta) << ", rotation: \n"
    //           << rotation << std::endl;
    const double cos = rotation(0, 0);
    const double sin = rotation(0, 1);
    dp_dt(0, 2) = p_x * sin - p_y * cos;
    dp_dt(1, 2) = p_x * cos + p_y * sin;

    Eigen::Vector3d jacobian = line_w_.transpose() * dp_dt;
    _jacobianOplusXi(0, 0) = jacobian(0);
    _jacobianOplusXi(0, 1) = jacobian(1);
    _jacobianOplusXi(0, 2) = jacobian(2);
  }

 private:
  // point in robot coordination
  Eigen::Vector2d point_r_;
  // line in world coordination
  Eigen::Vector3d line_w_;
};