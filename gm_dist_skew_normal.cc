#include "gm_dist_skew_normal.hh"
#include <cmath>

namespace gm {

dist_skew_normal_t::dist_t
dist_skew_normal_t::get_dist(const dist_param_t& x) {
  double l_xi = 0;
  double l_omega = 0;
  double l_alpha = 0;
  x.get_param_skew_norm_b(l_xi, l_omega, l_alpha);
  return dist_t(l_xi, l_omega, l_alpha);
}


double
dist_skew_normal_t::mean() const {
  const double l = _dist.location();
  const double s = _dist.scale();
  // const double a = _dist.shape();
  const double d = delta();
  const double f = std::sqrt(double{2} / M_PI);
  return (l + s * d * f);
}

double
dist_skew_normal_t::std_dev() const {
  return std::sqrt(variance());
}

double
dist_skew_normal_t::variance() const {
  // const double l = _dist.location();
  const double s = _dist.scale();
  // const double a = _dist.shape();
  const double d = delta();
  return (s * s * (1 - (2 * d * d / M_PI)));
}

double
dist_skew_normal_t::skewness() const {
  // const double l = _dist.location();
  // const double s = _dist.scale();
  // const double a = _dist.shape();
  const double d = delta();
  const double f = std::sqrt(double{2} / M_PI);
  return ( ((double{4} - M_PI) / double{2}) *
           (((d*f)*(d*f)*(d*f))  / std::pow( (1 - 2 * d * d / M_PI), double{3} / double{2})));
}

double
dist_skew_normal_t::ex_kurtosis() const {
  // const double l = _dist.location();
  // const double s = _dist.scale();
  // const double a = _dist.shape();
  const double d = delta();
  const double f = std::sqrt(double{2} / M_PI);
  const double x = (d * f) * (d * f) * (d * f) * (d * f); // (df)^4
  const double y = 1 - 2 * d * d / M_PI;
  return (2 * (M_PI - 3) * (x / (y * y)));
}

double
dist_skew_normal_t::delta() const {
  // const double l = _dist.location();
  // const double s = _dist.scale();
  const double a = _dist.shape();
  return (a / std::sqrt(1 + a * a));
}



} // end namespace

