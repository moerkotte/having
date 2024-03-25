#ifndef DISTRIBUTION_PARAMETER_HH
#define DISTRIBUTION_PARAMETER_HH
#pragma once

#include "util.hh"
#include <optional>

/*
 * struct dist_param_t
 *  used to store distribution parameters mean, std_dev, skew
 */

struct sexpr_t;

namespace gm {

/*
 *  struct dist_param_t
 */

struct dist_param_t {
  using double_ot = std::optional<double>;
  double    _mean;     // mean
  double    _std_dev;  // standard deviation
  double    _skewness; // skewness
  double_ot _kurtosis; // kurtosis
  dist_param_t() : _mean(0), _std_dev(0), _skewness(0), _kurtosis() {}
  dist_param_t(const double aMean,
               const double aStdDev,
               const double aSkew)
              : _mean(aMean), _std_dev(aStdDev), _skewness(aSkew), _kurtosis() {}
  dist_param_t(const double aMean,
               const double aStdDev,
               const double aSkew,
               const double aKurtosis)
              : _mean(aMean), _std_dev(aStdDev), _skewness(aSkew), _kurtosis(aKurtosis) {}

  inline double    mean()     const { return _mean; }
  inline double    std_dev()  const { return _std_dev; }
  inline double    variance() const { return (std_dev() * std_dev()); }
  inline double    var()      const { return variance(); }
  inline double    skew()     const { return _skewness; }
  inline double    skewness() const { return _skewness; }
  inline double    moment_2() const { return (std_dev() * std_dev() + mean() * mean()); }
         double    moment_3() const;
  inline bool      has_kurtosis() const { return _kurtosis.has_value(); }
  inline double_ot kurtosis() const { return _kurtosis; }

  inline double&    mean_nc()     { return _mean; }
  inline double&    std_dev_nc()  { return _std_dev; }
  inline double&    skew_nc()     { return _skewness; }
  inline double&    skewness_nc() { return _skewness; }
  inline double_ot& kurtosis_nc() { return _kurtosis; }


  // helper for beta-distribution with finite support
  // normalized to range [0,1] and denormalize to [lb,ub]
  inline double mean_01_1(const double aMin, const double aMax) const {
                   return ((mean() - aMin) / (aMax - aMin));
                 }
  inline double var_01_1(const double aMin, const double aMax) const {
                   const double d = aMax - aMin;
                   const double f = double{1} / d;
                   return (f * f * variance());
                 }

  // get parameters for skew-normal distribution 
  bool get_param_skew_norm_b(double& aXi, double& aEta, double& aLambda) const;
};

dist_param_t
operator*(const dist_param_t& x, const dist_param_t& y);

/*
 *  calculate yao function
 */

// k: no qual item, n: tot num item, m: num buck
double calc_yao(const double k, const double n, const double m);

} // end namespace

#endif

