#ifndef GM_DISTRIBUTION_NORMAL_HH
#define GM_DISTRIBUTION_NORMAL_HH
#pragma once

#include "util.hh"
#include "gm_dist_attr_stat.hh"
#include <boost/math/distributions/normal.hpp>

/*
 * struct dist_param_t
 *  used to store distribution parameters mean, std_dev, skew
 * attr_stat_t
 *  used to store dist_param_t for some attribute as well as
 *  other number of the simple profile as well as extended simple profile
 *
 * Wrapper for boost distributions 
 * to make them all instantiable from attr_stat_t
 */

namespace gm {

/*
 *  dist_normal_t
 */

class dist_normal_t {
  public:
    using dist_t = boost::math::normal_distribution<double>;
  public:
    dist_normal_t() : _min(0), _max(0), _dist() {}
    dist_normal_t(const double aMin, const double aMax, const double aMean, const double aStdDev)
                 : _min(aMin), _max(aMax), _dist(aMean, aStdDev) {}
  public:
    inline double min() const { return _min; }
    inline double max() const { return _max; }
  public:
    inline double get_cdf(const double x) const { return cdf(_dist, x); }
    inline double get_pdf(const double x) const { return pdf(_dist, x); }
  public:
    inline double mean() const { return _dist.mean(); }
    inline double std_dev() const { return _dist.standard_deviation(); }
    inline double variance() const { return (std_dev() * std_dev()); }
    inline double skewness() const { return 0; }
    inline double ex_kurtosis() const { return 0; }
  private:
    double _min;
    double _max;
    dist_t _dist;
};

} // end namespace

#endif

