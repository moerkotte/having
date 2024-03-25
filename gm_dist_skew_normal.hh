#ifndef GM_DISTRIBUTION_SKEW_NORMAL_HH
#define GM_DISTRIBUTION_SKEW_NORMAL_HH
#pragma once

#include "util.hh"
#include "gm_dist_attr_stat.hh"

#include <boost/math/distributions/skew_normal.hpp>

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
 *  dist_skew_normal_t
 */

class dist_skew_normal_t {
  public:
    using dist_t = boost::math::skew_normal_distribution<double>;
  public:
    dist_skew_normal_t() : _min(0), _max(0), _dist() {}
    dist_skew_normal_t(const double aMin, const double aMax, 
                       const double aLocation, const double aScale, const double aShape)
                      : _min(aMin), _max(aMax), _dist(aLocation, aScale, aShape) {}
    dist_skew_normal_t(const double aMin, const double aMax, const dist_param_t& aDistParam)
                      : _min(aMin), _max(aMax), _dist(get_dist(aDistParam)) {}
    dist_skew_normal_t(const attr_stat_t& aAttrStat)
                      : _min(aAttrStat.min()), _max(aAttrStat.max()),
                        _dist(get_dist(aAttrStat.dist_param())) {}
  public:
    inline double min() const { return _min; }
    inline double max() const { return _max; }
  public:         
    inline double get_cdf(const double x) const { return cdf(_dist, x); }
    inline double get_pdf(const double x) const { return pdf(_dist, x); }
  public:
    /*
    inline double mean() const { return _dist.mean(); }
    inline double std_dev() const { return _dist.standard_deviation(); }
    inline double variance() const { return _dist.variance(); }
    inline double skewness() const { return _dist.skewness(); }
    inline double ex_kurtosis() const { return _dist.excess_kurtosis(); }
    */
    double mean() const;
    double std_dev() const;
    double variance() const;
    double skewness() const;
    double ex_kurtosis() const;
    double delta() const;
  public:
    static dist_t get_dist(const dist_param_t& x);
  private:
    double _min;
    double _max;
    dist_t _dist;
};

} // end namespace

#endif

