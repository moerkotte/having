#ifndef GM_DISTRIBUTION_UNIFORM_HH
#define GM_DISTRIBUTION_UNIFORM_HH
#pragma once

#include "util.hh"
#include "gm_dist_attr_stat.hh"

#include <boost/math/distributions/uniform.hpp>

namespace gm {

/*
 * class dist_uniform_int_t
 */

class dist_uniform_int_t {
  public:
    using dist_t = boost::math::uniform_distribution<int64_t>;
  public:
    dist_uniform_int_t() : _dist() {}
    dist_uniform_int_t(const double aMin,
                       const double aMax) : _dist(aMin, aMax) {}
    dist_uniform_int_t(const attr_stat_t& aAttrStat) : _dist(aAttrStat.min(), aAttrStat.max()) {}
  public:
    inline double min()  const { return _dist.lower(); }
    inline double max()  const { return _dist.upper(); }
    inline double nodv() const { return (max() - min() + 1); }
  public:
    inline double get_cdf(const double x) const { return cdf(_dist, x); }
    inline double get_pdf(const double x) const { return pdf(_dist, x); }
  public:
    inline double mean() const { return ((min() + max()) / 2); }
    inline double std_dev() const { return std::sqrt(variance()); }
    inline double variance() const { return (((nodv() * nodv()) - 1) / double{12}); }
    inline double skewness() const { return 0; }
    inline double ex_kurtosis() const { return ((double{-6} * (nodv() * nodv() + 1)) / 
                                                (double{5} * (nodv() * nodv() - 1))); }
  private:
    dist_t _dist;
};


/*
 *  class dist_uniform_double_t
 */

class dist_uniform_double_t {
  public:
    using dist_t = boost::math::uniform_distribution<double>;
  public:
    dist_uniform_double_t() : _dist() {}
    dist_uniform_double_t(const double aMin,
                          const double aMax) : _dist(aMin, aMax) {}
    dist_uniform_double_t(const attr_stat_t& aAttrStat) : _dist(aAttrStat.min(), aAttrStat.max()) {}
  public:
    inline double min() const { return _dist.lower(); }
    inline double max() const { return _dist.upper(); }
  public:
    inline double get_cdf(const double x) const { return cdf(_dist, x); }
    inline double get_pdf(const double x) const { return pdf(_dist, x); }
  public:
    inline double mean() const { return ((min() + max()) / 2); }
    inline double std_dev() const { return std::sqrt(variance()); }
    inline double variance() const { return (((max() - min()) * (max() - min())) / 12); }
    inline double skewness() const { return 0; }
    inline double ex_kurtosis() const { return (double{-6} / double{5}); }
  private:
    dist_t _dist;
};



} // end namespace

#endif
