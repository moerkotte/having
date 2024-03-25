#ifndef GM_DISTRIBUTION_BETA_HH
#define GM_DISTRIBUTION_BETA_HH
#pragma once

#include "util.hh"
#include "gm_dist_attr_stat.hh"

#include <boost/math/distributions/beta.hpp>

namespace gm {

/*
 *  dist_beta_t
 */

class dist_beta_t {
  public:
    using dist_t = boost::math::beta_distribution<double>;
  public:
    dist_beta_t() : _min(0), _max(0), _dist() {}
    dist_beta_t(const double aMin, const double aMax, const double aMean, const double aVar) 
               : _min(aMin), _max(aMax), 
                 _dist(get_alpha(aMean, aVar, aMin, aMax),
                       get_beta(aMean, aVar, aMin, aMax)) {}
    dist_beta_t(const attr_stat_t& x)
               : _min(x.min()), _max(x.max()),
                 _dist(get_alpha(x.mean(), x.variance(), x.min(), x.max()),
                       get_beta(x.mean(), x.variance(), x.min(), x.max())) {}
  public:
    inline double min() const { return _min; }
    inline double max() const { return _max; }
    inline double alpha() const { return _dist.alpha(); }
    inline double beta() const { return _dist.beta(); }
  public:
    inline double to_01(const double x) const { return ((x - min()) / (max() - min())); }
    inline double to_rg(const double x) const { return (min() + x * (max() - min())); }
  public:
    inline double get_pdf(const double x) const {
                    return pdf(_dist, to_01(x));
                  }
    inline double get_cdf(const double x) const {
                    return cdf(_dist, to_01(x));
                  }
  public:
     // normalize to range [0,1] and denormalize to [lb,ub]
    static inline double mean_01_1(const double aMean, const double aMin, const double aMax) {
                     return ((aMean - aMin) / (aMax - aMin));
                   }
    static inline double var_01_1(const double aVariance, const double aMin, const double aMax) {
                     const double d = aMax - aMin;
                     const double f = double{1} / d;
                     return (f * f * aVariance);
                   }
    // calculate alpha, beta 
    static inline double get_alpha(const double aMean, const double aVar, const double aMin, const double aMax) {
                           return dist_t::find_alpha(mean_01_1(aMean, aMin, aMax), var_01_1(aVar, aMin, aMax));
                         }
    static inline double get_beta(const double aMean, const double aVar, const double aMin, const double aMax) {
                           return dist_t::find_beta(mean_01_1(aMean, aMin, aMax), var_01_1(aVar, aMin, aMax));
                         }
  private:
    double _min;
    double _max;
    dist_t _dist;
};

} // end namespace

#endif

