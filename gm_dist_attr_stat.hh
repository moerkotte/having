#ifndef GM_DISTRIBUTION_ATTR_STAT_HH
#define GM_DISTRIBUTION_ATTR_STAT_HH
#pragma once

#include "util.hh"
#include "gm_dist_param.hh"

namespace gm {

/*
 *  struct attr_stat_t
 *  similar to simple profile
 */

struct attr_stat_t {
  double       _min;
  double       _max;
  uint64_t     _nodv;
  dist_param_t _dist_param;
  bool         _is_int;
  bool         _is_ui_dense;
  attr_stat_t() : _min(0), _max(0), _nodv(0), _dist_param(), _is_int(false), _is_ui_dense(false) {}
  inline double   min() const { return _min; }
  inline double   max() const { return _max; }
  inline uint64_t nodv() const { return _nodv; }
  inline double   nodv_d() const { return ((double) _nodv); }
  inline double   mean() const { return _dist_param.mean(); }
  inline double   std_dev() const { return _dist_param.std_dev(); }
  inline double   variance() const { return _dist_param.variance(); }
  inline double   skew() const { return _dist_param.skew(); }
  inline double   avg_dist()  const { return ((max() - min()) / (nodv_d() - 1)); }
  inline bool     is_int() const { return _is_int; }
  inline bool     is_ui_dense() const { return _is_ui_dense; }
  inline const dist_param_t& dist_param() const { return _dist_param; }

  // generate new attr_stat_t after some selection with predicate p and selectivity aSelectivity
  // has been applied to relation R with |R| = aCardR
  // returns aSelectivity * |R|
  double apply_select(attr_stat_t& aRes, const double aCardR, const double aSelectivity) const;
};


} // end namespace

#endif

