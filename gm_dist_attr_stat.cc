#include "gm_dist_attr_stat.hh"
#include <cmath>

namespace gm {

double
attr_stat_t::apply_select(attr_stat_t& aRes, const double aCardR, const double aSelectivity) const {
  const double s = aSelectivity;
  const double lCardRSel = s * aCardR;
  aRes._min  = min();
  aRes._max  = max();
  aRes._nodv = calc_yao(lCardRSel, aCardR, nodv());
  aRes._is_int = is_int();
  aRes._is_ui_dense = is_ui_dense();
  aRes._dist_param = dist_param();
  return lCardRSel;
}

} // end namespace

