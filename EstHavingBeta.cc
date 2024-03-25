#include "EstHavingBeta.hh"


EstHavingBeta::EstHavingBeta()
              : _card_R(0), 
                _attr_stat_A(), 
                _attr_stat_B(), 
                _attr_stat_Cnt(),
                _dist_Cnt(),
                _dist_B(), 
                _dist_param_sum(),
                _dist_sum() {}
  
EstHavingBeta::EstHavingBeta(const uint64_t     aCardR,
                             const attr_stat_t& aAttrStatA,
                             const attr_stat_t& aAttrStatB,
                             const attr_stat_t& aAttrStatCnt)
              : _card_R(aCardR),
                _attr_stat_A(aAttrStatA), 
                _attr_stat_B(aAttrStatB), 
                _attr_stat_Cnt(aAttrStatCnt),
                _dist_Cnt(), 
                _dist_B(), 
                _dist_param_sum(),
                _dist_sum() {
  _dist_Cnt = get_dist_Cnt();
  _dist_B   = get_dist_B();
  _dist_sum = get_dist_sum();
}

EstHavingBeta::dist_t
EstHavingBeta::get_dist_Cnt() {
  return dist_t(_attr_stat_Cnt);
}

EstHavingBeta::dist_t
EstHavingBeta::get_dist_B() {
  return dist_t(_attr_stat_B);
}

EstHavingBeta::dist_t
EstHavingBeta::get_dist_sum() {
  _dist_param_sum = dist_param_Cnt() * dist_param_B();
  return dist_t(min_Cnt() * min_B(), max_Cnt() * max_B(), _dist_param_sum.mean(), _dist_param_sum.variance());
}

double
EstHavingBeta::estimate_cnt_eq(const uint aCnt) const {
  const double x = (double) aCnt;
  double l = 0, u = 0;
  if((min_Cnt() > aCnt) || (max_Cnt() < aCnt)) {
    return 0;
  }
  if(1 == aCnt) {
    l = _dist_Cnt.get_cdf(1.0);
    u = _dist_Cnt.get_cdf(1.5);
  } else
  if(max_Cnt() == aCnt) {
    l = _dist_Cnt.get_cdf(x - 0.5);
    u = _dist_Cnt.get_cdf(x);
  } else {
    l = _dist_Cnt.get_cdf(x - 0.5);
    u = _dist_Cnt.get_cdf(x + 0.5);
  }
  return ((double) nodv_A() * (u - l));
}

double
EstHavingBeta::estimate_cnt_between(const uint aLbCnt, const uint aUbCnt) const {
  double l = 0, u = 0;
  const double lLbCnt = std::max<double>(aLbCnt, min_Cnt());
  const double lUbCnt = std::min<double>(aUbCnt, max_Cnt());
  if(lLbCnt > lUbCnt) { return 0; }
  if(lLbCnt == lUbCnt) { return estimate_cnt_eq(lLbCnt); }

  if((min_Cnt() + 0.5) >= lLbCnt) {
    l = _dist_Cnt.get_cdf(min_Cnt());
  } else {
    l = _dist_Cnt.get_cdf(lLbCnt - 0.5);
  }

  if((max_Cnt() + 0.5) <= lUbCnt) {
    u = _dist_Cnt.get_cdf((double) lUbCnt);
  } else {
    u = _dist_Cnt.get_cdf((double) lUbCnt + 0.5);
  }
  return ((double) nodv_A() * (u - l));
}

double
EstHavingBeta::estimate_cnt_eq_sel(const uint aCnt, const double aSelectivity) const {
  return 1;
}

double
EstHavingBeta::estimate_cnt_between_sel(const uint   aLbCnt, 
                                        const uint   aUbCnt, 
                                        const double aSelectivity) const {
  return 1;
}

double
EstHavingBeta::estimate_sum_eq(const double aSumB) const {
  double l = 0, u = 0;
  l = _dist_sum.get_cdf(std::max(min_B(), aSumB - 0.5));
  u = _dist_sum.get_cdf(std::min(max_B() * max_Cnt(), aSumB + 0.5));
  return std::max<double>(1, ((double) nodv_A() * (u - l)));
}

double
EstHavingBeta::estimate_sum_eq_sel(const double aSumB, const double aSelectivity) const {
  return 1;
}

double
EstHavingBeta::estimate_sum_eq_cnt(const double aSumB, const uint aCntLb, const uint aCntUb) const {
  return 1;
}

double
EstHavingBeta::estimate_sum_between(const double aLbSumB, const double aUbSumB) const {
  const double l = _dist_sum.get_cdf(std::max(min_B(), aLbSumB - 0.5));
  const double u = _dist_sum.get_cdf(std::min(max_B() * max_Cnt(), aUbSumB + 0.5));
  return ((double) nodv_A() * (u - l));
}

double
EstHavingBeta::estimate_sum_between_sel(const double aLbSumB, 
                                        const double aUbSumB, 
                                        const double aSelectivity) const {
  return 1;
}

double
EstHavingBeta::estimate_sum_between_and_cnt_between(const double aLbSumB, 
                                                    const double aUbSumB,
                                                    const uint   aLbCnt,  
                                                        const uint   aUbCnt) const {
  return 1;
}



