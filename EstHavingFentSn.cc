#include "EstHavingFentSn.hh"


EstHavingFentSn::EstHavingFentSn()
                : _card_R(0), 
                  _attr_stat_A(), 
                  _attr_stat_B(), 
                  _dist_param_Cnt(),
                  _dist_Cnt(), 
                  _dist_B(), 
                  _dist_param_sum(),
                  _dist_sum() {}

EstHavingFentSn::EstHavingFentSn(const uint64_t     aCardR,
                                 const attr_stat_t& aAttrStatA,
                                 const attr_stat_t& aAttrStatB,
                                 const attr_stat_t& aAttrStatCnt)
                : _card_R(aCardR),
                  _attr_stat_A(aAttrStatA), 
                  _attr_stat_B(aAttrStatB), 
                  _dist_param_Cnt(),
                  _dist_Cnt(), 
                  _dist_B(), 
                  _dist_param_sum(),
                  _dist_sum() {
  _dist_Cnt = get_dist_Cnt();
  _dist_B   = get_dist_B();
  _dist_sum = get_dist_sum();
}

EstHavingFentSn::dist_t
EstHavingFentSn::get_dist_Cnt() {
  // 1. set _dist_param_Cnt
  const double lCardR = ((double) card_R());
  const double D      = ((double) nodv_A());
  const double p      = double{1} / D;
  dist_param_t lDp;
  _dist_param_Cnt._mean     =  lCardR / D;
  _dist_param_Cnt._std_dev  = std::sqrt(lCardR * p * (1 - p));
  _dist_param_Cnt._skewness = skew_A();
  return dist_t(0, D, _dist_param_Cnt);
}

EstHavingFentSn::dist_t
EstHavingFentSn::get_dist_B() {
  return dist_t(min_B(), max_B(), dist_param_B());
}

EstHavingFentSn::dist_t
EstHavingFentSn::get_dist_sum() {
  _dist_param_sum = dist_param_B() * _dist_param_Cnt;
  return dist_t(1 * min_B(), (((double) card_R()) / ((double) nodv_A())) * max_B(), _dist_param_sum);
}

double
EstHavingFentSn::estimate_cnt_eq(const uint aCnt) const {
  const double x = (double) aCnt;
  double l = 0, u = 0;
  if(1 == aCnt) {
    l = _dist_Cnt.get_cdf(1.0);
    u = _dist_Cnt.get_cdf(1.5);
  } else {
    l = _dist_Cnt.get_cdf(x - 0.5);
    u = _dist_Cnt.get_cdf(x + 0.5);
  }
  return ((double) nodv_A() * (u - l));
}

double
EstHavingFentSn::estimate_cnt_between(const uint aLbCnt, const uint aUbCnt) const {
  double l = 0, u = 0;
  if(1 >= aLbCnt) {
    l = _dist_Cnt.get_cdf(1.0);
    u = _dist_Cnt.get_cdf((double) aUbCnt + 0.5);
  } else {
    l = _dist_Cnt.get_cdf((double) aLbCnt - 0.5);
    u = _dist_Cnt.get_cdf((double) aUbCnt + 0.5);
  }
  return ((double) nodv_A() * (u - l));
}

double
EstHavingFentSn::estimate_cnt_eq_sel(const uint aCnt, const double aSelectivity) const {
  return 1;
}

double
EstHavingFentSn::estimate_cnt_between_sel(const uint aLbCnt, const uint aUbCnt, const double aSelectivity) const {
  return 1;
}

double
EstHavingFentSn::estimate_sum_eq(const double aSumB) const {
  double l = 0, u = 0;
  l = _dist_sum.get_cdf(aSumB - 0.5);
  u = _dist_sum.get_cdf(aSumB + 0.5);
  return std::max<double>(1, ((double) nodv_A() * (u - l)));
}

double
EstHavingFentSn::estimate_sum_eq_sel(const double aSumB, const double aSelectivity) const {
  return 1;
}

double
EstHavingFentSn::estimate_sum_eq_cnt(const double aSumB, const uint aCntLb, const uint aCntUb) const {
  return 1;
}

double
EstHavingFentSn::estimate_sum_between(const double aLbSumB, const double aUbSumB) const {
  const double l = _dist_sum.get_cdf(aLbSumB - 0.5);
  const double u = _dist_sum.get_cdf(aUbSumB + 0.5);
  return ((double) nodv_A() * (u - l));
}

double
EstHavingFentSn::estimate_sum_between_sel(const double aLbSumB, const double aUbSumB, const double aSelectivity) const {
  return 1;
}

double
EstHavingFentSn::estimate_sum_between_and_cnt_between(const double aLbSumB, const double aUbSumB,
                                                      const uint   aLbCnt,  const uint   aUbCnt) const {
  return 1;
}



