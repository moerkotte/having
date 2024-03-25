#include "EstHavingSqlServer.hh"


EstHavingSqlServer::EstHavingSqlServer()
                : _card_R(0), 
                  _attr_stat_A(), 
                  _attr_stat_B(), 
                  _dist_param_Cnt(),
                  _dist_Cnt(), 
                  _dist_B() {}

EstHavingSqlServer::EstHavingSqlServer(const uint64_t     aCardR,
                                 const attr_stat_t& aAttrStatA,
                                 const attr_stat_t& aAttrStatB,
                                 const attr_stat_t& aAttrStatCnt)
                : _card_R(aCardR),
                  _attr_stat_A(aAttrStatA), 
                  _attr_stat_B(aAttrStatB), 
                  _dist_param_Cnt(),
                  _dist_Cnt(), 
                  _dist_B() {
  _dist_Cnt = get_dist_Cnt();
  _dist_B   = get_dist_B();
}

EstHavingSqlServer::dist_t
EstHavingSqlServer::get_dist_Cnt() {
  // 1. set _dist_param_Cnt
  const double C = ((double) card_R());
  const double D = ((double) nodv_A());
  const double M = C/D;
  const double Mp = M * ((D-1)/D); // assumed mean
  const double S = std::sqrt(M);   // assumed standard deviation
  return dist_t(0, D, Mp, S);
}

EstHavingSqlServer::dist_t
EstHavingSqlServer::get_dist_B() {
  return dist_t(min_B(), max_B(), mean_B(), std_dev_B());
}

double
EstHavingSqlServer::estimate_cnt_eq(const uint aCnt) const {
  return ((double) nodv_A() * selectivity_cnt_eq(aCnt));
}

double
EstHavingSqlServer::estimate_cnt_between(const uint aLbCnt, const uint aUbCnt) const {
  return ((double) nodv_A() * selectivity_cnt_between(aLbCnt, aUbCnt));
}

double
EstHavingSqlServer::selectivity_cnt_eq(const uint aCnt) const {
  const double x = (double) aCnt;
  double l = 0, u = 0;
  if(1 == aCnt) {
    l = _dist_Cnt.get_cdf(1.0);
    u = _dist_Cnt.get_cdf(1.5);
  } else {
    l = _dist_Cnt.get_cdf(x - 0.5);
    u = _dist_Cnt.get_cdf(x + 0.5);
  }
  return (u - l);
}

double
EstHavingSqlServer::selectivity_cnt_between(const uint aLbCnt, const uint aUbCnt) const {
  double l = 0, u = 0;
  if(1 >= aLbCnt) {
    l = _dist_Cnt.get_cdf(1.0);
    u = _dist_Cnt.get_cdf((double) aUbCnt + 0.5);
  } else {
    l = _dist_Cnt.get_cdf((double) aLbCnt - 0.5);
    u = _dist_Cnt.get_cdf((double) aUbCnt + 0.5);
  }
  return (u - l);
}


