#ifndef SIMPLE_PROFILE_GROUPBY_HAVING_HH
#define SIMPLE_PROFILE_GROUPBY_HAVING_HH

#include "util.hh"
// #include "norm_dist.hh"

#include "gm_dist.hh"
#include "int_comp.hh"

/*
 *  class SimpleProfileGbH
 *  estimates the cardinality of
 *  select ...
 *  from   R
 *  group by A
 *  having sum(B) in [.,.]
 *  additionally, the having clause can be extended by
 *    and count(*) in [.,.]
 *  and a where clause with a predicate can be added
 *
 *  we use the central limit theorem and only very few numbers
 *  (details see papers/CardEstHaving)
 * the parameter aUiDense indicates that the dom(B) is [1,m] for some m
 * and that all numbers in this range occur
 */


class SimpleProfileGbH {
  public:
    using attr_stat_t = gm::attr_stat_t;
  public:
    SimpleProfileGbH();
    SimpleProfileGbH(const uint64_t     aCardR,
                     const attr_stat_t& aAttrStatA,
                     const attr_stat_t& aAttrStatB,
                     const attr_stat_t& aAttrStatCnt);
  public:
    double estimate_cnt_eq(const uint aCnt) const;
    double estimate_cnt_between(const uint aLbCnt, const uint aUbCnt) const;
    double estimate_cnt_eq_sel(const uint aCnt, const double aSelectivity) const;
    double estimate_cnt_between_sel(const uint aLbCnt, const uint aUbCnt, const double aSelectivity) const;
  public:
    double estimate_sum_eq(const double aSumB) const;
    double estimate_sum_eq_sel(const double aSumB, 
                               const double aSelectivity) const;
    double estimate_sum_eq_cnt(const double aSumB,
                               const uint aCntLb, const uint aCntUb) const;
    double estimate_sum_between(const double aLbSumB, const double aUbSumB) const;
    double estimate_sum_between_sel(const double aLbSumB, const double aUbSumB,
                                    const double aSelectivity) const;
    double estimate_sum_between_and_cnt_between(const double aLbSumB, const double aUbSumB, 
                                                const uint   aLbCnt,  const uint   aUbCnt) const;
    double estimate_sum_between_or_cnt_between(const double aLbSumB, const double aUbSumB, 
                                               const uint   aLbCnt,  const uint   aUbCnt) const;
  public:
    // special implementations using integer compositions (ic)
    // unlimited use of integer compositions (Sec. 4.4: Integer B: Counting Integer Compositions)
    double estimate_sum_eq_ic(const uint aSumB) const;
    // aLim: limit for special formulas (1 <= aLim <= 3)
    // for k > aLim: use general formula with normal distribution (Sec. 4.5 Eq. 30)
    double estimate_sum_eq_ic_lim(const uint aSumB, const uint aLim) const;
  public:
    double estimate_avg_between(const double aLbAvgB, const double aUbAvgB) const;
    double estimate_avg_between_sel(const double aLbAvgB, const double aUbAvgB,
                                    const double aSelectivity) const;
    double estimate_avg_between_and_cnt_between(const double aLbAvgB, const double aUbAvgB, 
                                                const uint   aLbCnt,  const uint   aUbCnt) const;
    double estimate_avg_between_or_cnt_between(const double aLbAvgB, const double aUbAvgB, 
                                               const uint   aLbCnt,  const uint   aUbCnt) const;
  public:
    double estimate_min_eq(const double aMinB) const;
    double estimate_min_eq_sel(const double aMinB,
                               const double aSelectivity) const;
    double estimate_min_eq_cnt(const double aMinB,
                               const uint aCntLb, const uint aCntUb) const;
    double estimate_min_between(const double aLbMinB, const double aUbMinB) const;
    double estimate_min_between_sel(const double aLbMinB, const double aUbMinB,
                                    const double aSelectivity) const;
    double estimate_min_between_and_cnt_between(const double aLbMinB, const double aUbMinB,
                                                const uint   aLbCnt,  const uint   aUbCnt) const;
  public:
    double estimate_max_eq(const double aMaxB) const;
    double estimate_max_eq_sel(const double aMaxB, const double aSelectivity) const;
    double estimate_max_eq_cnt(const double aMaxB, const uint aCntLb, const uint aCntUb) const;
    double estimate_max_between(const double aLbMaxB, const double aUbMaxB) const;
    double estimate_max_between_sel(const double aLbMaxB, const double aUbMaxB, const double aSelectivity) const;
    double estimate_max_between_and_cnt_between(const double aLbMaxB, const double aUbMaxB,
                                                const uint   aLbCnt,  const uint   aUbCnt) const;
  public:
    double est_sum_eq(const uint aK, const double aVal) const;
    double est_sum_1_eq(const double aVal) const;
    double est_sum_2_eq(const double aVal) const;
    double est_sum_3_eq(const double aVal) const;
    double est_sum_k_eq(const uint aK, const double aVal) const;
    double est_sum_between(const uint aK, const double aLb, const double aUb) const;
    double est_sum_1_between(const double aLb, const double aUb) const;
    double est_sum_2_between(const double aLb, const double aUb) const;
    double est_sum_k_between(const uint aK, const double aLb, const double aUb) const;
  public:
    inline uint64_t card_R()      const { return _card_R; }
    // count(*)
    inline const attr_stat_t& attr_stat_Cnt() const { return _attr_stat_Cnt; }
    inline uint64_t min_cnt()     const { return attr_stat_Cnt().min(); }
    inline uint64_t max_cnt()     const { return attr_stat_Cnt().max(); }
    inline uint64_t nodv_cnt()    const { return (max_cnt() - min_cnt() + 1); }
    // grouping attribute A
    inline const attr_stat_t& attr_stat_A() const { return _attr_stat_A; }
    inline uint64_t nodv_A()      const { return _attr_stat_A._nodv; }
    // aggregated attribute B
    inline const attr_stat_t& attr_stat_B() const { return _attr_stat_B; }
    inline uint64_t nodv_B()      const { return _attr_stat_B._nodv; }
    inline double   min_B()       const { return _attr_stat_B._min; }
    inline double   max_B()       const { return _attr_stat_B._max; }
    inline double   avg_dist_B()  const { return ((max_B() - min_B()) / (nodv_B() - 1)); }
    inline double   mean_B()      const { return _attr_stat_B.mean(); }
    inline double   std_dev_B()   const { return _attr_stat_B.std_dev(); }
    inline double   variance_B()  const { return _attr_stat_B.variance(); }
    inline bool     is_ui_dense_B() const { return _attr_stat_B.is_ui_dense(); }
  public:
    // for one particular counter value k
    inline double freq_tot_k_uda() const { return (((double) nodv_A()) / ((double) nodv_cnt())); }
    inline double freq_tot_k(const uint) const { return (((double) nodv_A()) / ((double) nodv_cnt())); }
    inline double min_B_k(const uint k) const { return (((double) k) * min_B()); }
    inline double max_B_k(const uint k) const { return (((double) k) * max_B()); }
    inline double mean_B_k(const uint k) const { return (((double) k) * mean_B()); }
    inline double std_dev_B_k(const uint k) const { return std::sqrt(variance_B_k(k)); }
    inline double variance_B_k(const uint k) const { return (((double) k) * variance_B()); }
  public:
    // probability of one given B-value assuming uniform distr. assumption
    inline double get_prob_B() const { return (double{1} / ((double) nodv_B())); }
    inline double get_p() const { return get_prob_B(); }
    inline double pdf_k(const uint k, const double a) const;
    inline double cdf_k(const uint k, const double a) const;
    // intersect range [a,b] for B with [k min_B, k max_B], a, b: in&out params
    void   cut_B_k(const uint k, double& a, double& b) const;
    // intersect range [a,b] for counts with [min_cnt(), max_cnt()]
    void   cut_cnt(uint& aLb, uint& aUb) const;
    // survival probability of a group with k elements after a selection with aSelectivity 
    double sp(const uint k, const double aSelectivity) const;
  public:
    std::ostream& print(std::ostream& os) const;
  private:
    uint64_t    _card_R;        // cardinality of grouped relation
    attr_stat_t _attr_stat_A;   // statistics of attribute in group by A
    attr_stat_t _attr_stat_B;   // statistics of attribute in having sum(B) theta k | between k1 and k2
    attr_stat_t _attr_stat_Cnt; // statistics of count(*)
};


#ifdef USE_BOOST
double
SimpleProfileGbH::pdf_k(const uint k, const double a) const {
  dist_t lDist(mean_B_k(k), std_dev_B_k(k));
   return pdf(lDist, a);
}

double
SimpleProfileGbH::cdf_k(const uint k, const double a) const {
  dist_t lDist(mean_B_k(k), std_dev_B_k(k));
  return cdf(lDist, a);
}
#else
double
SimpleProfileGbH::pdf_k(const uint k, const double a) const {
   return normal_pdf(a, mean_B_k(k), std_dev_B_k(k));
}

double
SimpleProfileGbH::cdf_k(const uint k, const double a) const {
   return normal_cdf(a, mean_B_k(k), std_dev_B_k(k));
}
#endif


#endif

