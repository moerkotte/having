#ifndef EST_HAVING_FENT_HH
#define EST_HAVING_FENT_HH

#include "util.hh"
#include "gm_dist.hh"


/*
 *  class EstHavingFentSn
 *  estimates the cardinality of
 *  select ...
 *  from   R
 *  group by A
 *  having sum(B) in [.,.]a
 *  using the method of Fent using a skew-normal distribution
 */

/*
  Query Pattern:
  select ...
  from   R
  [where p]
  group by A
  having  aggr(B) ...
*/

class EstHavingFentSn {
  public:
    using dist_param_t = gm::dist_param_t;
    using attr_stat_t  = gm::attr_stat_t;
    using dist_t       = gm::dist_skew_normal_t;
  public:
    EstHavingFentSn();
    EstHavingFentSn(const uint64_t     aCardR,
                    const attr_stat_t& aAttrStatA,
                    const attr_stat_t& aAttrStatB,
                    const attr_stat_t& aAttrStatCnt); // unused
  public:
    double estimate_cnt_eq(const uint aCnt) const;
    double estimate_cnt_between(const uint aLbCnt, const uint aUbCnt) const;
    double estimate_cnt_eq_sel(const uint aCnt, const double aSelectivity) const;
    double estimate_cnt_between_sel(const uint aLbCnt, const uint aUbCnt, const double aSelectivity) const;
  public:
    double estimate_sum_eq(const double aSumB) const;
    double estimate_sum_eq_sel(const double aSumB, const double aSelectivity) const;
    double estimate_sum_eq_cnt(const double aSumB, const uint aCntLb, const uint aCntUb) const;
    double estimate_sum_between(const double aLbSumB, const double aUbSumB) const;
    double estimate_sum_between_sel(const double aLbSumB, const double aUbSumB, const double aSelectivity) const;
    double estimate_sum_between_and_cnt_between(const double aLbSumB, const double aUbSumB, 
                                                const uint   aLbCnt,  const uint   aUbCnt) const;
  public:
    double estimate_avg_between(const double aLbAvgB, const double aUbAvgB) const;
    double estimate_avg_between_sel(const double aLbAvgB, const double aUbAvgB, const double aSelectivity) const;
    double estimate_avg_between_cnt(const double aLbAvgB, const double aUbAvgB, 
                                    const uint   aLbCnt,  const uint   aUbCnt) const;
  public:
    double estimate_min_eq(const double aMinB) const;
    double estimate_min_eq_sel(const double aMinB, const double aSelectivity) const;
    double estimate_min_eq_cnt(const double aMinB, const uint aCntLb, const uint aCntUb) const;
    double estimate_min_between(const double aLbMinB, const double aUbMinB) const;
    double estimate_min_between_sel(const double aLbMinB, const double aUbMinB, const double aSelectivity) const;
    double estimate_min_between_cnt(const double aLbMinB, const double aUbMinB,
                                    const uint   aLbCnt,  const uint   aUbCnt) const;
  public:
    double estimate_max_eq(const double aMaxB) const;
    double estimate_max_eq_sel(const double aMaxB, const double aSelectivity) const;
    double estimate_max_eq_cnt(const double aMaxB, const uint aCntLb, const uint aCntUb) const;
    double estimate_max_between(const double aLbMaxB, const double aUbMaxB) const;
    double estimate_max_between_sel(const double aLbMaxB, const double aUbMaxB, const double aSelectivity) const;
    double estimate_max_between_cnt(const double aLbMaxB, const double aUbMaxB,
                                    const uint   aLbCnt,  const uint   aUbCnt) const;
  public:
    inline uint64_t card_R()      const { return _card_R; }
    // grouping attribute A
    inline const attr_stat_t& attr_stat_A() const { return _attr_stat_A; }
    inline uint64_t nodv_A()      const { return _attr_stat_A.nodv(); }
    inline double   skew_A()      const { return _attr_stat_A.skew(); }
    // group sizes, i.e. count(*) values
    inline double mean_Cnt() const { return _dist_param_Cnt.mean(); }
    inline double var_Cnt()  const { return _dist_param_Cnt.variance(); }
    inline double skew_Cnt() const { return _dist_param_Cnt.skew(); }
    // aggregated attribute B
    inline const attr_stat_t&  attr_stat_B()  const { return _attr_stat_B; }
    inline const dist_param_t& dist_param_B() const { return _attr_stat_B.dist_param(); }
    inline uint64_t nodv_B()      const { return _attr_stat_B._nodv; }
    inline double   min_B()       const { return _attr_stat_B._min; }
    inline double   max_B()       const { return _attr_stat_B._max; }
    inline double   avg_dist_B()  const { return ((max_B() - min_B()) / (nodv_B() - 1)); }
    inline double   mean_B()      const { return _attr_stat_B.mean(); }
    inline double   std_dev_B()   const { return _attr_stat_B.std_dev(); }
    inline double   variance_B()  const { return _attr_stat_B.variance(); }
    inline bool     is_ui_dense_B() const { return _attr_stat_B.is_ui_dense(); }
  private:
    dist_t get_dist_Cnt();
    dist_t get_dist_B();
    dist_t get_dist_sum();
  public:
    std::ostream& print(std::ostream& os) const;
  private:
    uint64_t     _card_R;         // cardinality of grouped relation
    attr_stat_t  _attr_stat_A;    // statistics of attribute in group by A
    attr_stat_t  _attr_stat_B;    // statistics of attribute in having sum(B) theta k | between k1 and k2
    dist_param_t _dist_param_Cnt; // guessed parameters for distribution of group sizes
    dist_t       _dist_Cnt;       // distribution of group sizes i.e. count(*) values
    dist_t       _dist_B;         // distribution to model attribute B
    dist_param_t _dist_param_sum; // distribution parameters to model sum 
    dist_t       _dist_sum;       // _dist_G * _dist_B to model distribution of sum(B) values
};


#endif

