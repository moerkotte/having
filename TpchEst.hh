#ifndef TPCH_EST_HH
#define TPCH_EST_HH

#include "util.hh"
// #include "infra/tmath.hh"
#include "variance.hh"
#include <fstream>

/*
 * estimate result cardinality of
 * select ..
 * from R
 * group by R.A
 * having sum(B) in [.,.] 
 * and some extensions
 *
 * see papers/CardEstHaving
 *
 * contains also the materialization of 
create view no_sum_count as
select no_line, sum_quant, count(*) as cnt
from (select l_orderkey, count(*) as no_line, sum(l_quantity) as sum_quant
      from   lineitem
      group by l_orderkey)
group by no_line, sum_quant
order by no_line, sum_quant;
select * from no_sum_count;
  which is contained in nosumcount.dat
  to produce precise result cardinalities
 */


/*
 * Estimator using result of Q_c for estimation
 * specialized for TPC-H1 Q18' (SF=1) (see papers/CardEstHaving)
 * hardcoded for l_quantity in [1,50]
 */

typedef std::vector<double_vt> double_vvt;

class TpchEst {
  public:
    #ifdef USE_BOOST
    typedef boost::math::normal_distribution<double> dist_t;
    #endif
  public:
    TpchEst();
    ~TpchEst();
  public:
    double estimate_sum_eq(const uint aSumQuantity) const;
    double estimate_sum_between(const uint aSumQuantityLb, const uint aSumQuantityUb) const;
    double estimate_sum_eq_sel(const uint aSumQuantity, const double aSelectivity) const;
    double estimate_sum_between_sel(const uint aLb, const uint aUb, const double aSelectivity) const;
    double estimate_sum_eq_cnt(const uint aSumQuantity, const uint aCntLb, const uint aCntUb) const;
    double estimate_sum_between_and_cnt_between(const uint aLb, const uint aUb, const uint aCntLb, const uint aCntUb) const;
    double estimate_sum_between_or_cnt_between(const uint aLb, const uint aUb, const uint aCntLb, const uint aCntUb) const;
  public:
    double est_1_eq(const uint aVal) const;
    double est_2_eq(const uint aVal) const;
    double est_k_eq(const uint aK, const uint aVal) const;
    double est_1_between(const uint aLb, const uint aUb) const;
    double est_2_between(const uint aLb, const uint aUb) const;
    double est_k_between(const uint aK, const uint aLb, const uint aUb) const;
    void   cut_k(const uint k, uint& a, uint& b) const;
    uint   min_B_k(const uint aK) const;
    uint   max_B_k(const uint aK) const;
    uint   get_true_card_cnt_eq(const uint aVal) const;
    uint   get_true_card_cnt_between(const uint aLb, const uint aUb) const;
    uint   get_true_card_sum_eq(const uint aVal) const;
    uint   get_true_card_sum_eq_and_cnt_between(const uint aVal, const uint aCntLb, const uint aCntUb) const;
    uint   get_true_card_sum_eq_or_cnt_between(const uint aVal, const uint aCntLb, const uint aCntUb) const;
    uint   get_true_card_sum_between(const uint aLb, const uint aUb) const;
    uint   get_true_card_sum_between_and_cnt_between(const uint aLb, const uint aUb, const uint aCntLb, const uint aCntUb) const;
    uint   get_true_card_sum_between_or_cnt_between(const uint aLb, const uint aUb, const uint aCntLb, const uint aCntUb) const;
    uint   get_true_card_avg_between(const uint aLb, const uint aUb) const;
    uint   get_true_card_avg_between_cnt(const uint aLb, const uint aUb, const uint aCntLb, const uint aCntUb) const;
    uint   get_true_card_avg_between_one_count(const uint aLb, const uint aUb, const uint aCnt) const;
  public:
    // accessors
    inline uint64_t card_R()   const { return _card_R; }
    inline double   card_R_d() const { return ((double) _card_R); }
    inline double   freq_tot(const uint k) const { return _freq_tot[k]; }
    inline double   freq_sum(const uint k, const uint c) const { return _freq_sum[k][c]; }
    inline double   mean(const uint k) const { return _mean[k]; }
    inline double   std_dev(const uint k) const { return _std_dev[k]; }
  public:
    // probability of one value of l_quantity in [1,50] to occur (uniform distr. assumption)
    inline double get_p() const { return ((double) 1.0 / 50.0); }
    inline double pdf_k(const uint k, const double a) const;
    inline double cdf_k(const uint k, const double a) const;
  public:
    std::ostream& print_a(std::ostream& os) const;
  private:
    uint64_t   _card_R;   // cardinality of Lineitem
    double_vt  _freq_tot; // frequencies of group size: _freq_tot[k] = #groups of size k
    uint_vvt   _freq_sum; // frequencies for sum(l_quantity) = x, for all count-values, for all sum-values
    double_vvt _prob_sum; // _prob_sum[k][c] = _freq_sum[k][c] / _freq_tot[k]
    double_vt  _mean;     // mean   of _freq_sum[k][.]
    double_vt  _std_dev;  // stddev of _freq_sum[k][.]
};

#endif

