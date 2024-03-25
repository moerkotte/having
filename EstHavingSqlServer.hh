#ifndef EST_HAVING_SQL_SERVER_HH
#define EST_HAVING_SQL_SERVER_HH

#include "util.hh"
#include "gm_dist.hh"


/*
 *  class EstHavingSqlServer
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

class EstHavingSqlServer {
  public:
    using dist_param_t = gm::dist_param_t;
    using attr_stat_t  = gm::attr_stat_t;
    using dist_t       = gm::dist_normal_t;
  public:
    EstHavingSqlServer();
    EstHavingSqlServer(const uint64_t  aCardR,
                    const attr_stat_t& aAttrStatA,
                    const attr_stat_t& aAttrStatB,
                    const attr_stat_t& aAttrStatCnt);
  public:
    double estimate_cnt_eq(const uint aCnt) const;
    double estimate_cnt_between(const uint aLbCnt, const uint aUbCnt) const;
  public:
    double selectivity_cnt_eq(const uint aCnt) const;
    double selectivity_cnt_between(const uint aLbCnt, const uint aUbCnt) const;
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
  public:
    std::ostream& print(std::ostream& os) const;
  private:
    uint64_t     _card_R;         // cardinality of grouped relation
    attr_stat_t  _attr_stat_A;    // statistics of attribute in group by A
    attr_stat_t  _attr_stat_B;    // statistics of attribute in having sum(B) theta k | between k1 and k2
    dist_param_t _dist_param_Cnt; // guessed parameters for distribution of group sizes
    dist_t       _dist_Cnt;       // distribution of group sizes i.e. count(*) values
    dist_t       _dist_B;         // distribution to model attribute B
};


#endif

