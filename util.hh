#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdint>
#include <assert.h>
#include <cmath>

using uint = unsigned int;
using rc_t = int;
using uint_vt = std::vector<uint>;
using uint_vvt = std::vector<uint_vt>;

using double_vt = std::vector<double>;

double fBinomGamma(double n, double k);
double fYaoProbGamma(double N, double n, double k);
double normal_pdf(const double x, const double aMean, const double aStdDev);
double normal_cdf(const double x, const double aMean, const double aStdDev);

inline double
uni_dist_mean(const double aMin, const double aMax) {
  return ((aMin + aMax) / 2);
}

inline double
uni_dist_variance(const double aMin, const double aMax) {
  return ((((double) 1) / ((double) 12)) * ((aMax - aMin) * (aMax - aMin)));
}


enum aggrfun_et {
  kAggrCnt = 0,
  kAggrSum = 1,
  kAggrAvg = 2,
  kAggrMin = 3,
  kAggrMax = 4,
  kNoAggr  = 5
};

enum predkind_et {
  kPredEq = 0,
  kPredRg = 1,
  kNoPred = 2
};

enum andor_et {
  kAnd     = 0,
  kOr      = 1,
  kNoAndOr = 2
};

struct querytemplate_t {
  querytemplate_t();
  uint32_t _aggr_having       :  3;
  uint32_t _pred_having_aggr  :  2;
  uint32_t _pred_having_cnt   :  2; // additional count predicate in having clause
  uint32_t _bool_connector    :  2; // having aggr(B) [and/or count(*)] ...
  uint32_t _pred_where        :  2; // [where l_suppkey mod :e < :f]
  uint32_t _unused            : 21;

  uint32_t aggr_having()      const { return _aggr_having; }
  uint32_t pred_having_aggr() const { return _pred_having_aggr; }
  uint32_t pred_having_cnt()  const { return _pred_having_cnt; }
  uint32_t bool_connector()   const { return _bool_connector; }
  uint32_t pred_where()       const { return _pred_where; }

  // set to specific templates
  // template  1-17: true cardinality cannot be determined by TpchEst
  // template 18-25: true cardinality can    be determined by TpchEst
  static inline uint no_templates_to_gen() { return 17; }
  static inline uint no_templates() { return 25; }
  void set(int i);
  // true card from sqlite/DuckDb
  void set_1();  // count_eq_sel
  void set_2();  // count_rg_sel
  void set_3();  // sum_eq_sel
  void set_4();  // sum_rg_sel
  void set_5();  // avg_rg_sel
  void set_6();  // min_eq_sel
  void set_7();  // min_rg_sel
  void set_8();  // max_eq_sel
  void set_9();  // max_rg_sel
  void set_10(); // min_eq
  void set_11(); // min_eq_and_cnt_rg
  void set_12(); // min_rg
  void set_13(); // min_rg_and_cnt_rg
  void set_14(); // max_eq
  void set_15(); // max_eq_and_cnt_rg
  void set_16(); // max_rg
  void set_17(); // max_rg_and_cnt_rg
  // true card from TpchEst
  void set_18(); // cnt_eq
  void set_19(); // cnt_rg
  void set_20(); // sum_eq
  void set_21(); // sum_eq_and_cnt_rg
  void set_22(); // sum_rg
  void set_23(); // sum_rg_and_cnt_rg
  void set_24(); // avg_rg
  void set_25(); // avg_rg_and_cnt_rg

  // to generate a string for, e.g., filenames
  std::string to_string() const;
};
struct queryinstance_t {
  uint32_t        _params[7];
  queryinstance_t();
  uint32_t a() const { return _params[0]; } // having aggr(.) = a()
  uint32_t b() const { return _params[1]; } // having aggr(.) between a() and b()
  uint32_t c() const { return _params[2]; } // having ... and/or count(*) = c()
  uint32_t d() const { return _params[3]; } // having ... and/or count(*) between c() and d()
  uint32_t e() const { return _params[4]; } // where l_suppkey % e() < f()
  uint32_t f() const { return _params[5]; } // where l_suppkey % e() < f()
  uint32_t& get_a() { return _params[0]; } // having aggr(.) = a()
  uint32_t& get_b() { return _params[1]; } // having aggr(.) between a() and b()
  uint32_t& get_c() { return _params[2]; } // having ... and/or count(*) = c()
  uint32_t& get_d() { return _params[3]; } // having ... and/or count(*) between c() and d()
  uint32_t& get_e() { return _params[4]; } // where l_suppkey % e() < f()
  uint32_t& get_f() { return _params[5]; } // where l_suppkey % e() < f()
  queryinstance_t&  a(const uint x) { _params[0] = x; return (*this); }
  queryinstance_t&  b(const uint x) { _params[1] = x; return (*this); }
  queryinstance_t&  c(const uint x) { _params[2] = x; return (*this); }
  queryinstance_t&  d(const uint x) { _params[3] = x; return (*this); }
  queryinstance_t&  e(const uint x) { _params[4] = x; return (*this); }
  queryinstance_t&  f(const uint x) { _params[5] = x; return (*this); }

  // get selectivity of selection predicate [l_suppkey mod e() < f()]
  inline double sel() const { return (((double) f()) / ((double) e())); }

  std::ostream& print_param(std::ostream& os, const querytemplate_t& aQueryTemplate, const int aWidth) const;
  std::ostream& print_query(std::ostream& os, const querytemplate_t& aQueryTemplate) const;
  std::ostream& print_select(std::ostream& os, const querytemplate_t& aQueryTemplate) const;
  std::ostream& print_where(std::ostream& os, const querytemplate_t& aQueryTemplate) const;
  std::ostream& print_having(std::ostream& os, const querytemplate_t& aQueryTemplate) const;

};
typedef std::vector<queryinstance_t> queryinstance_vt;

static_assert(4 == sizeof(querytemplate_t), "sizeof(querytemplate_t)!" );
std::string to_string_aggr(const uint a);
std::string to_string_pred(const uint a);
std::string to_string_bool_connector(const uint a);
std::string to_string_pred_where(const uint a);

namespace mt {

template<class Tfloat>
inline Tfloat
roundt(const Tfloat);


template<>
inline double
roundt<double>(const double x) {
  return round(x);
}

template<>
inline float
roundt<float>(const float x) {
  return roundf(x);
}

template<class Tfloat>
inline Tfloat
roundXt(const Tfloat x) {
  return (roundt<Tfloat>(x * (Tfloat) 10) / (Tfloat) 10);
}   

template<class Tfloat>
inline Tfloat
roundXXt(const Tfloat x) {
  return (roundt<Tfloat>(x * (Tfloat) 100) / (Tfloat) 100);
}   

template<class Tfloat>
inline Tfloat
roundXXXt(const Tfloat x) {
  return (roundt<Tfloat>(x * (Tfloat) 1000) / (Tfloat) 1000);
}

}


namespace q {

template<class Tfloat>
Tfloat
qerror(const Tfloat x, const Tfloat y) {
  if(x == y) { return (Tfloat) 1; }
  if(x >= y) { 
    return (x/y);
  }
  return (y/x);
} 

template<class Tfloat>
inline Tfloat
q_error_safe(const Tfloat x, const Tfloat y) {
  const Tfloat a = (x == 0) ? 1 : x;
  const Tfloat b = (y == 0) ? 1 : y;
  return ((a < b) ? (b/a) : (a/b));
}   
  
}
