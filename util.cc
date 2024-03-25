#include "util.hh"

double
fLogFactorial(double x) {
  return std::lgamma(x+1);
}

double
fBinomGamma(double n, double k) {
  if(k <= 0.000001) { return 1.0; }
  return floor(0.5 + exp((fLogFactorial(n)-fLogFactorial(k))-fLogFactorial(n-k)));
}
double
fLogBinomGamma(double n, double k) {
  if(k <= 0.000001) { return 0.0; }
  return ((fLogFactorial(n)-fLogFactorial(k))-fLogFactorial(n-k));
}

double
fYaoProbGamma(double N, double n, double k) {
  if((N - n) < k) {
    return 1.0;
  }
  return (1.0 - exp(fLogBinomGamma(N-n,k) - fLogBinomGamma(N,k)));
}

double
dihrLower(double k, double n, double m) {
  const double B = n/m; 
  if(k > (n - B)) return m;
  double f = (1.0 - (k / (n - ((B - 1.0) / 2.0))));
  #ifdef DIHR_LOWER_POW
  const double p = pow(f, B);
  #else
  const double p = exp(log(f) * B); // pow(f, B);
  #endif
  const double lRes = (m * (1 - p));
  return (lRes > m) ? m : lRes;
}

namespace gm {
double
calc_yao(const double k, const double n, const double m) {
  double lRes = dihrLower(k, n, m);
  if((1.0 > lRes) && (0 < k)) {
    lRes = 1;
  }
  return lRes;
}
}

  
double
  normal_pdf(const double x, const double aMean, const double aStdDev) {
    const double lVariance = aStdDev * aStdDev;
    return (  (1.0 / (aStdDev * sqrt(2.0 * M_PI)))
            * (pow(M_E, -1.0 * ((pow((x - aMean), 2.0)) / ( lVariance * 2.0)))));
  }
double
  normal_cdf(const double x, const double aMean, const double aStdDev) {
    // const double z = (x - aMean) / (aStdDev * sqrt(2.0));
    // const double a = (0.5 + 0.5 * erf(z));
    const double y = (x - aMean) / aStdDev;
    const double b = (0.5 * erfc(-y * M_SQRT1_2));
    // printf("ab = %12f %12f\n", a, b);
    // assert(is_equal(a, b, 1e-7));
    return b;
  }


querytemplate_t::querytemplate_t() 
                : _aggr_having(kNoAggr),
                  _pred_having_aggr(kNoPred),
                  _pred_having_cnt(kNoPred),
                  _bool_connector(kNoAndOr),
                  _pred_where(kNoPred),
                  _unused(0) {
}

void
querytemplate_t::set(int i) {
  switch(i) {
    case  1: set_1();  break;
    case  2: set_2();  break;
    case  3: set_3();  break;
    case  4: set_4();  break;
    case  5: set_5();  break;
    case  6: set_6();  break;
    case  7: set_7();  break;
    case  8: set_8();  break;
    case  9: set_9();  break;
    case 10: set_10(); break;
    case 11: set_11(); break;
    case 12: set_12(); break;
    case 13: set_13(); break;
    case 14: set_14(); break;
    case 15: set_15(); break;
    case 16: set_16(); break;
    case 17: set_17(); break;
    case 18: set_18(); break;
    case 19: set_19(); break;
    case 20: set_20(); break;
    case 21: set_21(); break;
    case 22: set_22(); break;
    case 23: set_23(); break;
    case 24: set_24(); break;
    case 25: set_25(); break;
    default: assert(0 == 1);
  }
}

// template 1:
// where l_suppkey % :e < :f
// group by l_orderkey
// having count(l_quantity) = :a

void
querytemplate_t::set_1() {
  _aggr_having       = kAggrCnt; // having count(l_quantity) =
  _pred_having_aggr  = kPredEq;  // having count(l_quantity) =
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kPredRg;  // where l_suppkey % e < f
}

// template 2:
// where l_suppkey % :e < :f
// group by l_orderkey
// having count(l_quantity) between :a and :b

void
querytemplate_t::set_2() {
  _aggr_having       = kAggrCnt; // having count(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having count(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kPredRg;  // where l_suppkey % e < f
}

// template 3:
// where l_suppkey % :e < :f
// group by l_orderkey
// having sum(l_quantity) = :a

void
querytemplate_t::set_3() {
  _aggr_having       = kAggrSum; // having sum(l_quantity) =
  _pred_having_aggr  = kPredEq;  // having sum(l_quantity) =
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kPredRg;  // where l_suppkey % e < f
}

// template 4:
// where l_suppkey % :e < :f
// group by l_orderkey
// having sum(l_quantity) between :a and :b

void
querytemplate_t::set_4() {
  _aggr_having       = kAggrSum; // having sum(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having sum(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kPredRg;  // where l_suppkey % e < f
}

// template 5:
// where l_suppkey % :e < :f
// group by l_orderkey
// having avg(l_quantity) between :a and :b

void
querytemplate_t::set_5() {
  _aggr_having       = kAggrAvg; // having avg(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having avg(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kPredRg;  // where l_suppkey % e < f
}

// template 6:
// where l_suppkey % :e < :f
// group by l_orderkey
// having min(l_quantity) = :a

void
querytemplate_t::set_6() {
  _aggr_having       = kAggrMin; // having min(l_quantity) =
  _pred_having_aggr  = kPredEq;  // having min(l_quantity) =
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kPredRg;  // where l_suppkey % e < f
}

// template 7:
// where l_suppkey % :e < :f
// group by l_orderkey
// having min(l_quantity) between :a and :b

void
querytemplate_t::set_7() {
  _aggr_having       = kAggrMin; // having min(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having min(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kPredRg;  // where l_suppkey % e < f
}

// template 8:
// where l_suppkey % :e < :f
// group by l_orderkey
// having max(l_quantity) = :a

void
querytemplate_t::set_8() {
  _aggr_having       = kAggrMax; // having max(l_quantity) =
  _pred_having_aggr  = kPredEq;  // having max(l_quantity) =
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kPredRg;  // where l_suppkey % e < f
}

// template 9:
// where l_suppkey % :e < :f
// group by l_orderkey
// having max(l_quantity) between :a and :b

void
querytemplate_t::set_9() {
  _aggr_having       = kAggrMax; // having max(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having max(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kPredRg;  // where l_suppkey % e < f
}

// template 10:
// group by l_orderkey
// having min(l_quantity) = :a

void
querytemplate_t::set_10() {
  _aggr_having       = kAggrMin; // having min(l_quantity) between
  _pred_having_aggr  = kPredEq;  // having min(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kNoPred;  // empty where clause
}

// template 11:
// group by l_orderkey
// having min(l_quantity) = :a and count(*) between :c and :d

void
querytemplate_t::set_11() {
  _aggr_having       = kAggrMin; // having min(l_quantity) between
  _pred_having_aggr  = kPredEq;  // having min(l_quantity) between
  _pred_having_cnt   = kPredRg;  // and count(*) between :c and :d
  _bool_connector    = kAnd;     // and count(*) between :c and :d
  _pred_where        = kNoPred;  // empty where clause
}


// template 12:
// group by l_orderkey
// having min(l_quantity) between :a and :b

void
querytemplate_t::set_12() {
  _aggr_having       = kAggrMin; // having min(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having min(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kNoPred;  // empty where clause
}

// template 13:
// group by l_orderkey
// having min(l_quantity) between :a and :b and count(*) between :c and :d

void
querytemplate_t::set_13() {
  _aggr_having       = kAggrMin; // having min(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having min(l_quantity) between
  _pred_having_cnt   = kPredRg;  // and count(*) between :c and :d
  _bool_connector    = kAnd;     // and count(*) between :c and :d
  _pred_where        = kNoPred;  // empty where clause
}


// template 14:
// group by l_orderkey
// having max(l_quantity) = :a

void
querytemplate_t::set_14() {
  _aggr_having       = kAggrMax; // having max(l_quantity) between
  _pred_having_aggr  = kPredEq;  // having max(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kNoPred;  // empty where clause
}

// template 15:
// group by l_orderkey
// having max(l_quantity) = :a and count(*) between :c and :d

void
querytemplate_t::set_15() {
  _aggr_having       = kAggrMax; // having max(l_quantity) between
  _pred_having_aggr  = kPredEq;  // having max(l_quantity) between
  _pred_having_cnt   = kPredRg;  // and count(*) between :c and :d
  _bool_connector    = kAnd;     // and count(*) between :c and :d
  _pred_where        = kNoPred;  // empty where clause
}


// template 16:
// group by l_orderkey
// having max(l_quantity) between :a and :b

void
querytemplate_t::set_16() {
  _aggr_having       = kAggrMax; // having max(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having max(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kNoPred;  // empty where clause
}

// template 17:
// group by l_orderkey
// having max(l_quantity) between :a and :b and count(*) between :c and :d

void
querytemplate_t::set_17() {
  _aggr_having       = kAggrMax; // having max(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having max(l_quantity) between
  _pred_having_cnt   = kPredRg;  // and count(*) between :c and :d
  _bool_connector    = kAnd;     // and count(*) between :c and :d
  _pred_where        = kNoPred;  // empty where clause
}

// template 18:
// group by l_orderkey
// having count(*) = :a

void
querytemplate_t::set_18() {
  _aggr_having       = kAggrCnt; // having count(l_quantity) =
  _pred_having_aggr  = kPredEq;  // having count(l_quantity) =
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kNoPred;  // empty where clause
}


// template 19:
// group by l_orderkey
// having count(*) between :a and :b

void
querytemplate_t::set_19() {
  _aggr_having       = kAggrCnt; // having count(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having count(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kNoPred;  // empty where clause
}


// template 20:
// group by l_orderkey
// having sum(l_quantity) = :a

void
querytemplate_t::set_20() {
  _aggr_having       = kAggrSum; // having sum(l_quantity) =
  _pred_having_aggr  = kPredEq;  // having sum(l_quantity) =
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kNoPred;  // empty where clause
}

// template 21:
// group by l_orderkey
// having sum(l_quantity) = :a and count(*) between :c and :d

void
querytemplate_t::set_21() {
  _aggr_having       = kAggrSum; // having sum(l_quantity) =
  _pred_having_aggr  = kPredEq;  // having sum(l_quantity) =
  _pred_having_cnt   = kPredRg;  // count(*) between :c and :d
  _bool_connector    = kAnd;     // count(*) between :c and :d
  _pred_where        = kNoPred;  // empty where clause
}

// template 22:
// group by l_orderkey
// having sum(l_quantity) between :a and :b

void
querytemplate_t::set_22() {
  _aggr_having       = kAggrSum; // having sum(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having sum(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kNoPred;  // empty where clause
}


// template 23:
// group by l_orderkey
// having sum(l_quantity) between :a and :b and count(*) between :c and :d

void
querytemplate_t::set_23() {
  _aggr_having       = kAggrSum; // having sum(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having sum(l_quantity) between
  _pred_having_cnt   = kPredRg;  // count(*) between :c and :d
  _bool_connector    = kAnd;     // count(*) between :c and :d
  _pred_where        = kNoPred;  // empty where clause
}


// template 24:
// group by l_orderkey
// having avg(l_quantity) between :a and :b

void
querytemplate_t::set_24() {
  _aggr_having       = kAggrAvg; // having avg(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having avg(l_quantity) between
  _pred_having_cnt   = kNoPred;  // no additional count predicate
  _bool_connector    = kNoAndOr; // no additional count predicate
  _pred_where        = kNoPred;  // empty where clause
}


// template 25:
// group by l_orderkey
// having avg(l_quantity) between :a and :b and count(*) between :c and :d

void
querytemplate_t::set_25() {
  _aggr_having       = kAggrAvg; // having avg(l_quantity) between
  _pred_having_aggr  = kPredRg;  // having avg(l_quantity) between
  _pred_having_cnt   = kPredRg;  // count(*) between :c and :d
  _bool_connector    = kAnd;     // count(*) between :c and :d
  _pred_where        = kNoPred;  // empty where clause
}




std::string
querytemplate_t::to_string() const {
  std::string lRes;
  lRes += to_string_aggr(aggr_having());
  lRes += to_string_pred(pred_having_aggr());
  if(kNoAndOr != bool_connector()) {
    lRes += to_string_bool_connector(bool_connector());
    lRes += "_cnt";
    lRes += to_string_pred(pred_having_cnt());
  }
  lRes += to_string_pred_where(pred_where());
  return lRes;
}

queryinstance_t::queryinstance_t() : _params() {
  for(uint i = 0; i < 7; ++i) {
    _params[i] = 0;
  }
}

std::ostream&
queryinstance_t::print_param(std::ostream& os, const querytemplate_t& x, const int w) const {
  os << std::setw(w) << a();
  if(kPredRg == x.pred_having_aggr()) {
    os << ' ' << std::setw(w) << b();
  }
  if(kNoAndOr != x.bool_connector()) {
    os << ' ' << std::setw(w) << c() << ' ' << std::setw(w) << d();
  }
  if(kNoPred != x.pred_where()) {
    os << ' ' << std::setw(w) << e() << ' ' << std::setw(w) << f();
  }
  return os;
}

std::ostream&
queryinstance_t::print_query(std::ostream& os, const querytemplate_t& x) const {
  print_select(os, x);
  os << "from (select l_orderkey" << std::endl
     << "      from   Lineitem" << std::endl;
  print_where(os, x);
  os << "      group by l_orderkey" << std::endl;
  print_having(os, x);
  return os;
}

std::ostream&
queryinstance_t::print_select(std::ostream& os, const querytemplate_t& x) const {
  os << "select ";
  os << a();
  if(kPredRg == x.pred_having_aggr()) {
    os << ", " << b();
  }
  if(kNoAndOr != x.bool_connector()) {
    os << ", " << c();
    if(kPredRg == x.pred_having_cnt()) {
      os << ", " << d();
    }
  }
  if(kNoPred != x.pred_where()) {
    os << ", " << e() << ", " << f();
  }
  os << ", count(*)";
  os << std::endl;
  return os;
}

std::ostream&
queryinstance_t::print_where(std::ostream& os, const querytemplate_t& x) const {
  if(kNoPred != x.pred_where()) {
    os << "      where l_suppkey % " << e() << " < " << f() << std::endl;
  }
  return os;
}

std::ostream&
queryinstance_t::print_having(std::ostream& os, const querytemplate_t& x) const {
  os << "      having " 
     << to_string_aggr(x.aggr_having()) << "(l_quantity)";
  if(kPredEq == x.pred_having_aggr()) {
    os << " = " << a();
  } else {
    os << " between " << a() << " and " << b();
  }
  if(kNoAndOr != x.bool_connector()) {
    if(kAnd == x.bool_connector()) {
      os << " and ";
    } else {
      os << " or ";
    }
    if(kPredEq == x.pred_having_cnt()) {
      os << " count(*) = " << c();
    } else {
      os << " count(*) between " << c() << " and " << d();
    }
  }
  os << ");" << std::endl << std::endl;
  return os;
}

std::string 
to_string_aggr(const uint a) {
  std::string lRes;
  switch(a) {
    case kAggrCnt: lRes = "count"; break;
    case kAggrSum: lRes = "sum";   break;
    case kAggrAvg: lRes = "avg";   break;
    case kAggrMin: lRes = "min";   break;
    case kAggrMax: lRes = "max";   break;
    default: lRes = "aggr";
  }
  return lRes;
}

std::string 
to_string_pred(const uint a) {
  std::string lRes;
  switch(a) {
    case kPredEq: lRes = "_eq"; break;
    case kPredRg: lRes = "_rg";  break;
  }
  return lRes;
}

std::string 
to_string_bool_connector(const uint a) {
  std::string lRes;
  switch(a) {
    case kAnd: lRes = "_and"; break;
    case kOr:  lRes = "_or";  break; 
  }
  return lRes;
}

std::string 
to_string_pred_where(const uint a) {
  std::string lRes;
  if(kNoPred != a) {
    lRes = "_sel";
  }
  return lRes;
}


