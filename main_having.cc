#include "util.hh"
#include "variance.hh"
#include <fstream>
#include <cstring>

#include "arg.hh"
#include "cb.hh"

#include "loop.hh"
#include "SimpleProfileGbH.hh"
#include "TpchEst.hh"
#include "EstHavingSqlServer.hh"
#include "EstHavingFentSn.hh"
#include "EstHavingBeta.hh"

/*
 * estimate result cardinality of
 * select ..
 * from R
 * group by R.A
 * having sum(B) in [.,.] 
 * and some extensions
 *
 * see papers/CardEstHaving
 */


/*
 * compare two different possibilities to calculate the number of integer compositions
 */

void
test0() {
  std::cout << "=== test0 ===" << std::endl;
  std::vector<uint64_t> N = {1, 17, 19, 25, 50, 77, 88, 99};
  const uint64_t m = 50;
  std::cout << "k   n       #comp1       #comp2" << std::endl;
  for(uint k = 1; k < 4; ++k) {
    for(auto n : N) {
       const uint64_t lNoDecomp1 = no_comp_1(k, n, m);
       const uint64_t lNoDecomp2 = no_comp_2(k, n, m);
       std::cout << k
                 << ' ' << std::setw( 3) << n
                 << ' ' << std::setw(12) << lNoDecomp1
                 << ' ' << std::setw(12) << lNoDecomp2
                 << std::endl;
    }
  }
  std::cout << std::endl;
}

/*
 *  free functions to produce estimates for
 *  having sum(B) = n
 *  p = p_B = 1/d_B
 *  B in [1,m], m upper bound for B values
 *  F frequencies F_c of count(*) as C values, either estimated or true values
 *  using method 4.4 equation 24
 *  they use different implementations of calculating the number of integer compositions
 *  these estimators are also implemented in SimpleProfileGbh
 */

double
f_est_sum_eq_1(const double p, const uint n, const uint m, const uint_vt& F) {
  double lRes = 0;
  for(uint k = 1; k < F.size(); ++k) {
    lRes += std::pow(p, k) * no_comp_1(k, n, m) * F[k];
  }
  return lRes;
}

double
f_est_sum_eq_2(const double p, const uint n, const uint m, const uint_vt& F) {
  assert(n <= m);
  double lRes = 0;
  for(uint k = 1; k < F.size(); ++k) {
    lRes += std::pow(p, k) * no_comp_2(k, n, m) * F[k];
  }
  return lRes;
}

double
f_est_sum_eq_3(const double p, const uint n, const uint m, const uint_vt& F) {
  double lRes = 0;
  for(uint k = 1; k < F.size(); ++k) {
    lRes += std::pow(p, k) * f_int_comp_csd_rec(k, 1, m, n) * F[k];
  }
  return lRes;
}

/*
 * produce some estimates for Q18'
 * having sum(l_orderkey) [=<|>] <const>
 * estimation method: via integer composition
 * (see Section 4.4 Integer B: Counting Integer Compositions)
 * todo: integrate into SimpleProfile as a new estimation method
 */

void
test1(const TpchEst& aTpchEst) {
  std::cout << "=== test1 ===" << std::endl;
  std::cout << "sum(l_orderkey) [=<|>] b" << std::endl;
  std::cout << "Estimation method: integer compositions" << std::endl;
  bool lEstFreq = true; // use estimate under uda instead of true values
  // table with true frequencies (overwritten if F_c is estimated)
  // #lineitem/order:    1      2      3        4       5       6        7
  uint_vt lFreq = {0, 214172, 214434, 214379, 213728, 214217, 214449, 214621};
  if(lEstFreq) {
    for(uint i = 1; i <= 7; ++i) {
       lFreq[i] = double{1500000} / double{7}; // nodv(l_orderkey) / nodv(C)
    }
  }
  const double p = 1.0/50.0; // probability for each l_quantity to occur
  std::cout << " b   true card        est1     q-error   [est2]" << std::endl;
  for(uint b = 5; b <= 50; b += 5) {
    uint64_t lTrue = aTpchEst.get_true_card_sum_eq(b);
    uint64_t lEst1 = ((uint64_t) std::round(f_est_sum_eq_1(p, b, 50, lFreq)));
    uint64_t lEst2 = ((uint64_t) std::round(f_est_sum_eq_2(p, b, 50, lFreq)));
    std::cout << std::setw(2) << b
              << "  " << std::setw(10) << lTrue
              << "  " << std::setw(10) << lEst1
              << "  " << std::setw(10) << q::qerror<double>(lEst1, lTrue);
    if(lEst2 != lEst1) {
      std::cout << "  " << std::setw(8) << lEst2 << " <DIFFER>";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  double lRes = 0;
  for(uint n = 1; n <= 50; ++n) {
    lRes += f_est_sum_eq_2(p, n, 50, lFreq);
  }
  std::cout << "having sum(l_orderkey) <=  50: " << ((uint64_t) std::round(lRes)) << std::endl;

  lRes = 0;
  for(uint n = 1; n < 50; ++n) {
    lRes += f_est_sum_eq_2(p, n, 50, lFreq);
  }
  std::cout << "having sum(l_orderkey) <   50: " << ((uint64_t) std::round(lRes)) << std::endl;

  lRes = 0;
  for(int n = 301; n <= 350; ++n) {
    lRes += f_est_sum_eq_3(p, n, 50, lFreq);
  }
  std::cout << "having sum(l_orderkey) >  300: " << ((uint64_t) std::round(lRes)) << std::endl;
}

/*
 *  produce estimates using TpchEst
 */

double
run_loop_sum_eq_a(const TpchEst& aTpchEst, 
                  const uint aBegin, const uint aEnd, const uint aStep, 
                  const double aQErrLim, const bool aTrace) {
  double lQErrMax = 0;
  if(aTrace) {
    std::cout << "estimates for equality queries:" << std::endl;
    std::cout << "  b  estimate true-card   q-error" << std::endl;
  }
  for(uint n = aBegin; n < aEnd; n += aStep) {
     const double lEst = aTpchEst.estimate_sum_eq(n);
     const double lTru = aTpchEst.get_true_card_sum_eq(n);
     const double lQEr = q::qerror(std::max<double>(1,lEst), std::max<double>(1, lTru));
     lQErrMax = std::max<double>(lQEr, lQErrMax);
     if(aTrace && (aQErrLim < lQEr)) {
       std::cout << std::setw(3) << n
                 << "  " << std::setw(8) << ((uint64_t) std::round(lEst))
                 << "  " << std::setw(8) << ((uint64_t) lTru)
                 << "  " << std::setw(8) << mt::roundXXt<double>(lQEr)
                 << std::endl;
    }
  }
  if(aTrace) {
    std::cout << "max-qerror: " << lQErrMax << std::endl;
    std::cout << std::endl;
  }
  return lQErrMax;
}

/*
 * produce probability table z_probtable.dat
 * and some estimates using TpchEst and compare
 * them with the true results, also produced by TpchEst
 */

void
test2(const bool aProduceProbTable) {
  std::cout << "=== test2 ===" << std::endl;

  // calculate mean/variance of l_quantity
  Variance<double> lVar;
  lVar.init();
  for(uint i = 1; i <= 50; ++i) {
    lVar.step(i);
  }
  lVar.fin();
  std::cout << "l_quantity in [1,50]:"
            << "  mean = " << lVar.mean() 
            << ", var = " << lVar.variance() 
            << ", stddev = " << lVar.standardDeviation() 
            << std::endl << std::endl;

  TpchEst lTpchEst;

  if(aProduceProbTable) {
    std::ofstream lOsPt("z_probtable.dat");
    lOsPt << "#Probability table:" << std::endl;
    lTpchEst.print_a(lOsPt);
  }

  std::cout << "checking sum(l_quantity) = b" << std::endl;
  std::cout << "----------------------------" << std::endl;
  std::cout << "max-qerr in loop b in 1-199: "
            << run_loop_sum_eq_a(lTpchEst, 1, 200, 1, 7.77, false)
            << std::endl;
  std::cout << "max-qerr in loop b in 200-249: "
            << run_loop_sum_eq_a(lTpchEst, 200, 250, 1, 7.77, false)
            << std::endl;
  std::cout << "max-qerr in loop b in 250-300: "
            << run_loop_sum_eq_a(lTpchEst, 250, 299, 1, 7.77, false)
            << std::endl;

  const double lMaxQErrAllowed = 1.5;
  std::cout << "cases where q-error > " << lMaxQErrAllowed << " for b in [1,350]:" << std::endl;
  run_loop_sum_eq_a(lTpchEst, 1, 351, 1, true, lMaxQErrAllowed);

  double lEst = 0, lTru = 0;

  std::cout << "two half-open range queries:" << std::endl;
  lEst = lTpchEst.estimate_sum_between(0, 50);
  lTru = lTpchEst.get_true_card_sum_between(0, 50);
  std::cout << "having sum(l_orderkey) <= 50: " 
            << "  " << std::setw(8) << ((uint64_t) std::round(lEst)) 
            << "  " << std::setw(8) << ((uint64_t) std::round(lTru)) 
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst, lTru))
            << std::endl;

  lEst = lTpchEst.estimate_sum_between(301, 350);
  lTru = lTpchEst.get_true_card_sum_between(301, 350);
  std::cout << "having sum(l_orderkey) > 300: " 
            << "  " << std::setw(8) << ((uint64_t) std::round(lEst)) 
            << "  " << std::setw(8) << ((uint64_t) std::round(lTru)) 
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst, lTru))
            << std::endl;


  lEst = lTpchEst.estimate_sum_between(1, 19);
  std::cout << "est[sum(qty) < 20, ----] = " << std::setw(8) << (uint64_t) lEst << std::endl;
  lEst = lTpchEst.estimate_sum_between_and_cnt_between(1, 19, 1, 4);
  std::cout << "est[sum(qty) < 20, 1..4] = " << std::setw(8) << (uint64_t) lEst << std::endl;
  lEst = lTpchEst.estimate_sum_between_and_cnt_between(1, 19, 4, 7);
  std::cout << "est[sum(qty) < 20, 4..7] = " << std::setw(8) << (uint64_t) lEst << std::endl;



  lEst = lTpchEst.estimate_sum_between(1, 49);
  std::cout << "est[nosel  , sum(qty) < 50] = " << std::setw(8) << (uint64_t) lEst << std::endl;
  lEst = lTpchEst.estimate_sum_between_sel(1, 49, 0.5);
  std::cout << "est[sel=0.5, sum(qty) < 50] = " << std::setw(8) << (uint64_t) lEst << std::endl;

  lEst = lTpchEst.estimate_sum_between(51, 350);
  std::cout << "est[nosel  , sum(qty) > 50] = " << std::setw(8) << (uint64_t) lEst << std::endl;
  lEst = lTpchEst.estimate_sum_between_sel(51, 350, 0.5);
  std::cout << "est[sel=0.5, sum(qty) > 50] = " << std::setw(8) << (uint64_t) lEst << std::endl;
  std::cout << std::endl;
}

/*
 *  check whether variance of sum is variance of
 *  the variables in the sum times the number of variables
 *  [see central limit theorem]
 *  here: calculate the variance of the l_quantity values
 *  test2 gives observed mean and variance for sums
 */
void
test3() {
  std::cout << "=== test3 ===" << std::endl;
  Variance<double> lVar;
  lVar.init();
  for(uint i = 1; i <= 50; ++i) {
    lVar.step(i);
  }
  lVar.fin();
  std::cout << "[1,50]: mean = " << lVar.mean() 
            << ", var = " << lVar.variance() 
            << ", stddef = " << lVar.standardDeviation() 
            << std::endl;
  std::cout << std::endl;
}

/*
 *
 *
 */

void
mk_table_count_sqlserver(const EstHavingSqlServer& aSqlServerEst, const TpchEst& aTpchEst) {
  std::cout << "=== mk_table_count_sqlserver ===" << std::endl;
  std::cout << "having count(*) ..." << std::endl;
  std::cout << "estimation method: SqlServer normal distribution" << std::endl;
  const int w = 12;
  std::cout << " c"
            << " & " << std::setw(w) << "true card"
            // << " & " << std::setw(w) << "selectivity"
            << " & " << std::setw(w) << "estimate"
            << " & " << std::setw(w) << "q-error"
            << std::endl;
  for(uint k = 1; k <=10; ++k) {
    const uint64_t lTrue = aTpchEst.get_true_card_cnt_eq(k);
    // const double   lSeli = aSqlServerEst.selectivity_cnt_eq(k);
    const uint64_t lEst1 = std::round(aSqlServerEst.estimate_cnt_eq(k));
    const double   lQErr = q::qerror<double>(lEst1, lTrue);
    std::cout << std::setw(2) << k
              << " & " << std::setw(w) << lTrue
              // << " & " << std::setw(w) << lSeli
              << " & " << std::setw(w) << lEst1
              << " & " << std::setw(w) << mt::roundXXt(lQErr)
              << "  \\\\" << std::endl;
   
  }
}


/*
 *  struct testcase_eq_t
 *  for test cases with selectivity
 *  which cannot be derived from TpchEst
 *  and are thus hardcoded here
 *  the query being:
 *  select ...
 *  from   Lineitem
 *  where l_suppkey % :_mod = 0 
 *  group by l_orderkey
 *  having aggr(B) = :_val
 */

struct testcase_eq_t {
  uint _mod;
  uint _val;
  uint _card; // true result cardinality;
  inline uint   mod()  const { return _mod; }
  inline uint   val()  const { return _val; }
  inline uint   card() const { return _card; }
  inline double sel()  const { return (((double) 1) / ((double) mod())); }
};
typedef std::vector<testcase_eq_t> testcase_eq_vt;

struct testcase_between_t {
  uint _mod;
  uint _lb;
  uint _ub;
  uint _card; // true result cardinality;
  inline uint   mod()  const { return _mod; }
  inline uint   lb()   const { return _lb; }
  inline uint   ub()   const { return _ub; }
  inline uint   card() const { return _card; }
  inline double sel()  const { return (((double) 1) / ((double) mod())); }
};
typedef std::vector<testcase_between_t> testcase_between_vt;

/*
 * compare estimates produced using SimpleProfileGbH
 * with those produced by TpchEst
 */

// check having count(*) ...
void
test4(const SimpleProfileGbH&   aSimpEst, 
      const TpchEst&            aTpchEst, 
      const EstHavingSqlServer& aSqlServerEst,
      const EstHavingFentSn&    aEstFentSn,
      const EstHavingBeta       aEstBeta,
      const Cb& aCb) {
  std::cout << "=== test4 ===" << std::endl;
  std::cout << "|-----------------------------|" << std::endl
            << "|q-errors for selected queries|" << std::endl
            << "|-----------------------------|" << std::endl
            << std::endl;
  std::cout << "having count(*) = b:" << std::endl;
  std::cout << "---------------------------" << std::endl;
  double lTru  = 0;
  double lEst[4];
  double lQEr[4];
  /*
  std::cout << "  c true-card"
            << "  SqlServer  q-error" 
            << "  Fent       q-error" 
            << "  Beta       q-error" 
            << "  SimProfile q-error" 
            << std::endl;
  */
  std::cout << "  c"
            << "  SqlServer " 
            << "  Fent      " 
            << "  Beta      " 
            << "  SimProfile" 
            << std::endl;

  for(uint lCnt = 0; lCnt <= 8; ++lCnt) {
    lTru = aTpchEst.get_true_card_cnt_eq(lCnt);

    lEst[0] = aSqlServerEst.estimate_cnt_eq(lCnt);
    lEst[1] = aEstFentSn.estimate_cnt_eq(lCnt);
    lEst[2] = aEstBeta.estimate_cnt_eq(lCnt);
    lEst[3] = aSimpEst.estimate_cnt_eq(lCnt);

    lQEr[0] = q::qerror(lEst[0], lTru);
    lQEr[1] = q::qerror(lEst[1], lTru);
    lQEr[2] = q::qerror(lEst[2], lTru);
    lQEr[3] = q::qerror(lEst[3], lTru);

    std::cout << std::setw(3) << lCnt
              // << "  " << std::setw(8) << ((uint64_t) lTru)

              // << "  " << std::setw(9) << ((uint64_t) std::round(lEst[0]))
              << "  " << std::setw(7) << mt::roundXXXt<double>(lQEr[0])     // SqlServer

              // << "  " << std::setw(9) << ((uint64_t) std::round(lEst[1]))
              << "  " << std::setw(7) << mt::roundXXXt<double>(lQEr[1])     // Fent

              // << "  " << std::setw(9) << ((uint64_t) std::round(lEst[2]))   // Beta-Distribution
              << "  " << std::setw(7) << mt::roundXXXt<double>(lQEr[2])

              // << "  " << std::setw(9) << ((uint64_t) std::round(lEst[3]))
              << "  " << std::setw(7) << mt::roundXXXt<double>(lQEr[3])     // SimpleProfile

              << std::endl;
  }
  std::cout << std::endl;
  testcase_eq_vt lCases = { { 2, 1, 413450},
                            { 3, 1, 517230},
                            { 4, 1, 542115},
                            {10, 1, 399855},
                            { 2, 4, 155767},
                            { 3, 4,  56765},
                            { 4, 4,  23179},
                            {10, 4,    933} };
  std::cout << "... with selectivity:" << std::endl;
  std::cout << "  c  mod true-card    SimplEst q-error" << std::endl;
  for(uint i = 0; i < lCases.size(); ++i) {
    lTru = lCases[i].card();
    lEst[0] = aSimpEst.estimate_cnt_eq_sel(lCases[i].val(), lCases[i].sel());
    lQEr[0] = q::qerror(lEst[0], lTru);
    std::cout << std::setw(3) << lCases[i].val()
              << "  " << std::setw(3) << lCases[i].mod()
              << "  " << std::setw(8) << ((uint64_t) lTru)
              << "  " << std::setw(8) << ((uint64_t) std::round(lEst[0]))
              << "  " << std::setw(8) << mt::roundXXXt<double>(lQEr[0])
              << std::endl;
  }
  std::cout << std::endl;
}


// check having sum(l_quantity)...
double_vt
run_loop_sum_eq_b(const SimpleProfileGbH& aSimpleEst,
                  const TpchEst&          aTpchEst,
                  const EstHavingFentSn&  aEstFentSn,
                  const EstHavingBeta&    aEstBeta,
                  const uint aBegin, const uint aEnd, const uint aStep, 
                  const bool aTrace, const double aQErrLim = 0) {
  double    lEst[4];
  double    lQEr[4];
  double_vt lQErrMax = {0, 0, 0, 0};
  if(aTrace) {
    std::cout << "estimates for equality queries with q-error > " << aQErrLim << ':' << std::endl;
    std::cout << "  b true-card"
              << "   TpchEst   q-error"
              << "  SimplEst   q-error"
              << " FentSnEst   q-error" 
              << "      Beta   q-error" 
              << std::endl;
  }
  for(uint n = aBegin; n < aEnd; n += aStep) {
     const double lTru = aTpchEst.get_true_card_sum_eq(n);
     lEst[0] = aTpchEst.estimate_sum_eq(n);
     lEst[1] = aSimpleEst.estimate_sum_eq(n);
     lEst[2] = aEstFentSn.estimate_sum_eq(n);
     lEst[3] = aEstBeta.estimate_sum_eq(n);
     lQEr[0] = q::qerror(std::max<double>(1,lEst[0]), std::max<double>(1, lTru));
     lQEr[1] = q::qerror(std::max<double>(1,lEst[1]), std::max<double>(1, lTru));
     lQEr[2] = q::qerror(std::max<double>(1,lEst[2]), std::max<double>(1, lTru));
     lQEr[3] = q::qerror(std::max<double>(1,lEst[3]), std::max<double>(1, lTru));
     lQErrMax[0] = std::max<double>(lQEr[0], lQErrMax[0]);
     lQErrMax[1] = std::max<double>(lQEr[1], lQErrMax[1]);
     lQErrMax[2] = std::max<double>(lQEr[2], lQErrMax[2]);
     lQErrMax[3] = std::max<double>(lQEr[3], lQErrMax[3]);
     if(aTrace && (aQErrLim < lQEr[1])) {
       std::cout << std::setw(3) << n
                 << "  " << std::setw(8) << ((uint64_t) lTru)
                 << "  " << std::setw(8) << ((uint64_t) std::round(lEst[0]))
                 << "  " << std::setw(8) << mt::roundXXt<double>(lQEr[0])
                 << "  " << std::setw(8) << ((uint64_t) std::round(lEst[1]))
                 << "  " << std::setw(8) << mt::roundXXt<double>(lQEr[1])
                 << "  " << std::setw(8) << ((uint64_t) std::round(lEst[2]))
                 << "  " << std::setw(8) << mt::roundXXt<double>(lQEr[2])
                 << "  " << std::setw(8) << ((uint64_t) std::round(lEst[3]))
                 << "  " << std::setw(8) << mt::roundXXt<double>(lQEr[3])
                 << std::endl;
    }
  }
  if(aTrace) {
    std::cout << "max-qerror: " 
              << mt::roundXXt(lQErrMax[0]) << ' '
              << mt::roundXXt(lQErrMax[1]) << ' '
              << mt::roundXXt(lQErrMax[2]) << ' '
              << mt::roundXXt(lQErrMax[3])
              << std::endl
              << std::endl;
  }
  return lQErrMax;
}

// check having sum(l_quantity)...
double_vt
run_loop_sum_eq_c(const SimpleProfileGbH& aSimpleEst,
                  const TpchEst&          aTpchEst,
                  const EstHavingFentSn&  aEstFentSn,
                  const EstHavingBeta&    aEstBeta,
                  const uint aBegin, const uint aEnd, const uint aStep,
                  const double aQErrLim, const bool aTrace) {
  const uint lNoEst = 7;
  double    lEst[lNoEst];
  double    lQEr[lNoEst];
  double_vt lQErrMax = {0, 0, 0, 0, 0, 0, 0};
  assert(lQErrMax.size() == lNoEst);
  if(aTrace) {
    std::cout << "estimates for equality queries with q-error > " << aQErrLim << ':' << std::endl;
    std::cout << "  b true-card"
              << "      Fent   q-error"
              << "      Beta   q-error"
              << "     SP(1)   q-error"
              << "     SP(1)   q-error" 
              << "     SP(2)   q-error" 
              << "     SP(3)   q-error" 
              << "        IC   q-error"
              << std::endl;
  }
  for(uint n = aBegin; n < aEnd; n += aStep) {
     const double lTru = aTpchEst.get_true_card_sum_eq(n);
     lEst[0] = aEstFentSn.estimate_sum_eq(n);
     lEst[1] = aEstBeta.estimate_sum_eq(n);
     lEst[2] = aSimpleEst.estimate_sum_eq(n);
     lEst[3] = aSimpleEst.estimate_sum_eq_ic_lim(n, 1);
     lEst[4] = aSimpleEst.estimate_sum_eq_ic_lim(n, 2);
     lEst[5] = aSimpleEst.estimate_sum_eq_ic_lim(n, 3);
     lEst[6] = aSimpleEst.estimate_sum_eq_ic(n);
     for(uint i = 0; i < lNoEst; ++i) {
       lQEr[i] = q::qerror(std::max<double>(1,lEst[i]), std::max<double>(1, lTru));
       lQErrMax[i] = std::max<double>(lQEr[i], lQErrMax[i]);
     }
     if(aTrace && (aQErrLim < lQEr[1])) {
       std::cout << std::setw(3) << n << "  " << std::setw(8) << ((uint64_t) lTru);
       for(uint i = 0; i < lNoEst; ++i) {
         std::cout << "  " << std::setw(8) << ((uint64_t) std::round(lEst[i]))
                   << "  " << std::setw(8) << mt::roundXXt<double>(lQEr[i]);
       }
       std::cout << std::endl;
    }
  }
  if(aTrace) {
    std::cout << "max-qerror:";
    for(uint i = 0; i < lNoEst; ++i) {
      std::cout << ' ' << mt::roundXXt(lQErrMax[i]);
    }
    std::cout << std::endl;
    std::cout << std::endl;
  }
  return lQErrMax;
}


struct range_t {
  uint _min;
  uint _max;
  range_t() : _min(0), _max(0) {}
  range_t(const uint aMin, const uint aMax) : _min(aMin), _max(aMax) {}
  inline uint min() const { return _min; }
  inline uint max() const { return _max; }
  std::ostream& print(std::ostream& os, const int aWidth) const;
};

std::ostream&
range_t::print(std::ostream& os, const int aWidth) const {
  os << std::setw(0)
     << '[' << std::setw(aWidth) << min()
     << ',' << std::setw(aWidth) << max()
     << ']';
  return os;
}

std::ostream&
operator<<(std::ostream& os, const range_t& x) {
  return x.print(os, os.width());
}

using range_vt = std::vector<range_t>;

// check having sum() ...
void
test5(const SimpleProfileGbH& aSimpEst, 
      const TpchEst&          aTpchEst, 
      const EstHavingFentSn&  aEstFentSn, 
      const EstHavingBeta     aEstBeta,
      const Cb&               aCb) {
  std::cout << "=== test5 ===" << std::endl;
  // A: group by attribute(s)
  // B: aggregated attribute as in having sum(B) between 13 and 17
  std::cout << "|-----------------------------|" << std::endl
            << "|q-errors for selected queries|" << std::endl
            << "|-----------------------------|" << std::endl
            << std::endl;
  std::cout << "having sum(l_quantity) = b:" << std::endl;
  std::cout << "---------------------------" << std::endl;
  double    lQErr0 = 0;
  double_vt lQErr1;
  lQErr0 = run_loop_sum_eq_a(          aTpchEst,             1, 200, 1, 7.77, false);
/*
  const bool lTraceLoopB = false;
  if(lTraceLoopB) {
    std::cout << "BEGIN TRACE" << std::endl;
  }
  lQErr1 = run_loop_sum_eq_b(aSimpEst, aTpchEst, aEstFentSn, aEstBeta, 1, 200, 1, 7.77, lTraceLoopB);
  if(lTraceLoopB) {
    std::cout << "END TRACE" << std::endl;
  }
*/

  std::cout << std::endl;

  range_vt lRangesB = {{1, 200}, {200, 249}, {250, 300}};


  std::cout << "                               "
            << ' ' << std::setw(8) << "TpchEst1"
            << ' ' << std::setw(8) << "TpchEst2"
            << ' ' << std::setw(8) << "SimpEst"
            << ' ' << std::setw(8) << "Fent"
            << ' ' << std::setw(8) << "Beta"
            << std::endl;

  for(const auto& lRange : lRangesB) {
    lQErr0 = run_loop_sum_eq_a(          aTpchEst,             lRange.min(), lRange.max(), 1, 7.77, false);
    lQErr1 = run_loop_sum_eq_b(aSimpEst, aTpchEst, aEstFentSn, aEstBeta, lRange.min(), lRange.max(), 1, 7.77, false);
    std::cout << "  max-qerr for b in " << std::setw(3) << lRange << ": "
              << ' ' << std::setw(8) << mt::roundXXt(lQErr0)
              << ' ' << std::setw(8) << mt::roundXXt(lQErr1[0])
              << ' ' << std::setw(8) << mt::roundXXt(lQErr1[1])
              << ' ' << std::setw(8) << mt::roundXXt(lQErr1[2])
              << ' ' << std::setw(8) << mt::roundXXt(lQErr1[3])
              << std::endl;
  }

  std::cout << std::endl;
  std::cout << std::endl;

  std::cout << "                              "
            << ' ' << std::setw(8) << "Fent"   // 0
            << ' ' << std::setw(8) << "Beta"   // 1
            << ' ' << std::setw(8) << "SP"     // 2 
            << ' ' << std::setw(8) << "SP(1)"  // 3 // must be equal to SP
            << ' ' << std::setw(8) << "SP(2)"  // 4
            << ' ' << std::setw(8) << "SP(3)"  // 5
            << ' ' << std::setw(8) << "IC"     // 6
            << std::endl;


  using double_vt = std::vector<double>;
  double_vvt lTable(lRangesB.size());

  // for(const auto& lRange : lRangesB)
  for(uint i = 0; i < lRangesB.size(); ++i) {
    const range_t& lRange = lRangesB[i];
    lQErr1 = run_loop_sum_eq_c(aSimpEst, aTpchEst, aEstFentSn, aEstBeta, lRange.min(), lRange.max(), 1, 7.77, false);
    lTable[i] = lQErr1;
    std::cout << "  max-qerr for b in " << std::setw(3) << lRange << ":";
    for(uint i = 0; i < lQErr1.size(); ++i) {
      std::cout << ' ' << std::setw(8) << mt::roundXXt(lQErr1[i]);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  using string_vt = std::vector<std::string>;
  const string_vt lNames = {"Fent", "Beta", "SP", "SP(1)", "SP(2)", "SP(3)", "IC"};
  assert(lNames.size() == lQErr1.size());
  
  // print table pivoted
  std::cout << "=== table for paper ===" << std::endl;
  std::cout << std::setw(10) << "estimator";
  for(uint i = 0; i < lRangesB.size(); ++i) {
    const range_t& lRange = lRangesB[i];
    std::cout << " & " << std::setw(3) << lRange;
  }
  std::cout << "  \\\\\\hline" << std::endl;
  for(uint i : {0, 1, 3, 4, 5, 6}) {
    std::cout << std::setw(10) << lNames[i];
    for(uint j = 0; j < lRangesB.size(); ++j) {
      std::cout << " & " << mt::roundXXXt(lTable[j][i]);
    }
    std::cout << "  \\\\" << std::endl;
  }


  std::cout << std::endl;
  std::cout << "=== endtable for paper ===" << std::endl;
  std::cout << std::endl;


  const double lMaxQErrAllowed = 1.5;
  std::cout << "  all cases where q-error > " << lMaxQErrAllowed << " for b in [1,350]:" << std::endl;
  run_loop_sum_eq_b(aSimpEst, aTpchEst, aEstFentSn, aEstBeta, 1, 351, 1, lMaxQErrAllowed, true);


  double    lTrue = 0;
  double_vt lEst  = {0, 0, 0, 0};

  std::cout << "half-open range queries:" << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << std::string(30, ' ')
            << "  " << "true card"
            << "  " << "tpch-est"
            << "  " << "q-err"
            << "  " << "simp-est"
            << "  " << "q-err"
            << "  " << "Fent(SN)"
            << "  " << "q-err"
            << "  " << "    Beta"
            << "  " << "q-err"
            << std::endl;
  lTrue = aTpchEst.get_true_card_sum_between(0, 50);
  lEst[0] = aTpchEst.estimate_sum_between(0, 50);
  lEst[1] = aSimpEst.estimate_sum_between(0, 50);
  lEst[2] = aEstFentSn.estimate_sum_between(0, 50);
  lEst[3] = aEstBeta.estimate_sum_between(0, 50);
  std::cout << "  sum(l_quantity) <=  50      :"
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << ((uint64_t) std::round(lEst[0]))
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[0], lTrue))
            << "  " << std::setw(8) << ((uint64_t) std::round(lEst[1]))
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << "  " << std::setw(8) << ((uint64_t) std::round(lEst[2]))
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[2], lTrue))
            << "  " << std::setw(8) << ((uint64_t) std::round(lEst[3]))
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[3], lTrue))
            << std::endl;

  lTrue = aTpchEst.get_true_card_sum_between(301, 350);
  lEst[0] = aTpchEst.estimate_sum_between(301, 350);
  lEst[1] = aSimpEst.estimate_sum_between(301, 350);
  lEst[2] = aEstFentSn.estimate_sum_between(301, 350);
  lEst[3] = aEstBeta.estimate_sum_between(301, 350);
  std::cout << "  sum(l_quantity) >  300      :"
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << ((uint64_t) std::round(lEst[0]))
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[0], lTrue))
            << "  " << std::setw(8) << ((uint64_t) std::round(lEst[1]))
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << "  " << std::setw(8) << ((uint64_t) std::round(lEst[2]))
            << "  " << std::setw(5) << mt::roundXt(q::qerror(lEst[2], lTrue))
            << "  " << std::setw(8) << ((uint64_t) std::round(lEst[3]))
            << "  " << std::setw(5) << mt::roundXt(q::qerror(lEst[3], lTrue))
            << std::endl;
  std::cout << std::endl;


  std::cout << "... with count restrictions (and):" << std::endl;
  lTrue = aTpchEst.get_true_card_sum_between(1, 19);
  lEst[0] = aTpchEst.estimate_sum_between(1, 19);
  lEst[1] = aSimpEst.estimate_sum_between(1, 19);
  lEst[2] = aEstFentSn.estimate_sum_between(1, 19);
  lEst[3] = aEstBeta.estimate_sum_between(1, 19);
  std::cout << "  sum(l_quantity) <   20, ----:"
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[0]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[0], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[1]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[2]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[2], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[3]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[3], lTrue))
            << std::endl;
  lTrue   = aTpchEst.get_true_card_sum_between_and_cnt_between(1, 19, 1, 4);
  lEst[0] = aTpchEst.estimate_sum_between_and_cnt_between(1, 19, 1, 4);
  lEst[1] = aSimpEst.estimate_sum_between_and_cnt_between(1, 19, 1, 4);
  // lEst[2] = aEstFentSn.estimate_sum_between_and_cnt_between(1, 19, 1, 4); // does not exist
  std::cout << "  sum(l_quantity) <   20, 1..4:"
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[0]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[0], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[1]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << std::endl;
  lTrue   = aTpchEst.get_true_card_sum_between_and_cnt_between(1, 19, 4, 7);
  lEst[0] = aTpchEst.estimate_sum_between_and_cnt_between(1, 19, 4, 7);
  lEst[1] = aSimpEst.estimate_sum_between_and_cnt_between(1, 19, 4, 7);
  // lEst[2] = aEstFentSn.estimate_sum_between_and_cnt_between(1, 19, 4, 7); // does not exist
  std::cout << "  sum(l_quantity) <   20, 4..7:"
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[0]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[0], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[1]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << std::endl;
  std::cout << std::endl;

  std::cout << "... with count restrictions (or):" << std::endl;
  lTrue   = aTpchEst.get_true_card_sum_between(1, 19);
  lEst[1] = aSimpEst.estimate_sum_between(1, 19);
  std::cout << "  sum(l_quantity) <   20, ----:" 
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << "---"
            << "  " << std::setw(5) << "---"
            << "  " << std::setw(8) << (uint64_t) lEst[1]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << std::endl;

  lTrue   = aTpchEst.get_true_card_sum_between_or_cnt_between(1, 19, 1, 4);
  lEst[1] = aSimpEst.estimate_sum_between_or_cnt_between(1, 19, 1, 4);
  std::cout << "  sum(l_quantity) <   20, 1..4:"
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << "---"
            << "  " << std::setw(5) << "---"
            << "  " << std::setw(8) << (uint64_t) lEst[1]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << std::endl;

  lTrue   = aTpchEst.get_true_card_sum_between_or_cnt_between(1, 19, 4, 7);
  lEst[1] = aSimpEst.estimate_sum_between_or_cnt_between(1, 19, 4, 7);
  std::cout << "  sum(l_quantity) <   20, 4..7:"
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << "---"
            << "  " << std::setw(5) << "---"
            << "  " << std::setw(8) << (uint64_t) lEst[1]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << std::endl;



  

  std::cout << std::endl;

  std::cout << "... with selection predicate on Lineitem:" << std::endl;
  lTrue = aTpchEst.get_true_card_sum_between(1, 49); // 351209
  lEst[0] = aTpchEst.estimate_sum_between(1, 49);
  lEst[1] = aSimpEst.estimate_sum_between(1, 49);
  lEst[2] = aEstFentSn.estimate_sum_between(1, 49);
  lEst[3] = aEstBeta.estimate_sum_between(1, 49);
  std::cout << "  sum(l_quantity) <   50, ----:" 
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[0]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[0], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[1]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[2]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[2], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[3]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[3], lTrue))
            << std::endl;
  lTrue = 623204; // must be faked, cannot use aTpchEst
  lEst[0] = aTpchEst.estimate_sum_between_sel(1, 49, 0.5);
  lEst[1] = aSimpEst.estimate_sum_between_sel(1, 49, 0.5);
  // lEst[2] = aEstFentSn.estimate_sum_between_sel(1, 49, 0.5); // does not exist
  std::cout << "  sum(l_quantity) <   50,  0.5:"
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[0]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[0], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[1]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << std::endl;

  lTrue = aTpchEst.get_true_card_sum_between(51, 350); // 1137386
  lEst[0] = aTpchEst.estimate_sum_between(51, 350);
  lEst[1] = aSimpEst.estimate_sum_between(51, 350);
  lEst[2] = aEstFentSn.estimate_sum_between(51, 350);
  lEst[3] = aEstBeta.estimate_sum_between(51, 350);
  std::cout << "  sum(l_quantity) >   50, ----:"
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[0]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[0], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[1]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[2]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[2], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[3]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[3], lTrue))
            << std::endl;
  lTrue = 645251; // must be faked, cannot use aTpchEst
  lEst[0] = aTpchEst.estimate_sum_between_sel(51, 350, 0.5);
  lEst[1] = aSimpEst.estimate_sum_between_sel(51, 350, 0.5);
  // lEst[2] = aEstFentSn.estimate_sum_between_sel(51, 350, 0.5); // does not exist
  std::cout << "  sum(l_quantity) >   50,  0.5:"
            << "  " << std::setw(8) << ((uint64_t) std::round(lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[0]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[0], lTrue))
            << "  " << std::setw(8) << (uint64_t) lEst[1]
            << "  " << std::setw(5) << mt::roundXXXt(q::qerror(lEst[1], lTrue))
            << std::endl;
}


void
test6(const SimpleProfileGbH& aSimpEst, const TpchEst& aTpchEst, const Cb& aCb) {
  std::cout << "=== test6 ===" << std::endl;
  std::cout << "|-----------------------------|" << std::endl
            << "|q-errors for selected queries|" << std::endl
            << "|-----------------------------|" << std::endl
            << std::endl;
  std::cout << "having avg(l_quantity) between lb and ub:" << std::endl;
  std::cout << "---------------------------" << std::endl;

  double lTru = 0;
  double lEst = 0;
  double lQEr = 0;
  std::cout << "lb,ub true-card  simp-est    qerror" << std::endl;
  for(uint lLb = 0; lLb <= 50; lLb += 5) {
    const uint lUb = lLb + 5;
    lTru = aTpchEst.get_true_card_avg_between(lLb, lUb);
    lEst = aSimpEst.estimate_avg_between(lLb, lUb);
    lQEr = q::qerror(lEst, lTru);
    std::cout << std::setw(2) << lLb << ',' << std::setw(2) << lUb
              << "  " << std::setw(8) << ((uint64_t) lTru)
              << "  " << std::setw(8) << ((uint64_t) std::round(lEst))
              << "  " << std::setw(8) << mt::roundXXt<double>(lQEr)
              << std::endl;
  }
}

void
test7a(const SimpleProfileGbH& aSimpEst, const TpchEst& aTpchEst, const Cb& aCb) {
  std::cout << "=== test7a ===" << std::endl;
  testcase_eq_vt lCases {
                 {0,  5, 83067},
                 {0, 10, 55725},
                 {0, 15, 37057},
                 {0, 20, 24829},
                 {0, 25, 17031},
                 {0, 30, 12141},
                 {0, 35,  8981},
                 {0, 40,  6752},
                 {0, 45,  5315},
                 {0, 50,  4441},
              };

  std::cout << "|-----------------------------|" << std::endl
            << "|q-errors for selected queries|" << std::endl
            << "|-----------------------------|" << std::endl
            << std::endl;
  std::cout << "having min(l_quantity) = val" << std::endl;
  std::cout << "----------------------------" << std::endl;

  double lTru = 0;
  double lEst = 0;
  double lQEr = 0;

  std::cout << "val  true-card    SimplEst q-error" << std::endl;
  for(uint i = 0; i < lCases.size(); ++i) {
    lTru = lCases[i].card();
    lEst = aSimpEst.estimate_min_eq(lCases[i].val());
    lQEr = q::qerror(lEst, lTru);
    std::cout << std::setw(3) << lCases[i].val()
              << "  " << std::setw(8) << ((uint64_t) lTru)
              << "  " << std::setw(8) << ((uint64_t) std::round(lEst))
              << "  " << std::setw(8) << mt::roundXXXt<double>(lQEr)
              << std::endl;
  }
}

void          
test7b(const SimpleProfileGbH& aSimpEst, const TpchEst& aTpchEst, const Cb& aCb) {
  std::cout << "=== test7b ===" << std::endl;
  testcase_eq_vt lCases {
                 {0,  1,   4372},
                 {0,  4,   4978},
                 {0,  7,   5641},
                 {0, 10,   6663},
                 {0, 13,   7567},
                 {0, 16,   9066},
                 {0, 19,  10814},
                 {0, 22,  13046},
                 {0, 25,  16006},
                 {0, 28,  19909},
                 {0, 31,  24782},
                 {0, 34,  31464},
                 {0, 37,  40102},
                 {0, 40,  51202},
                 {0, 43,  65035},
                 {0, 46,  83601},
                 {0, 49, 106214}
              };

  std::cout << "|-----------------------------|" << std::endl
            << "|q-errors for selected queries|" << std::endl
            << "|-----------------------------|" << std::endl
            << std::endl;
  std::cout << "having max(l_quantity) = val" << std::endl;
  std::cout << "----------------------------" << std::endl;
    
  double lTru = 0;
  double lEst = 0;
  double lQEr = 0;
    
  std::cout << "val  true-card    SimplEst q-error" << std::endl;
  for(uint i = 0; i < lCases.size(); ++i) {
    lTru = lCases[i].card();
    lEst = aSimpEst.estimate_max_eq(lCases[i].val());
    lQEr = q::qerror(lEst, lTru);
    std::cout << std::setw(3) << lCases[i].val()
              << "  " << std::setw(8) << ((uint64_t) lTru)
              << "  " << std::setw(8) << ((uint64_t) std::round(lEst))
              << "  " << std::setw(8) << mt::roundXXXt<double>(lQEr)
              << std::endl;
  }
}

double
get_estimate(const SimpleProfileGbH& aSp, const queryinstance_t& aInst, const querytemplate_t& aQueryTemplate) {
  double lRes = 0;
  if(kPredEq == aQueryTemplate.pred_having_aggr()) {
    // having aggr(.) = ...
    if(kNoAndOr == aQueryTemplate.bool_connector()) {
      if(kNoPred == aQueryTemplate.pred_where()) {
        switch(aQueryTemplate.aggr_having()) {
          case kAggrCnt: lRes = aSp.estimate_cnt_eq(aInst.a()); break;
          case kAggrSum: lRes = aSp.estimate_sum_eq(aInst.a()); break;
          case kAggrAvg: std::cout << "NYI: avg(.) = " << std::endl; assert(0 == 3); break;
          case kAggrMin: lRes = aSp.estimate_min_eq(aInst.a()); break;
          case kAggrMax: lRes = aSp.estimate_max_eq(aInst.a()); break;
        }
      } else {
        switch(aQueryTemplate.aggr_having()) {
          case kAggrCnt: lRes = aSp.estimate_cnt_eq_sel(aInst.a(), aInst.sel()); break;
          case kAggrSum: lRes = aSp.estimate_sum_eq_sel(aInst.a(), aInst.sel()); break;
          case kAggrAvg: std::cout << "NYI: avg(.) = " << std::endl; assert(0 == 4); break;
          case kAggrMin: lRes = aSp.estimate_min_eq_sel(aInst.a(), aInst.sel()); break;
          case kAggrMax: lRes = aSp.estimate_max_eq_sel(aInst.a(), aInst.sel()); break;
        }
      }
    } else
    if(kAnd == aQueryTemplate.bool_connector()) {
      // ... and count(*) between ...
      if(kNoPred == aQueryTemplate.pred_where()) {
        switch(aQueryTemplate.aggr_having()) {
          case kAggrCnt: assert(0 == 1); break;
          case kAggrSum: lRes = aSp.estimate_sum_eq_cnt(aInst.a(), aInst.c(), aInst.d()); break;
          case kAggrAvg: std::cout << "NYI: avg(.) = " << std::endl; assert(0 == 5); break;
          case kAggrMin: lRes = aSp.estimate_min_eq_cnt(aInst.a(), aInst.c(), aInst.d()); break;
          case kAggrMax: lRes = aSp.estimate_max_eq_cnt(aInst.a(), aInst.c(), aInst.d()); break;
        }
      } else {
        std::cout << "NYI: where p having aggr(.) ... and count(*) between ..." << std::endl;
        assert(0 == 2);
      }
    } else
    if(kOr == aQueryTemplate.bool_connector()) {
      std::cout << "NYI: or" << std::endl;
    } else {
      assert(0 == 7);
    }
  } else {
    if(kNoAndOr == aQueryTemplate.bool_connector()) {
      if(kNoPred == aQueryTemplate.pred_where()) {
        switch(aQueryTemplate.aggr_having()) {
          case kAggrCnt: lRes = aSp.estimate_cnt_between(aInst.a(), aInst.b()); break;
          case kAggrSum: lRes = aSp.estimate_sum_between(aInst.a(), aInst.b()); break;
          case kAggrAvg: lRes = aSp.estimate_avg_between(aInst.a(), aInst.b()); break;
          case kAggrMin: lRes = aSp.estimate_min_between(aInst.a(), aInst.b()); break;
          case kAggrMax: lRes = aSp.estimate_max_between(aInst.a(), aInst.b()); break;
        }
      } else {
        switch(aQueryTemplate.aggr_having()) {
          case kAggrCnt: lRes = aSp.estimate_cnt_between_sel(aInst.a(), aInst.b(), aInst.sel()); break;
          case kAggrSum: lRes = aSp.estimate_sum_between_sel(aInst.a(), aInst.b(), aInst.sel()); break;
          case kAggrAvg: lRes = aSp.estimate_avg_between_sel(aInst.a(), aInst.b(), aInst.sel()); break;
          case kAggrMin: lRes = aSp.estimate_min_between_sel(aInst.a(), aInst.b(), aInst.sel()); break;
          case kAggrMax: lRes = aSp.estimate_max_between_sel(aInst.a(), aInst.b(), aInst.sel()); break;
        }
      }
    } else
    if(kAnd == aQueryTemplate.bool_connector()) {
      // ... and count(*) between ...
      if(kNoPred == aQueryTemplate.pred_where()) {
        switch(aQueryTemplate.aggr_having()) {
          case kAggrCnt: assert(0 == 1); break;
          case kAggrSum: lRes = aSp.estimate_sum_between_and_cnt_between(aInst.a(), aInst.b(), aInst.c(), aInst.d()); break;
          case kAggrAvg: lRes = aSp.estimate_avg_between_and_cnt_between(aInst.a(), aInst.b(), aInst.c(), aInst.d()); break;
          case kAggrMin: lRes = aSp.estimate_min_between_and_cnt_between(aInst.a(), aInst.b(), aInst.c(), aInst.d()); break;
          case kAggrMax: lRes = aSp.estimate_max_between_and_cnt_between(aInst.a(), aInst.b(), aInst.c(), aInst.d()); break;
        }
      } else {
        std::cout << "NYI: where p having aggr(.) ... and count(*) between ..." << std::endl;
        assert(0 == 2);
      }
    } else 
    if(kOr == aQueryTemplate.bool_connector()) {
      std::cout << "NYI: or" << std::endl;
    } else {
      assert(0 == 7);
    }
  }
  return lRes;
}

bool
read_data(double_vt& aCardTru, queryinstance_vt&  aQueries, const querytemplate_t& aTemp) {
  aCardTru.clear();
  aQueries.clear();
  std::string   lFilename("data/" + aTemp.to_string() + ".dat");
  std::ifstream lIs(lFilename);
  if(!lIs) {
    std::cout << "can't open file '" << lFilename << "'." << std::endl;
    return false;
  }

  queryinstance_t lQuery;
  uint            lCard = 0;
  while(!lIs.eof()) {
    lIs >> lQuery.get_a();
    if(kPredRg == aTemp.pred_having_aggr()) {
      lIs >> lQuery.get_b();
    }
    if(kNoAndOr != aTemp.bool_connector()) {
      lIs >> lQuery.get_c();
      if(kPredRg == aTemp.pred_having_cnt()) {
        lIs >> lQuery.get_d();
      }
    }
    if(kNoPred != aTemp.pred_where()) {
      lIs >> lQuery.get_e() >> lQuery.get_f();
    }
    lIs >> lCard;
 
    if(!lIs.eof()) {
      aCardTru.push_back(lCard);
      aQueries.push_back(lQuery);
    }
  }
  return true;
}

double
run_estimates(const SimpleProfileGbH& aSimpEst, 
              const double_vt&        aCardTru, 
              const queryinstance_vt& aQueries, 
              const querytemplate_t&  aTemp,
              const Cb&               aCb) {
  assert(aCardTru.size() == aQueries.size());
  double lQErrMax = 0;
  double lTru     = 0;
  double lEst     = 0;
  double lQEr     = 0;
  if(aCb.trace() || aCb.printBad()) {
      std::cout << " "
                << ' ' << std::setw(3) << 'a';
      if(kPredRg == aTemp.pred_having_aggr()) {
        std::cout << ' ' << std::setw(3) << 'b';
      }
      if(kNoAndOr != aTemp.bool_connector()) {
        std::cout << ' ' << std::setw(3) << 'c' << ' ' << std::setw(3) << 'd';
      }
      if(kNoPred != aTemp.pred_where()) {
        std::cout << ' ' << std::setw(3) << 'e' << ' ' << std::setw(3) << 'f';
      }
      std::cout << "  " << std::setw(8) << "true"
                << "  " << std::setw(8) << "est"
                << "  " << std::setw(8) << "qerr"
                << std::endl;

  }
  for(uint i = 0; i < aCardTru.size(); ++i) {
    lTru = aCardTru[i];
    lEst = get_estimate(aSimpEst, aQueries[i], aTemp);
    lQEr = q::qerror<double>(std::max<double>(1,lEst), std::max<double>(1, lTru));
    if(lQEr > lQErrMax) {
      lQErrMax = lQEr; 
    }
    if(((aCb.theta() < lTru) || (aCb.theta() < lEst)) && (aCb._theta_qerr_max < lQEr)) {
      aCb._theta_qerr_max = lQEr;
    }
    
    if(aCb.trace() || (aCb.printBad() && (aCb.lim_qerr() < lQEr))) {
      std::cout << "  ";
      aQueries[i].print_param(std::cout, aTemp, 3);
      std::cout << "  " << std::setw(8) << ((uint) lTru)
                << "  " << std::setw(8) << ((uint) std::round(lEst))
                << "  " << std::setw(8) << lQEr
                << std::endl;
    }
  }
  return lQErrMax;
}

/*
 *  queries with true cardinalities calculated by sqlite
 */

void
test8(const SimpleProfileGbH& aSimpEst, const TpchEst& aTpchEst, const Cb& aCb) {
  if(aCb.trace()) {
    std::cout << "=== test8 ===" << std::endl;
  }
  
  querytemplate_t  lTemp;
  double_vt        lCardTru;
  queryinstance_vt lQueries;
  const uint       lCaseMin = std::max<uint>(1, aCb.caseNo());
  const uint       lCaseMax = (0 == aCb.caseNo() ? querytemplate_t::no_templates_to_gen() : aCb.caseNo());
  std::cout << std::string(28, ' ')
            << "#query  maxqerr     qxerr (theta)"
            << std::endl;
  for(uint i = lCaseMin; i <= lCaseMax; ++i) {
    lTemp.set(i);
    std::cout << "case " << std::setw(2) << i 
              << "  "    << std::left << std::setw(17) << lTemp.to_string() << std::right
              << ": ";
    if(read_data(lCardTru, lQueries, lTemp)) {
      std::cout << std::setw(6) << lCardTru.size();
      if(aCb.trace() || aCb.printBad()) {
        std::cout << std::endl;
      }
      aCb._theta_qerr_max = 0;
      const double lQErrMax = run_estimates(aSimpEst, lCardTru, lQueries, lTemp, aCb);
      std::cout << std::fixed;
      std::cout.precision(2);
      std::cout << ' ' << std::setw(8) << mt::roundXXt(lQErrMax) 
                << ' ' << std::setw(8) << mt::roundXXt(aCb._theta_qerr_max)
                << " (" << (uint) mt::roundXXt(aCb.theta()) << ')'
                << std::endl;
    } else {
      std::cout << "  failure reading data." << std::endl;
    }
  }
}

/*
 *  queries with true cardinalities from TpchEst
 */

double
get_true_card(const TpchEst& aTpchEst,
     queryinstance_t aInst, const querytemplate_t aQueryTemplate) {
  uint lRes = 0;
  if(kPredEq == aQueryTemplate.pred_having_aggr()) {
    // having aggr(.) = ...
    if(kNoAndOr == aQueryTemplate.bool_connector()) {
      if(kNoPred == aQueryTemplate.pred_where()) {
        switch(aQueryTemplate.aggr_having()) {
          case kAggrCnt: lRes = aTpchEst.get_true_card_cnt_eq(aInst.a()); break;
          case kAggrSum: lRes = aTpchEst.get_true_card_sum_eq(aInst.a()); break;
          case kAggrAvg: assert(0 == 3); break;
          case kAggrMin: assert(0 == 4); break;
          case kAggrMax: assert(0 == 5); break;
        }
      } else {
        // where p
        assert(0 == 1);
      }
    } else
    if(kAnd == aQueryTemplate.bool_connector()) {
      // ... and count(*) between ...
      if(kNoPred == aQueryTemplate.pred_where()) {
        switch(aQueryTemplate.aggr_having()) {
          case kAggrCnt: assert(0 == 1); break;
          case kAggrSum: lRes = aTpchEst.get_true_card_sum_eq_and_cnt_between(aInst.a(), aInst.c(), aInst.d()); break;
          case kAggrAvg: assert(0 == 6); break;
          case kAggrMin: assert(0 == 7); break;
          case kAggrMax: assert(0 == 8); break;
        }
      } else {
        std::cout << "NYI: where p having aggr(.) ... and count(*) between ..." << std::endl;
        assert(0 == 9);
      }
    } else
    if(kOr == aQueryTemplate.bool_connector()) {
      std::cout << "NYI: or" << std::endl;
    } else {
      assert(0 == 10);
    }
  } else {
    // having aggr(.) between
    if(kNoAndOr == aQueryTemplate.bool_connector()) {
      if(kNoPred == aQueryTemplate.pred_where()) {
        switch(aQueryTemplate.aggr_having()) {
          case kAggrCnt: lRes = aTpchEst.get_true_card_cnt_between(aInst.a(), aInst.b()); break;
          case kAggrSum: lRes = aTpchEst.get_true_card_sum_between(aInst.a(), aInst.b()); break;
          case kAggrAvg: lRes = aTpchEst.get_true_card_avg_between(aInst.a(), aInst.b()); break;
          case kAggrMin: assert(0 == 11); break;
          case kAggrMax: assert(0 == 12); break;
        }
      } else {
        assert(0 == 13);
      }
    } else
    if(kAnd == aQueryTemplate.bool_connector()) {
      // ... and count(*) between ...
      if(kNoPred == aQueryTemplate.pred_where()) {
        switch(aQueryTemplate.aggr_having()) {
          case kAggrCnt: assert(0 == 1); break;
          case kAggrSum: lRes = aTpchEst.get_true_card_sum_between_and_cnt_between(aInst.a(), aInst.b(), aInst.c(), aInst.d()); break;
          case kAggrAvg: lRes = aTpchEst.get_true_card_avg_between_cnt(aInst.a(), aInst.b(), aInst.c(), aInst.d()); break;
          case kAggrMin: assert(0 == 14); break;
          case kAggrMax: assert(0 == 14); break;
        }
      } else {
        assert(0 == 15);
      }
    } else
    if(kOr == aQueryTemplate.bool_connector()) {
      assert(0 == 16);
    } else {
      assert(0 == 17);
    }
  }
  return lRes;
}


void
eval_query(const SimpleProfileGbH& aSimpEst,       const TpchEst&         aTpchEst,
                 queryinstance_t&  aQueryInstance, const querytemplate_t& aQueryTemplate, const Cb& aCb) {
  const uint   lTru = get_true_card(aTpchEst, aQueryInstance, aQueryTemplate);
  const double lEst = get_estimate(aSimpEst, aQueryInstance, aQueryTemplate);
  const double lQEr = q::q_error_safe<double>(lEst, lTru);
  ++(aCb._count);
  if(aCb._qerr_max < lQEr) {
    aCb._qerr_max = lQEr;
  }
  if(((aCb.theta() < lTru) || (aCb.theta() < lEst)) && (aCb._theta_qerr_max < lQEr)) {
    aCb._theta_qerr_max = lQEr;
  }
  if(aCb.trace() || (aCb.printBad() && (aCb.lim_qerr() < lQEr))) {
    std::cout << ' ';
    aQueryInstance.print_param(std::cout, aQueryTemplate, 3);
    std::cout << "  " << std::setw(8) << lTru
              << "  " << std::setw(8) << ((uint64_t) std::round(lEst))
              << "  " << std::setw(8) << mt::roundXXt(lQEr)
              << std::endl;
  }
}

void
eval_having_count(const SimpleProfileGbH& aSimpEst, const TpchEst& aTpchEst,
                  queryinstance_t& aQueryInstance, const querytemplate_t& aQueryTemplate, const Cb& aCb) {
  loopspec_t lLoop;
  if(kNoAndOr == aQueryTemplate.bool_connector()) {
    eval_query(aSimpEst, aTpchEst, aQueryInstance, aQueryTemplate, aCb);
    return;
  }
  if(kPredEq == aQueryTemplate.pred_having_cnt()) {
    lLoop.set(1, 7, 1);
  } else
  if(kPredRg == aQueryTemplate.pred_having_cnt()) {
    lLoop.set(1, 7, 1, 2, 5, 1);
  }

  lLoop.eval_estimates(&eval_query, aQueryInstance.get_c(), aQueryInstance.get_d(),
                       aSimpEst, aTpchEst, aQueryInstance, aQueryTemplate, aCb);
}

void
eval_having_aggr(const SimpleProfileGbH& aSimpEst, const TpchEst& aTpchEst,
                  queryinstance_t& aQueryInstance, const querytemplate_t& aQueryTemplate, const Cb& aCb) {
  loopspec_t lLoop;
  if(kAggrCnt == aQueryTemplate.aggr_having()) {
    if(kPredEq == aQueryTemplate.pred_having_aggr()) {
      lLoop.set(1, 7, 1);
    } else
    if(kPredRg == aQueryTemplate.pred_having_aggr()) {
       lLoop.set(1, 7, 1, 2, 5, 1);
    }
  } else
  if(kAggrSum == aQueryTemplate.aggr_having()) {
    if(kPredEq == aQueryTemplate.pred_having_aggr()) {
      lLoop.set(1, 300, 1);
    } else
    if(kPredRg == aQueryTemplate.pred_having_aggr()) {
      lLoop.set(1, 300, 1, 2, 10, 1);
    }
  } else
  if(kAggrAvg == aQueryTemplate.aggr_having()) {
    if(kPredEq == aQueryTemplate.pred_having_aggr()) {
      assert(0 == 1);
    } else {
      lLoop.set(1, 300, 1, 2, 10, 1);
    }
  } else
  if((kAggrMin == aQueryTemplate.aggr_having()) ||
     (kAggrMax == aQueryTemplate.aggr_having())) {
    if(kPredEq == aQueryTemplate.pred_having_aggr()) {
      lLoop.set(1, 50, 1);
    } else
    if(kPredRg == aQueryTemplate.pred_having_aggr()) {
      lLoop.set(1, 50, 1, 2, 10, 1);
    }
  } else {
    assert(0 == 1);
  }

  lLoop.eval_estimates(&eval_having_count, aQueryInstance.get_a(), aQueryInstance.get_b(),
                       aSimpEst, aTpchEst, aQueryInstance, aQueryTemplate, aCb);
}

void
eval_queries(const querytemplate_t& aQueryTemplate, 
             const SimpleProfileGbH& aSimpEst, const TpchEst& aTpchEst, const Cb& aCb) {
  queryinstance_t lQueryInstance;
  eval_having_aggr(aSimpEst, aTpchEst, lQueryInstance, aQueryTemplate, aCb);
}

void
test9(const SimpleProfileGbH& aSimpEst, const TpchEst& aTpchEst, const Cb& aCb) {
  if(aCb.trace()) {
    std::cout << "=== test9 ===" << std::endl;
  }
  const uint lBegin = std::max<uint>(querytemplate_t::no_templates_to_gen() + 1, aCb.caseNo());
  const uint lEnd   = (0 == aCb.caseNo() ? querytemplate_t::no_templates() : aCb.caseNo());
  querytemplate_t lQueryTemplate;
  for(uint i = lBegin; i <= lEnd; ++i) {
    aCb._qerr_max       = 0;
    aCb._theta_qerr_max = 0;
    aCb._count          = 0;
    lQueryTemplate.set(i);
    std::cout << "case " << std::setw(2) << i
              << "  " << std::left << std::setw(17) << lQueryTemplate.to_string() << std::right
              << ": ";
    if(aCb.trace() || aCb.printBad()) {
      std::cout << std::endl;
    }
    eval_queries(lQueryTemplate, aSimpEst, aTpchEst, aCb);
    std::cout << ' ' << std::setw(5)  << aCb._count
              << ' ' << std::setw(8)  << mt::roundXXt(aCb._qerr_max)
              << ' ' << std::setw(8)  << mt::roundXXt(aCb._theta_qerr_max)
              << ' ' << '(' << (uint) mt::roundXXt(aCb.theta()) << ')'
              << std::endl;
  }
}

int
main(const int argc, const char* argv[]) {
  Cb lCb;
  argdesc_vt lArgDesc;
  construct_arg_desc(lArgDesc);

  if(!parse_args<Cb>(1, argc, argv, lArgDesc, lCb, &Cb::pat_or_aggr)) {
    std::cerr << "error while parsing arguments" << std::endl;
    print_usage(std::cerr, argv[0], lArgDesc);
    return -1;
  }

  if(lCb.help()) {
    std::cout << "|--" << std::string(strlen(argv[0]), '-') << "---" << std::endl
              << "|  " << argv[0] << "  |" << std::endl
              << "|--" << std::string(strlen(argv[0]), '-') << "---" << std::endl
              << "| produces estimates and q-errors for different estimation tasks." << std::endl
              << "| The estimates given in the paper are produced by" << std::endl
              << "|  " << argv[0] << " [cnt|sum|avg|min] " << std::endl
              << "| Other estimates can be produced by" << std::endl
              << "|  " << argv[0] << " pat" << std::endl
              << "| Here, there are 25 test cases available, corresponding to 25 different" << std::endl
              << "| estimation procedures in SimpleProfileGbH." << std::endl
              << "| cases 1-17 require prior evaluation of SQL queries generated via main_gen." << std::endl
              << "| For the remaining cases TpchEst is able to give the correct result cardinality." << std::endl
              << "| If only a single test case is to be run, use -n"
              << std::endl;
    print_usage(std::cout, argv[0], lArgDesc);
    return 0;
  }

  // get instance of TpchEst
  TpchEst lTpchEst;

  // test0(); // checks on number of compositions, deprecated (see infra/integer_compositions.hh/cc)
  if(1 == argc) {
    test1(lTpchEst);  // having sum(B) ...; using Sec. 4.4 Integer B: counting integer compositions
    std::cout << std::endl;
  }
  // test2(false);  // check estimates of TpchEst, generate probability table z_probtable.dat
  // test3();       // check Variance<> for sums

  const double lMean = (lCb.calcMeanVar() ? uni_dist_mean(1, 50) : 25.507967136655);
  const double lVari = (lCb.calcMeanVar() ? uni_dist_variance(1, 50) : 208.117016107785);
  const double lStdD = std::sqrt(lVari);
  const double lSkew = 0; // skew of uniform distribution on [a,b] is always zero

  if(lCb.calcMeanVar() && lCb.trace()) {
    std::cout << "calculated mean(l_quantity)     = " << lMean << std::endl;
    std::cout << "calculated variance(l_quantity) = " << lVari << std::endl;
  }

  // relation R: lineitem 
  // A = l_orderkey
  // B = l_quantity
  const uint lCardR = 6'001'215; 
  gm::attr_stat_t lAttrStatA;

  // l_orderkey
  lAttrStatA._min      = 1;
  lAttrStatA._max      = 6'000'000;
  lAttrStatA._nodv     = 1'500'000;
  lAttrStatA._dist_param._mean     = 3000279.604204982;
  lAttrStatA._dist_param._std_dev  = 1732187.8734803419;
  lAttrStatA._dist_param._skewness = -0.00017891030171438376;
  lAttrStatA._is_int   = true;
  lAttrStatA._is_ui_dense = false;

  // count(*)
  gm::attr_stat_t lAttrStatCnt;
  lAttrStatCnt._min      = 1;
  lAttrStatCnt._max      = 7;
  lAttrStatCnt._nodv     = 7;
  lAttrStatCnt._dist_param._mean     = 4.0;
  lAttrStatCnt._dist_param._std_dev  = 2.00;
  lAttrStatCnt._dist_param._skewness = 0;
  lAttrStatCnt._is_int      = true;
  lAttrStatCnt._is_ui_dense = true;


  // l_quantity
  gm::attr_stat_t lAttrStatB;
  lAttrStatB._min      = 1;
  lAttrStatB._max      = 50;
  lAttrStatB._nodv     = 50;
  lAttrStatB._dist_param._mean     = lCb.calcMeanVar() ? lMean : 25.507967136654827;
  lAttrStatB._dist_param._std_dev  = lCb.calcMeanVar() ? lStdD : 14.426262537016854;
  lAttrStatB._dist_param._skewness = lCb.calcMeanVar() ? lSkew : -0.0010518682122090267;
  lAttrStatA._is_int   = true;
  lAttrStatA._is_ui_dense = lCb.isUiDense();

  // instantiate simple profile
  SimpleProfileGbH lSimpleProfile(lCardR,
                                  lAttrStatA,
                                  lAttrStatB,
                                  lAttrStatCnt);

  // get instance of SqlServerEst
  EstHavingSqlServer lSqlServerEst(lCardR,
                                   lAttrStatA,
                                   lAttrStatB,
                                   lAttrStatCnt);

  // get instance of EstFentSn
  EstHavingFentSn lEstFentSn(lCardR,
                             lAttrStatA,
                             lAttrStatB,
                             lAttrStatCnt);

  // gest instance of EstBeta
  EstHavingBeta lEstBeta(lCardR,
                         lAttrStatA,
                         lAttrStatB,
                         lAttrStatCnt);
  // some tracing

  if(lCb.trace()) {
    std::cout << "# main_having"
              << " -i  " << lCb.isUiDense()
              << " -c  " << lCb.calcMeanVar()
              << " -n  " << lCb.caseNo()
              << " -q  " << lCb.lim_qerr()
              << " -theta " << lCb.theta()
              << ' ' << lCb.pat_or_aggr()
              << std::endl;
  }

  if(lCb.trace()) {
    lSimpleProfile.print(std::cout);
  }

  // make the table for SQL-Server estimation procedure only
  if(1 == argc) {
    mk_table_count_sqlserver(lSqlServerEst, lTpchEst);
  }

  // main work, depending on the aggregate function or 'pat'

  if("cnt" == lCb.pat_or_aggr()) {
    test4(lSimpleProfile, lTpchEst, lSqlServerEst, lEstFentSn, lEstBeta, lCb);
  }

  if("sum" == lCb.pat_or_aggr()) {
    test5(lSimpleProfile, lTpchEst, lEstFentSn, lEstBeta, lCb);
  }

  if("avg" == lCb.pat_or_aggr()) {
    test6(lSimpleProfile, lTpchEst, lCb);
  }

  if("min" == lCb.pat_or_aggr()) {
    test7a(lSimpleProfile, lTpchEst, lCb);
  }

  if("max" == lCb.pat_or_aggr()) {
    test7b(lSimpleProfile, lTpchEst, lCb);
  }

  if("pat" == lCb.pat_or_aggr()) {
    if(lCb.caseNo() <= querytemplate_t::no_templates_to_gen()) {
      test8(lSimpleProfile, lTpchEst, lCb);
    }
    if((0 == lCb.caseNo()) || (querytemplate_t::no_templates_to_gen() < lCb.caseNo())) {
      test9(lSimpleProfile, lTpchEst, lCb);
    }
  }

};
