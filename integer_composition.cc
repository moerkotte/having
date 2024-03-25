#include "integer_composition.hh"
#include <cmath>

#define NYI 1

/*
 * functions to calculate integer decompositions
 * a_1 + ... + a_k = n
 * where a_i in [lb,ub]
 * (see also papers/CardEstHaving/intcomp.tex)
 */

// general recurrence, slow for k > 5
uint64_t
f_int_comp_rec(const uint k, const uint lb, const uint ub, const uint n) {
  uint64_t lRes = 0;
  // recursion end
  // for 1 == k the solution, if it exists, is uniquely determined
  if(1 == k) {
    if((lb <= n) && (n <= ub)) {
      return 1;
    } else {
      return 0;
    }
  }

  // some shortcuts
  if((k * lb) > n) { return 0; }
  if((k * ub) < n) { return 0; }

  // can fill in zero's if 0 == lb and 0 == n
  if((0 == n) && (0 == lb) && (0 < k)) {
    return 1;
  }

  if(0 == k) { return 0; }

  for(uint i = lb; i <= std::min<uint>(ub,n); ++i) {
    lRes += f_int_comp_rec(k - 1, lb, ub, n - i);
  }
  return lRes;
}

uint64_t
f_int_comp_csd_rec(const uint k, const uint lb, const uint ub, const uint n) {
  if((k * lb) > n) { return 0; }
  if((k * ub) < n) { return 0; }
  uint64_t lRes = 0;
  if(1 == lb) {
    lRes = f_int_comp_1_m(k, 1, ub, n);
  } else
  if(0 == lb) {
    // TODO
    assert(NYI);
  } else {
    lRes = f_int_comp_1_m(k, 1, ub - (lb - 1), n - (k * (lb - 1)));
  }
  return lRes;
}


// a_i in [1,m], n <= m
uint64_t
f_int_comp_1_m_geq_n(const uint k, const uint lb, const uint ub, const uint n) {
  assert(1 == lb);
  assert(n <= ub);
  return ((uint64_t) std::round(fBinomGamma(n - 1, k - 1)));
}

// a_i in [0,m], n <= m;
uint64_t 
f_int_comp_0_m_geq_n(const uint k, const uint lb, const uint ub, const uint n) {
  assert(0 == lb);
  assert(n <= ub);
  return ((uint64_t) std::round(fBinomGamma(n + k - 1, k - 1)));
}

uint64_t
f_gauss(const uint64_t n) {
  return ((n * (n + 1)) / 2);
}


// number of integer compositions for the special case
// k = 2, a_i in [1,ub]
uint64_t
f_int_comp_k_eq_2_lb_eq_1(const uint ub, const uint n) {
  // constexpr uint64_t lb = 1;
  constexpr uint64_t k = 2;
  if(n < k) { return 0; }
  if(n > (k * ub)) { return 0; }
  if(ub >= n) {
    return (n - 1);
  } 
  return (2 * ub - n + 1);
  // return (ub - (n - ub) + 1);
}

// summation form
uint64_t
f_int_comp_k_eq_3_lb_eq_1_sum(const int ub, const int n) {
  constexpr bool lTrace = false;
  const int m = ub;
  const int lb = 1;
  if constexpr(lTrace) {
    std::cout << "f_int_comp_k_eq_3_lb_eq_1_sum(" << ub << ',' << n << "):" << std::endl;
  }
  if(3 > n)     { return 0; }
  if((3*m) < n) { return 0; }
  int lBegin = std::max<int>(n - m, 2*lb); // ok
  int lEnd   = std::min<int>(n - 1, 2*m);  // ok
  uint64_t lRes = 0;
  if constexpr(lTrace) {
    const int lLimit = std::min<int>(ub, n - 2 * lb); // limit for i
    const int lStart = (n > (2 * ub)) ? std::max<uint>(lb, n - 2*ub) : lb;
    std::cout << "  lBegin     = " << std::setw(3) << lBegin << std::endl;
    std::cout << "  lEnd       = " << std::setw(3) << lEnd   << std::endl;
    std::cout << "  lLimit     = " << std::setw(3) << lLimit << std::endl;
    std::cout << "  lStart     = " << std::setw(3) << lStart << std::endl;
    std::cout << "  n - lLimit = " << std::setw(3) << (n - lLimit) << std::endl;
    std::cout << "  case: " << ( (m < lBegin) ? "low" :
                                ((m >= lEnd)  ? "hi" : "mid"))
              << std::endl;
  }
  for(int np = lBegin; np <= lEnd; ++np) {
    lRes += f_int_comp_k_eq_2_lb_eq_1(ub, np);
  }
  return lRes;
}

// n' in [d,e], m >= e
inline uint64_t
f_int_comp_k_eq_3_case_1(const uint64_t m, // unused for calculations
                         const uint64_t d,
                         const uint64_t e,
                         const uint64_t n) {
  return ( (e - d + 1) * (2*m + 1) - f_gauss(e) + f_gauss(d - 1) );
}

// n' in [d,e], m < d
inline uint64_t
f_int_comp_k_eq_3_case_2(const uint64_t m,
                         const uint64_t d,
                         const uint64_t e,
                         const uint64_t n) {
  return f_gauss(e) - f_gauss(d - 1) - (e - d + 1);
}

// closed form (see papers/EstCardHaving/intcomp.tex)
uint64_t
f_int_comp_k_eq_3_lb_eq_1_csd(const int ub, const int n) {
  constexpr bool lTrace = false;
  const int m = ub;
  const int lb = 1;
  if constexpr(lTrace) {
    std::cout << "f_int_comp_k_eq_3_lb_eq_1_csd(" << ub << ',' << n << "):" << std::endl;
  }
  if(3 > n)     { return 0; }
  if((3*m) < n) { return 0; }
  int lBegin = std::max<int>(n - m, 2*lb); // ok
  int lEnd   = std::min<int>(n - 1, 2*m);  // ok
  uint64_t lRes = 0;
  if constexpr(lTrace) {
    std::cout << "  lBegin     = " << std::setw(3) << lBegin << std::endl;
    std::cout << "  lEnd       = " << std::setw(3) << lEnd   << std::endl;
    std::cout << "  case: " << ( (m < lBegin) ? "low" :
                                ((m >= lEnd)  ? "hi" : "mid"))
              << std::endl;
  }
  if(m < lBegin) {
    // lRes = (lEnd - lBegin + 1) * (2*m + 1) - f_gauss(lEnd) + f_gauss(lBegin - 1);
    lRes = f_int_comp_k_eq_3_case_1(m, lBegin, lEnd, n);
  } else
  if(m >= lEnd) {
    // lRes = f_gauss(lEnd) - f_gauss(lBegin - 1) - (lEnd - lBegin + 1);
    lRes = f_int_comp_k_eq_3_case_2(m, lBegin, lEnd, n);
  } else {
    lRes  = f_int_comp_k_eq_3_case_2(m, lBegin, m, n);
    lRes += f_int_comp_k_eq_3_case_1(m, m + 1,  lEnd,  n);
  }
  return lRes;
}




// a_i in [1,m], m <  n;
uint64_t
f_int_comp_1_m(const uint k, const uint lb, const uint ub, const uint n) {
  // constexpr bool lTrace = false;
  assert(1 == lb);
  assert(lb < ub);
  // assert(ub < n);
  if(n < (k * lb)) { return 0; }
  if(n > (k * ub)) { return 0; }
  uint64_t lRes = 0;
  if(1 == k) {
    return 1;
  } else
  if(2 == k) {
    lRes = f_int_comp_k_eq_2_lb_eq_1(ub, n);
  } else
  if(3 == k) {
    lRes = f_int_comp_k_eq_3_lb_eq_1_csd(ub, n);
  } else {
    if(ub >= n) {
      lRes = f_int_comp_1_m_geq_n(k, lb, ub, n);
    } else {
      for(uint i = lb; i <= ub; ++i) {
        lRes += f_int_comp_1_m(k - 1, lb, ub, n - i);
      }
    }
  }
  return lRes;
}



// generic implementation
uint64_t
f_int_decomp(const uint k, const uint lb, const uint ub, const uint n) {
  return 0;
}

/*
 * functions to calculate integer decompositions less than n
 * a_1 + ... + a_k <= n
 */

// a_i in [0,m], n <= m; [Stanley Vol1 p15]
uint64_t 
f_int_comp_0_m_geq_n_leq(const uint k, const uint lb, const uint ub, const uint n) {
  assert(0 == lb);
  assert(n <= ub);
  return  ((uint64_t) std::round(fBinomGamma(n + k, k)));
}
