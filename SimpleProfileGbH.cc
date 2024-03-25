#include "SimpleProfileGbH.hh"

#include "float_ops.hh"

// good: USE_SEL_YAO_APPROX
#define USE_SEL_A
// #define USE_SEL_YAO_APPROX

// good: ALT_MINMAX_A
#define ALT_MINMAX_A 
// #define ALT_MINMAX_B  // use binomial coefficient for min/max calculation

SimpleProfileGbH::SimpleProfileGbH()
                 : _card_R(0),
                   _attr_stat_A(),
                   _attr_stat_B(),
                   _attr_stat_Cnt() {
}

SimpleProfileGbH::SimpleProfileGbH(const uint64_t      aCardR,
                                   const attr_stat_t& aAttrStatA,
                                   const attr_stat_t& aAttrStatB,
                                   const attr_stat_t& aAttrStatCnt)
                 : _card_R(aCardR),
                   _attr_stat_A(aAttrStatA),
                   _attr_stat_B(aAttrStatB),
                   _attr_stat_Cnt(aAttrStatCnt) {
}


/*
 *  having count(*) ...
 */

double
SimpleProfileGbH::estimate_cnt_eq(const uint aCnt) const {
  if((aCnt < min_cnt()) || (max_cnt() < aCnt)) {
    return 0;
  }
  return freq_tot_k(aCnt);
}

double
SimpleProfileGbH::estimate_cnt_between(const uint aLbCnt, const uint aUbCnt) const {
  const uint lLb = std::max<uint>(min_cnt(), aLbCnt);
  const uint lUb = std::min<uint>(max_cnt(), aUbCnt);
  if(lLb > lUb) {
    return 0;
  }
  return ((lUb - lLb + 1) * freq_tot_k(1)); // argument of freq_tot_k ignored anyway
}

double
SimpleProfileGbH::estimate_cnt_eq_sel(const uint aCnt, const double aSelectivity) const {
  if(max_cnt() < aCnt) {
    return 0;
  }
  double lRes1 = 0;
  for(uint k = aCnt; k <= max_cnt(); ++k) {
    lRes1 +=   fBinomGamma(k, aCnt)
            * std::pow(1 - aSelectivity, k - aCnt)
            * std::pow(    aSelectivity, aCnt)
            * freq_tot_k(aCnt);
  }
  double lRes2 = 0;
  for(uint k = aCnt; k <= max_cnt(); ++k) {
    lRes2 +=   fBinomGamma(k, aCnt) * std::pow(1 - aSelectivity, k - aCnt);
  }
  lRes2 *= freq_tot_k(aCnt) * std::pow(aSelectivity, aCnt);
  assert(is_equal(lRes1, lRes2, 1e-6));
  // TODO: add f_prob_cnt_eq_sel_3 from main_simplify.cc
  return lRes2; 
}

double
SimpleProfileGbH::estimate_cnt_between_sel(const uint aLbCnt, const uint aUbCnt, const double aSelectivity) const {
  const uint lLb = std::max<uint>(min_cnt(), aLbCnt);
  const uint lUb = std::max<uint>(max_cnt(), aUbCnt);
  double lRes = 0;
  for(uint k = lLb; k <= lUb; ++k) {
    lRes += estimate_cnt_eq_sel(k, aSelectivity);
  }
  return lRes;
}

/*
 * sum
 */

double
SimpleProfileGbH::estimate_sum_eq(const double aSumB) const {
  double lRes = 0;
  if((min_cnt() <= 1) && (1 <= max_cnt())) {
    lRes += est_sum_1_eq(aSumB);
  }
  // if((min_cnt() <= 2) && (2 <= max_cnt())) {
  //   lRes += est_sum_2_eq(aSumB);
  // }
  for(uint k = std::max<uint>(2, min_cnt()); k <= max_cnt(); ++k) {
    lRes += est_sum_k_eq(k, aSumB);
  }
  return lRes;
}

/*
 *  still sum(B) = b
 *  but using integer compositions either limited or unlimited
 */

double
SimpleProfileGbH::estimate_sum_eq_ic(const uint aSumB) const {
  const double p = get_p();
  const double F = freq_tot_k_uda();
  double lRes = 0;
  const uint l = min_B();
  const uint u = max_B();
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
     lRes += std::pow(p, k) * f_int_comp_csd_rec(k, l, u, aSumB) * F;
  }
  return lRes;
}

double
SimpleProfileGbH::estimate_sum_eq_ic_lim(const uint aSumB, const uint aLim) const {
  const double p = get_p();
  const double F = freq_tot_k_uda(); 
  double lRes = 0;
  const uint l = min_B();
  const uint u = max_B();
  assert(1 <= aLim);
  assert(3 >= aLim);
  for(uint k = min_cnt(); (k <= max_cnt()) && (k <= aLim); ++k) {
    lRes += std::pow(p, k) * f_int_comp_csd_rec(k, l, u, aSumB) * F;
  }
  for(uint k = aLim + 1; k <= max_cnt(); ++k) {
    lRes += est_sum_k_eq(k, aSumB);
  }
  return lRes;
}

/*
 *  sum(B) between aLb and aUb
 */

double
SimpleProfileGbH::estimate_sum_between(const double aLb, const double aUb) const {
  double lRes = 0;
  if((min_cnt() <= 1) && (1 <= max_cnt())) {
    lRes += est_sum_1_between(aLb, aUb);
  }
  if((min_cnt() <= 2) && (2 <= max_cnt())) {
    lRes += est_sum_2_between(aLb, aUb);
  }
  for(uint k = std::max<uint>(3, min_cnt()); k <= max_cnt(); ++k) {
    lRes += est_sum_k_between(k, aLb, aUb);
  }
  return lRes;
}


// if we add a restriction such that the count is in [lb,ub]
// that is, we have a having clause of the form
// having sum(A) in [lb,ub] and count(*) in [aLbCnt,aUbCnt]

double
SimpleProfileGbH::estimate_sum_eq_cnt(const double aSumB, const uint aLbCnt, const uint aUbCnt) const {
  double lRes = 0;
  const uint lCntLb = std::max<uint>(min_cnt(), aLbCnt);
  const uint lCntUb = std::min<uint>(max_cnt(), aUbCnt);
  if((lCntLb <= 1) && (1 <= lCntUb)) {
    lRes += est_sum_1_eq(aSumB);
  }
  if((lCntLb <= 2) && (2 <= lCntUb)) {
    lRes += est_sum_2_eq(aSumB);
  }
  for(uint k = std::max<uint>(3,lCntLb); k <= lCntUb; ++k) {
    lRes += est_sum_k_eq(k, aSumB);
  }
  return lRes;

}




/*
 *  sum between
 */


double
SimpleProfileGbH::estimate_sum_between_and_cnt_between(const double aLbSum, const double aUbSum, 
                                                       const uint   aLbCnt, const uint   aUbCnt) const {
  double lRes = 0;
  const uint lCntLb = std::max<uint>(min_cnt(), aLbCnt);
  const uint lCntUb = std::min<uint>(max_cnt(), aUbCnt);
  if((lCntLb <= 1) && (1 <= lCntUb)) {
    lRes += est_sum_1_between(aLbSum, aUbSum);
  }
  if((lCntLb <= 2) && (2 <= lCntUb)) {
    lRes += est_sum_2_between(aLbSum, aUbSum);
  }
  for(uint k = std::max<uint>(3,lCntLb); k <= lCntUb; ++k) {
    lRes += est_sum_k_between(k, aLbSum, aUbSum);
  }
  return lRes;
}

// TODO
double
SimpleProfileGbH::estimate_sum_between_or_cnt_between(const double aLbSum, const double aUbSum,
                                                      const uint   aLbCnt, const uint   aUbCnt) const {
  assert(0 < aLbCnt);
  const uint a = min_cnt();
  const uint b = max_cnt();
  // [a, b] setminus [aLbCnt, aUbCnt]
  // = [a,aLbCnt - 1], [aUbCnt + 1, b]
  const uint lLb1 = a;
  const uint lUb1 = aLbCnt - 1;
  const uint lLb2 = aUbCnt + 1;
  const uint lUb2 = b;
  double lRes = 0;
  lRes = estimate_cnt_between(aLbCnt, aUbCnt);
  if(lLb1 <= lUb1) {
    lRes += estimate_sum_between_and_cnt_between(aLbSum, aUbSum, lLb1, lUb1);
  }
  if(lLb2 <= lUb2) {
    lRes += estimate_sum_between_and_cnt_between(aLbSum, aUbSum, lLb2, lUb2);
  }

  return lRes;
}   




// if a selectivity predicate is applied before the grouping
// and this selectivity predicate has the selectivity aSelectivity
double
SimpleProfileGbH::estimate_sum_eq_sel(const double aSumB, const double aSelectivity) const {
  double lRes = 0;
#ifdef USE_SEL_A
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
    // std::cout << "     lRes(1) = " << lRes << std::endl;
    for(uint j = 1; j <= k; ++j) {
      // std::cout << "     lRes(2) = " << lRes << std::endl;
      lRes +=   fBinomGamma(k, j)
              * std::pow(1 - aSelectivity, k - j)
              * std::pow(    aSelectivity, j)
              * est_sum_eq(j, aSumB);
    }
  }
#else
  double lSumB = aSumB * (1 / aSelectivity);
  if(is_ui_dense_B()) {
     lSumB = (uint) std::max<double>(1, std::round((double) aSumB * (1/aSelectivity)));
  }
  std::cout << "estimate_sum_eq_sel: lSumB = " << lSumB << std::endl;
  if((lSumB < min_B()) || (max_B() < lSumB)) {
    return 0;
  }
  #ifdef USE_SEL_YAO_APPROX
    double lProd = (1 - aSelectivity); // from 1 - (1 - s)^k this contains (1 - s)^k
    if((min_cnt() <= 1) && (1 <= max_cnt())) {
      lRes += (1 - lProd) * est_sum_1_eq(lSumB);
    }
    lProd *= (1 - aSelectivity);
    if((min_cnt() <= 2) && (2 <= max_cnt())) {
      lRes += (1 - lProd) * est_sum_2_eq(lSumB);
    }
    for(uint k = std::max<uint>(3, min_cnt()); k <= max_cnt(); ++k) {
      lProd *= (1 - aSelectivity);
      lRes  += (1 - lProd) * est_sum_k_eq(k, lSumB);
    }
  #else
    const double lProd = (1 - aSelectivity); // only for case 1
    if((min_cnt() <= 1) && (1 <= max_cnt())) {
      lRes += (1 - lProd) * est_sum_1_eq(lSumB);
    }
    if((min_cnt() <= 2) && (2 <= max_cnt())) {
      lRes += sp(2, aSelectivity) * est_sum_2_eq(lSumB);
    }
    for(uint k = std::max<uint>(3, min_cnt()); k <= max_cnt(); ++k) {
      lRes  += sp(k, aSelectivity) * est_sum_k_eq(k, lSumB);
    }
  #endif
#endif
    return std::max<double>(1, lRes);
}

double
SimpleProfileGbH::estimate_sum_between_sel(const double aLb, const double aUb, const double aSelectivity) const {
  if(aLb > aUb) {
    return 0;
  }
  double lRes = 0;
  #ifdef USE_SEL_A
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
    // std::cout << "     lRes(1) = " << lRes << std::endl;
    for(uint j = 1; j <= k; ++j) {
      // std::cout << "     lRes(2) = " << lRes << std::endl;
      lRes +=   fBinomGamma(k, j)
              * std::pow(1 - aSelectivity, k - j)
              * std::pow(    aSelectivity, j)
              * est_sum_between(j, aLb, aUb);
    }
  }
#else
  double lLb = aLb * (1 / aSelectivity);
  double lUb = aUb * (1 / aSelectivity);
  if(is_ui_dense_B()) {
    lLb = (uint) std::max<double>(1, std::round((double) aLb * (1/aSelectivity)));
    lUb = (uint) std::max<double>(1, std::round((double) aUb * (1/aSelectivity)));
  }
  #ifdef USE_SEL_YAO_APPROX
    double lProd = (1 - aSelectivity); // from 1 - (1 - s)^c this contains (1 - s)^c
    if((min_cnt() <= 1) && (1 <= max_cnt())) {
      lRes += (1 - lProd) * est_sum_1_between(lLb, lUb);
    }
    lProd *= (1 - aSelectivity);
    if((min_cnt() <= 2) && (2 <= max_cnt())) {
      lRes += (1 - lProd) * est_sum_2_between(lLb, lUb);
    }
    for(uint k = std::max<uint>(3, min_cnt()); k <= max_cnt(); ++k) {
      lProd *= (1 - aSelectivity);
      lRes += (1 - lProd) * est_sum_k_between(k, lLb, lUb);
    }
  #else
    const double lProd = (1 - aSelectivity);
    if((min_cnt() <= 1) && (1 <= max_cnt())) {
      lRes += (1 - lProd) * est_sum_1_between(lLb, lUb);
    }
    if((min_cnt() <= 2) && (2 <= max_cnt())) {
      lRes += sp(2, aSelectivity) * est_sum_2_between(lLb, lUb);
    }
    for(uint k = std::max<uint>(3, min_cnt()); k <= max_cnt(); ++k) {
      lRes += sp(k, aSelectivity) * est_sum_k_between(k, lLb, lUb);
    }
  #endif
#endif
  return lRes;
}

/*
 *  avg
 */

double 
SimpleProfileGbH::estimate_avg_between(const double aLbAvgB, const double aUbAvgB) const {
  double lRes = 0;
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
    if(1 == k) {
      lRes += est_sum_1_between(aLbAvgB, aUbAvgB);
    }
    if(2 == k) {
      lRes += est_sum_2_between(2 * aLbAvgB, 2 * aUbAvgB);
    } else {
      lRes += est_sum_k_between(k, k * aLbAvgB, k * aUbAvgB);
    }
  }
  return lRes;
}

// TODO
double 
SimpleProfileGbH::estimate_avg_between_and_cnt_between(const double aLbAvgB, const double aUbAvgB,
                                                       const uint   aLbCnt,  const uint   aUbCnt) const {
  const uint lLbCnt = std::max<uint>(aLbCnt, min_cnt());
  const uint lUbCnt = std::min<uint>(aUbCnt, max_cnt());
  double lRes = 0;
  for(uint k = lLbCnt; k <= lUbCnt; ++k) {
    lRes += estimate_sum_between_and_cnt_between(k * aLbAvgB, k * aUbAvgB, k, k);
  }
  return lRes; 
}

// TODO
double
SimpleProfileGbH::estimate_avg_between_or_cnt_between(const double aLbAvgB, const double aUbAvgB,
                                                      const uint   aLbCnt,  const uint   aUbCnt) const {
  // const uint lLbCnt = std::max<uint>(aLbCnt, min_cnt());
  // const uint lUbCnt = std::min<uint>(aUbCnt, max_cnt());
  double lRes = 0;
  // for(uint k = lLbCnt; k <= lUbCnt; ++k) {
  //   lRes += estimate_sum_between_and_cnt_between(k * aLbAvgB, k * aUbAvgB, k, k);
  // }
  return lRes;
}



// TODO
double
SimpleProfileGbH::estimate_avg_between_sel(const double aLb, const double aUb,
                                           const double aSelectivity) const {
  double lRes = 0;
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
    for(uint j = 1; j <= k; ++j) {
      lRes +=   fBinomGamma(k, j)
              * std::pow(1 - aSelectivity, k - j)
              * std::pow(    aSelectivity, j)
              * est_sum_between(j, j * aLb, j * aUb);
    }
  }
  return lRes;
}

/*
 *  min
 */


double 
SimpleProfileGbH::estimate_min_eq(const double aVal) const {
  if((aVal < min_B()) || (max_B() < aVal)) {
    return 0;
  }
  #ifdef ALT_MINMAX_A
  double lRes = 0;
  const double b =  (aVal - min_B()) / (max_B() - min_B()); // prob one in [min_b,aVal)
    for(uint k = min_cnt(); k <= max_cnt(); ++k) {
       if(1 == k) {
         lRes += get_p() * freq_tot_k(1);
       } else {
         lRes += k * get_p() * std::pow(1 - b, k - 1) * freq_tot_k(k); // clearly too high
       }
    }
  #endif
  #ifdef ALT_MINMAX_B
  double lRes = 0;
  double lFak = 0;
  if(aVal < (max_B() - avg_dist_B())) {
    for(uint k = min_cnt(); k <= max_cnt(); ++k) {
      if(1 == k) {
        lRes += get_p() * freq_tot_k(1);
        continue;
      }
      if(max_B() < (aVal + k * avg_dist_B())) { break; }
        lFak =   mt::fLogBinomGamma<double>(((max_B() - aVal + avg_dist_B()) / avg_dist_B()), k - 1)
               - mt::fLogBinomGamma<double>(nodv_B(), k);
      // std::cout << "k = " << k 
      //           << ", fak = " << std::exp(lFak) 
      //           << ", inc = " << std::exp(lFak) * freq_tot_k(k)
      //           << std::endl;
        lRes += std::exp(lFak) * freq_tot_k(k);
      }
  } else {
    for(uint k = min_cnt(); k <= max_cnt(); ++k) {
      lRes += std::pow(get_p(), k) * freq_tot_k(k);
    }
  }
  #endif
  return lRes;
}

double 
SimpleProfileGbH::estimate_min_eq_cnt(const double aVal, const uint aLbCnt, const uint aUbCnt) const {
  if((aVal < min_B()) || (max_B() < aVal)) {
    return 0;
  }
  double lRes = 0;
  uint lLbCnt = aLbCnt;
  uint lUbCnt = aUbCnt;
  cut_cnt(lLbCnt, lUbCnt);
  // std::cout << "estimate_min_eq_cnt k in [" << aLbCnt << ',' << aUbCnt << "]"
  //           << ", [" << lLbCnt << ',' << lUbCnt << "]:" << std::endl;
  if(lUbCnt < lLbCnt) {
    return 0;
  }
  const double b =  (aVal - min_B()) / (max_B() - min_B()); // prob one in [min_b,aVal)
  for(uint k = lLbCnt; k <= lUbCnt; ++k) {
     if(1 == k) {
       lRes += get_p() * freq_tot_k(1);
     } else {
       lRes += k * get_p() * std::pow(1 - b, k - 1) * freq_tot_k(k);
     }
  }
  return lRes;
}

double 
SimpleProfileGbH::estimate_min_eq_sel(const double aVal, const double aSelectivity) const {
  double lRes = 0;
  const double b =  (aVal - min_B()) / (max_B() - min_B()); // prob one in [min_b,aVal)
  double lEstMinEq = 0;
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
    // std::cout << "     lRes(1) = " << lRes << std::endl;
    for(uint j = 1; j <= k; ++j) {
      // std::cout << "     lRes(2) = " << lRes << std::endl;
      if(1 == j) {
        lEstMinEq = get_p() * freq_tot_k(1);
      } else {
        lEstMinEq = j * get_p() * std::pow(1 - b, j - 1) * freq_tot_k(j);
      }
      lRes +=   fBinomGamma(k, j)
              * std::pow(1 - aSelectivity, k - j)
              * std::pow(    aSelectivity, j)
              * lEstMinEq;
    }
  }
  return lRes;
}

double 
SimpleProfileGbH::estimate_min_between(const double aLb, const double aUb) const {
  const double lLb = std::max(aLb, min_B());
  const double lUb = std::min(aUb, max_B());
  if(lUb < lLb) {
    return 0;
  }
  if(lLb == lUb) {
    return estimate_min_eq(lLb);
  }

  double lRes = 0;
  const double a =  (lUb - lLb     + avg_dist_B()) / nodv_B(); // prob one in [lb,ub]
  const double b =  (lLb - min_B() + avg_dist_B()) / nodv_B(); // prob one in [min,lb]
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
    if(1 == k) {
      lRes += a * freq_tot_k(1);
    } else {
      lRes += k * a * std::pow(1-b, k - 1) * freq_tot_k(k);
    }
  }
  return lRes;
}

double 
SimpleProfileGbH::estimate_min_between_and_cnt_between(const double aLbMinB, const double aUbMinB,
                                                       const uint   aLbCnt,  const uint   aUbCnt) const {
  const double lLb = std::max(aLbMinB, min_B());
  const double lUb = std::min(aUbMinB, max_B());
  if(lUb < lLb) {
    return 0;
  }
  if(lLb == lUb) {
    return estimate_min_eq_cnt(lLb, aLbCnt, aUbCnt);
  }

  uint lLbCnt = aLbCnt;
  uint lUbCnt = aUbCnt;
  cut_cnt(lLbCnt, lUbCnt);

  double lRes = 0;
  const double a =  (lUb - lLb     + avg_dist_B()) / nodv_B(); // prob one in [lb,ub]
  const double b =  (lLb - min_B() + avg_dist_B()) / nodv_B(); // prob one in [min,lb]
  for(uint k = lLbCnt; k <= lUbCnt; ++k) {
    if(1 == k) {
      lRes += a * freq_tot_k(1);
    } else {
      lRes += k * a * std::pow(1-b, k - 1) * freq_tot_k(k); 
    }
  } 
  return lRes;
}

double 
SimpleProfileGbH::estimate_min_between_sel(const double aLbMinB, const double aUbMinB, 
                                           const double aSelectivity) const {
  const double lLb = std::max(aLbMinB, min_B());
  const double lUb = std::min(aUbMinB, max_B());
  if(lUb < lLb) {
    return 0;
  }
  if(lLb == lUb) {
    return estimate_min_eq_sel(lLb, aSelectivity);
  }

  double lRes = 0;
  const double a =  ((lUb - lLb     + avg_dist_B()) / avg_dist_B()) / nodv_B(); // prob one in [lb,ub]
  const double b =  ((lLb - min_B() + avg_dist_B()) / avg_dist_B()) / nodv_B(); // prob one in [min,lb]
  double lEstMinRg = 0;
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
    for(uint j = 1; j <= k; ++j) {
      if(1 == j) {
        lEstMinRg = a * freq_tot_k(1);
      } else {
        lEstMinRg = j * a * std::pow(1 - b, j - 1) * freq_tot_k(j);
      }
      lRes +=   fBinomGamma(k, j)
              * std::pow(1 - aSelectivity, k - j)
              * std::pow(    aSelectivity, j)
              * lEstMinRg;
    }
  }
  return lRes;
}


/*
 *  max
 */

double 
SimpleProfileGbH::estimate_max_eq(const double aVal) const {
  if((aVal < min_B()) || (max_B() < aVal)) {
    return 0;
  }
  #ifdef ALT_MINMAX_A
  double lRes = 0;
  const double a =  (aVal - min_B()) / (max_B() - min_B()); // prob one in [min_b, b]
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
     if(1 == k) {
       lRes += get_p() * freq_tot_k(1);
     } else {
       lRes += k * get_p() * std::pow(a, k - 1) * freq_tot_k(k); // clearly too high
     }
  }
  #endif
  #ifdef ALT_MINMAX_B
  double lRes = 0;
  double lFak = 0;
  if(aVal > (min_B() + avg_dist_B())) {
    for(uint k = min_cnt(); k <= max_cnt(); ++k) {
      if(1 == k) {
        lRes += get_p() * freq_tot_k(1);
        continue;
      }
      lFak =   mt::fLogBinomGamma<double>(((aVal - min_B() + avg_dist_B()) / avg_dist_B()), k - 1)
             - mt::fLogBinomGamma<double>(nodv_B(), k);
      // std::cout << "k = " << k 
      //           << ", fak = " << std::exp(lFak) 
      //           << ", inc = " << std::exp(lFak) * freq_tot_k(k)
      //           << std::endl;
        lRes += std::exp(lFak) * freq_tot_k(k);
      }
  } else {
    for(uint k = min_cnt(); k <= max_cnt(); ++k) {
      lRes += std::pow(get_p(), k) * freq_tot_k(k);
    }
  }
  #endif
  return lRes;
}

double 
SimpleProfileGbH::estimate_max_eq_cnt(const double aVal, const uint aLbCnt, const uint aUbCnt) const {
  if((aVal < min_B()) || (max_B() < aVal)) {
    return 0;
  }
  uint lLbCnt = aLbCnt;
  uint lUbCnt = aUbCnt;
  cut_cnt(lLbCnt, lUbCnt);
  if(lUbCnt < lLbCnt) {
    return 0;
  }
  double lRes = 0;
  const double a =  (aVal - min_B()) / (max_B() - min_B()); // prob one in [min_b, b]
  for(uint k = lLbCnt; k <= lUbCnt; ++k) {
     if(1 == k) {
       lRes += get_p() * freq_tot_k(1);
     } else {
       lRes += k * get_p() * std::pow(a, k - 1) * freq_tot_k(k); // clearly too high
     }
  }
  return lRes;
}

double 
SimpleProfileGbH::estimate_max_eq_sel(const double aVal, const double aSelectivity) const {
  double lRes = 0;
  const double a =  (aVal - min_B()) / (max_B() - min_B()); // prob one in [min_b, b]
  double lEstMaxEq = 0;
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
    for(uint j = 1; j <= k; ++j) {
      if(1 == j) {
        lEstMaxEq = get_p() * freq_tot_k(1);
      } else {
        lEstMaxEq = j * get_p() * std::pow(a, j - 1) * freq_tot_k(j);
      }
      lRes +=   fBinomGamma(k, j)
              * std::pow(1 - aSelectivity, k - j)
              * std::pow(    aSelectivity, j)
              * lEstMaxEq;
    }
  }
  return lRes;
}

double 
SimpleProfileGbH::estimate_max_between(const double aLbMaxB, const double aUbMaxB) const {
  const double lLb = std::max(aLbMaxB, min_B());
  const double lUb = std::min(aUbMaxB, max_B());
  if(lUb < lLb) {
    return 0;
  }
  if(lLb == lUb) {
    return estimate_max_eq(lLb);
  }

  double lRes = 0;
  const double a =  ((lUb - lLb     + avg_dist_B()) / avg_dist_B()) / nodv_B(); // prob one in [lb,ub]
  const double b =  ((max_B() - lUb + avg_dist_B()) / avg_dist_B()) / nodv_B(); // prob one in [ub,max]
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
    if(1 == k) {
      lRes += a * freq_tot_k(1);
    } else {
      lRes += k * a * std::pow(1-b, k - 1) * freq_tot_k(k); 
    }
  } 
  return lRes;
}

double 
SimpleProfileGbH::estimate_max_between_and_cnt_between(const double aLbMaxB, const double aUbMaxB,
                                                       const uint   aLbCnt,  const uint   aUbCnt) const {
  const double lLb = std::max(aLbMaxB, min_B());
  const double lUb = std::min(aUbMaxB, max_B());
  if(lUb < lLb) {
    return 0;
  }
  if(lLb == lUb) {
    return estimate_max_eq_cnt(lLb, aLbCnt, aUbCnt);
  }

  uint lLbCnt = aLbCnt;
  uint lUbCnt = aUbCnt;
  cut_cnt(lLbCnt, lUbCnt);

  double lRes = 0;
  const double a =  (lUb - lLb     + avg_dist_B()) / nodv_B(); // prob one in [lb,ub]
  const double b =  (max_B() - lUb + avg_dist_B()) / nodv_B(); // prob one in [ub,max]
  for(uint k = lLbCnt; k <= lUbCnt; ++k) {
    if(1 == k) {
      lRes += a * freq_tot_k(1);
    } else {
      lRes += k * a * std::pow(1-b, k - 1) * freq_tot_k(k);
    }
  }
  return lRes;
}

double 
SimpleProfileGbH::estimate_max_between_sel(const double aLbMaxB, const double aUbMaxB, const double aSelectivity) const {
  const double lLb = std::max(aLbMaxB, min_B());
  const double lUb = std::min(aUbMaxB, max_B());
  if(lUb < lLb) {
    return 0;
  }
  if(lLb == lUb) {
    return estimate_max_eq_sel(lLb, aSelectivity);
  }

  double lRes = 0;  
  double lEstMaxRg = 0;
  const double a =  (lUb - lLb     + avg_dist_B()) / nodv_B(); // prob one in [lb,ub]
  const double b =  (max_B() - lUb + avg_dist_B()) / nodv_B(); // prob one in [ub,max]
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
    for(uint j = 1; j <= k; ++j) {
      if(1 == j) {
        lEstMaxRg = a * freq_tot_k(1);
      } else {
        lEstMaxRg = j * a * std::pow(1 - b, j - 1) * freq_tot_k(j);
      }
      lRes +=   fBinomGamma(k, j)
              * std::pow(1 - aSelectivity, k - j)
              * std::pow(    aSelectivity, j)
              * lEstMaxRg;
    }
  }
  return lRes;

}

// helper routines for sum

double
SimpleProfileGbH::est_sum_eq(const uint aK, const double aVal) const {
  if(1 == aK) { return est_sum_1_eq(aVal); }
  if(2 == aK) { return est_sum_2_eq(aVal); }
  if(3 == aK) { return est_sum_3_eq(aVal); }
  return est_sum_k_eq(aK, aVal);
}

double
SimpleProfileGbH::est_sum_1_eq(const double aVal) const {
  if((min_B_k(1) > aVal) || (max_B_k(1) < aVal)) {
    return 0;
  }
  return (get_p() * freq_tot_k(1));
};

double
SimpleProfileGbH::est_sum_2_eq(const double aVal) const {
  if((min_B_k(2) > aVal) || (max_B_k(2) < aVal)) {
    return 0;
  }
  // std::cout << std::endl;
  // std::cout << "no_comp_2(2, aVal, 50) = " << no_comp_2(2, aVal, 50) << std::endl;
  double lRes = (get_p() * get_p() * no_comp_2(2, aVal, 50) * freq_tot_k(2));
  // std::cout << "est_2_eq(" << aVal << ") = " << lRes << std::endl;
  return lRes;
};

double
SimpleProfileGbH::est_sum_3_eq(const double aVal) const {
  if((min_B_k(2) > aVal) || (max_B_k(2) < aVal)) {
    return 0;
  } 
  // std::cout << std::endl;
  // std::cout << "no_comp_2(3, aVal, 50) = " << no_comp_2(3, aVal, 50) << std::endl;
  double lRes = (get_p() * get_p() * no_comp_2(3, aVal, 50) * freq_tot_k(3));
  // std::cout << "est_3_eq(" << aVal << ") = " << lRes << std::endl;
  return lRes;
};

double
SimpleProfileGbH::est_sum_k_eq(const uint aK, const double aVal) const {
  // constexpr double DELTA = 0.5;
  constexpr double DELTA = 0.4;

  if((min_B_k(aK) > aVal) || (max_B_k(aK) < aVal)) {
    return 0;
  }
  const double v = aVal;
  double lRes = 0;
  const double lCdfDiff = cdf_k(aK, v + DELTA) - cdf_k(aK, v - DELTA);
  lRes = lCdfDiff * freq_tot_k(aK);
  return lRes;
};

double
SimpleProfileGbH::est_sum_between(const uint aK, const double aLb, const double aUb) const {
  if(1 == aK) { return est_sum_1_between(aLb, aUb); }
  if(2 == aK) { return est_sum_2_between(aLb, aUb); }
  return est_sum_k_between(aK, aLb, aUb);
}


double
SimpleProfileGbH::est_sum_1_between(const double aLb, const double aUb) const {
  double a = aLb;
  double b = aUb;
  cut_B_k(1, a, b);
  if(a > b) { return 0; }
  if(a == b) { return est_sum_1_eq(a); }
  double lRes = 0;
  if(is_ui_dense_B()) {
    lRes = ((double) (b - a + 1) / (max_B_k(1) - min_B_k(1))) * freq_tot_k(1);
  } else {
    lRes = ((double) (b - a) / (max_B_k(1) - min_B_k(1))) * freq_tot_k(1);
  }
  return lRes;
};

double
SimpleProfileGbH::est_sum_2_between(const double aLb, const double aUb) const {
  double a = aLb;
  double b = aUb;
  cut_B_k(2, a, b);
  if(a > b) { return 0; }
  if(a == b) { return est_sum_2_eq(a); }
  double lRes = 0;

  if(is_ui_dense_B() && (10 >= (b - a))) {
    for(uint x = a; x <= b; ++x) {
      lRes += est_sum_2_eq(x);
    }
  } else {
    const double lCdfDiff = cdf_k(2, b + 0.4) - cdf_k(2, a - 0.4);
    lRes = (lCdfDiff * freq_tot_k(2));
  }
  return lRes;
};

double
SimpleProfileGbH::est_sum_k_between(const uint aK, const double aLb, const double aUb) const {
  double a = aLb;
  double b = aUb;
  cut_B_k(aK, a, b);
  if(b < a) { return 0; }
  if(a == b) { return est_sum_k_eq(aK, a); }
  const double lLb = a;
  const double lUb = b;
  const double lCdfDiff = cdf_k(aK, lUb + 0.4) - cdf_k(aK, lLb - 0.4);
  return (lCdfDiff * freq_tot_k(aK));
};

void
SimpleProfileGbH::cut_B_k(const uint k, double& a, double& b) const {
  a = std::max<double>(a, min_B_k(k));
  b = std::min<double>(b, max_B_k(k));
};

void
SimpleProfileGbH::cut_cnt(uint& a, uint& b) const {
  a = std::max<uint>(a, min_cnt());
  b = std::min<uint>(b, max_cnt());
}

// survival probability of a group of size k
// if a selection predicate with aSelectivity has been applied
// either by approximation
// or by Yao
double
SimpleProfileGbH::sp(const uint k, const double aSelectivity) const {
  double lRes = 0;
#ifdef USE_SEL_YAO_APPROX
    lRes = (1 - std::pow(1 - aSelectivity, k));
#else
    lRes = fYaoProbGamma(card_R(), k, card_R() * aSelectivity);
#endif
  return lRes;
}


std::ostream& 
SimpleProfileGbH::print(std::ostream& os) const {
  os << "SimpleEstimatorGbH:" << std::endl;
  os << "  |R|        = " << std::setw(8) << card_R()     << std::endl
     << "  nodv(A)    = " << std::setw(8) << nodv_A()     << std::endl
     << "  min_cnt    = " << std::setw(8) << min_cnt()    << std::endl
     << "  max_cnt    = " << std::setw(8) << max_cnt()    << std::endl
     << "  ndv_cnt    = " << std::setw(8) << nodv_cnt()   << std::endl
     << "  min_B      = " << std::setw(8) << min_B()      << std::endl
     << "  max_B      = " << std::setw(8) << max_B()      << std::endl
     << "  ndv_B      = " << std::setw(8) << nodv_B()     << std::endl
     << "  dist_B     = " << std::setw(8) << avg_dist_B() << std::endl
     << "  mean_B     = " << std::setw(8) << mean_B()     << std::endl
     << "  variance_B = " << std::setw(8) << variance_B() << std::endl
     << "  std_dev_B  = " << std::setw(8) << std_dev_B()  << std::endl
     << "  pr(b = c)  = " << std::setw(8) << get_p()      << std::endl
     << "  is_ui_dns  = " << std::setw(8) << (is_ui_dense_B() ? "true" : "false") << std::endl
     << std::endl;
  os << "     k    min_k    max_k   freq_k   mean_k stddev_k    var_k" << std::endl;
  for(uint k = min_cnt(); k <= max_cnt(); ++k) {
    os << "    " << std::setw(2) << k
       << ' '    << std::setw(8) << min_B_k(k)
       << ' '    << std::setw(8) << max_B_k(k)
       << ' '    << std::setw(8) << freq_tot_k(k)
       << ' '    << std::setw(8) << mean_B_k(k)
       << ' '    << std::setw(8) << std_dev_B_k(k)
       << ' '    << std::setw(8) << variance_B_k(k)
       << std::endl;
  }
  os << "----" << std::endl;
  os << std::endl;
  return os;
}


