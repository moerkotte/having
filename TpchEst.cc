#include "TpchEst.hh"
#include "SimpleProfileGbH.hh"


/*
 % result of query Q_c
 \begin{tabular}{cr}
c & freq ($F_{c}$)   \\\hline
1 & 214172  \\
2 & 214434  \\
3 & 214379  \\
4 & 213728  \\
5 & 214217  \\
6 & 214449  \\
7 & 214621
\end{tabular}
*/

TpchEst::TpchEst() : _card_R(6'001'215), _freq_tot(8), _freq_sum(8), _prob_sum(8), _mean(8), _std_dev(8) {
  const bool lTrace = false;
  for(uint k = 0; k < 8; ++k) {
    _freq_sum[k].resize(351);
    _prob_sum[k].resize(351);
  }

  // read result of query Q_d (see ~papers/CardEstHaving/main.pdf)
  // select no_line as a, sum_quant as b, count(*) as c
  // from (select l_orderkey, count(*) as no_line, sum(l_quantity) as sum_quant
  //       from   lineitem
  //       group by l_orderkey)
  // group by no_line, sum_quant
  // order by no_line, sum_quant;

  std::string lFilename("data/nosumcount.dat");
  std::ifstream lIs(lFilename);
  if(!lIs) {
    std::cout << "can't open '" << lFilename << "'." << std::endl;
    assert(0 == 1);
    return;
  }
  uint64_t a = 0, b = 0, c = 0;
  while(!lIs.eof()) {
     lIs >> a >> b >> c;
     if(lIs.eof()) { break; }
     _freq_sum[a][b] = c;
     _freq_tot[a] += c;
  }
  if(lTrace) {
    std::cout << "total frequency for summing up $k$ elements:" << std::endl;
    for(uint k = 1; k <= 7; ++k) {
      std::cout << std::setw(3) << k
                << "  " << std::setw(8) << _freq_tot[k]
                << std::endl;
    }
    std::cout << std::endl;
  }
  // calculate probabilities by dividing by total frequency
  for(uint k = 1; k <= 7; ++k) {
    for(uint c = 0; c <= 350; ++c) {
      _prob_sum[k][c] = _freq_sum[k][c] / _freq_tot[k];
    }
  }
  // calculate mean and standard deviation (goal: approximate by normal distribution)
  if(lTrace) {
    std::cout << "Mean, variance, and standard deviation:" << std::endl;
  }
  for(uint k = 1; k <= 7; ++k) {
    Variance<double> lVar;
    lVar.init();
    for(uint c = 1; c < 350; ++c) {
      lVar.step(c, _freq_sum[k][c]);
    }
    lVar.fin();
    if(lTrace) {
      std::cout << std::setw(3) << k
                << ' ' << std::setw(8) << lVar.mean()
                << ' ' << std::setw(8) << lVar.variance()
                << ' ' << std::setw(8) << lVar.standardDeviation()
                << std::endl;
    }
    _mean[k]    = lVar.mean();
    _std_dev[k] = lVar.standardDeviation();
  }
  if(lTrace) {
    std::cout << std::endl;
  }
}

TpchEst::~TpchEst() {
}

#ifdef USE_BOOST
double
TpchEst::pdf_k(const uint k, const double a) const {
  dist_t lDist(mean(k), std_dev(k));
   return pdf(lDist, a); 
}

double
TpchEst::cdf_k(const uint k, const double a) const {
  dist_t lDist(mean(k), std_dev(k));
  return cdf(lDist, a); 
}
#else
double
TpchEst::pdf_k(const uint k, const double a) const {
   return normal_pdf(a, mean(k), std_dev(k));
}

double
TpchEst::cdf_k(const uint k, const double a) const {
   return normal_cdf(a, mean(k), std_dev(k));
}
#endif

std::ostream&
TpchEst::print_a(std::ostream& os) const {
  for(uint c = 0; c <= 350; ++c) {
    os << std::setw(3) << c;
    for(uint k = 1; k <= 7; ++k) {
      os << ' ' << std::setw(12) << _prob_sum[k][c];
    }
    for(uint k = 1; k <= 7; ++k) {
      os << ' ' << pdf_k(k, c);
    }
    os << std::endl;
  }
  return os;
}

double
TpchEst::estimate_sum_eq(const uint aSumQuantity) const {
  double lRes = 0;
  lRes += est_1_eq(aSumQuantity);
  lRes += est_2_eq(aSumQuantity);
  for(uint k = 3; k <= 7; ++k) {
    lRes += est_k_eq(k, aSumQuantity);
  }
  return lRes;
}

double
TpchEst::estimate_sum_between(const uint aSumQuantityLb, const uint aSumQuantityUb) const {
  double lRes = 0;
  lRes += est_1_between(aSumQuantityLb, aSumQuantityUb);
  lRes += est_2_between(aSumQuantityLb, aSumQuantityUb);
  for(uint k = 3; k <= 7; ++k) {
    lRes += est_k_between(k, aSumQuantityLb, aSumQuantityUb);
  }
  return lRes;
}

// if we add a restriction such that the count is in [lb,ub]
// that is, we have a having clause of the form
// having sum(A) in [lb,ub] and count(*) in [aCntLb,aCntUb]

double
TpchEst::estimate_sum_eq_cnt(const uint aSumQuantity, const uint aCntLb, const uint aCntUb) const {
  double lRes = 0;
  if((aCntLb <= 1) && (1 <= aCntUb)) {
    lRes += est_1_eq(aSumQuantity);
  }
  if((aCntLb <= 2) && (2 <= aCntUb)) {
    lRes += est_2_eq(aSumQuantity);
  }
  for(uint k = std::max<uint>(3,aCntLb); k <= std::min<uint>(7,aCntUb); ++k) {
    lRes += est_k_eq(k, aSumQuantity);
  }
  return lRes;

}

double
TpchEst::estimate_sum_between_and_cnt_between(const uint aLb, const uint aUb, const uint aCntLb, const uint aCntUb) const {
  double lRes = 0;
  if((aCntLb <= 1) && (1 <= aCntUb)) {
    lRes += est_1_between(aLb, aUb);
  }
  if((aCntLb <= 2) && (2 <= aCntUb)) {
    lRes += est_2_between(aLb, aUb);
  }
  for(uint k = std::max<uint>(3,aCntLb); k <= std::min<uint>(7,aCntUb); ++k) {
    lRes += est_k_between(k, aLb, aUb);
  }
  return lRes;
}

// if a selectivity predicate is applied before the grouping
// and this selectivity predicate has the selectivity aSelectivity
double
TpchEst::estimate_sum_eq_sel(const uint aSumQuantity, const double aSelectivity) const {
  double lRes = 0;
  uint lSumQuantity = (uint) std::max<double>(1, std::round((double) aSumQuantity * (1/aSelectivity)));
  double lProd = (1 - aSelectivity); // from 1 - (1 - s)^c this contains (1 - s)^c
  lRes += (1 - lProd) * est_1_eq(lSumQuantity);
  lProd *= (1 - aSelectivity);
  lRes += (1 - lProd) * est_2_eq(lSumQuantity);
  for(uint k = 3; k <= 7; ++k) {
    lProd *= (1 - aSelectivity);
    lRes += (1 - lProd) * est_k_eq(k, lSumQuantity);
  }
  return lRes;
}

double
TpchEst::estimate_sum_between_sel(const uint aLb, const uint aUb, const double aSelectivity) const {
  double lRes = 0;
  uint lLb = (uint) std::max<double>(1, std::round((double) aLb * (1/aSelectivity)));
  uint lUb = (uint) std::max<double>(1, std::round((double) aUb * (1/aSelectivity)));
  double lProd = (1 - aSelectivity); // from 1 - (1 - s)^c this contains (1 - s)^c
  lRes += (1 - lProd) * est_1_between(lLb, lUb);
  lProd *= (1 - aSelectivity);
  lRes += (1 - lProd) * est_2_between(lLb, lUb);
  for(uint k = 3; k <= 7; ++k) {
    lProd *= (1 - aSelectivity);
    lRes += (1 - lProd) * est_k_between(k, lLb, lUb);
  }
  return lRes;
}


double 
TpchEst::est_1_eq(const uint aVal) const {
  if((min_B_k(1) > aVal) || (max_B_k(1) < aVal)) {
    return 0;
  }
  return (get_p() * freq_tot(1));
};

double 
TpchEst::est_2_eq(const uint aVal) const {
  if((min_B_k(2) > aVal) || (max_B_k(2) < aVal)) {
    return 0;
  }
  // std::cout << std::endl;
  // std::cout << "no_comp_2(2, aVal, 50) = " << no_comp_2(2, aVal, 50) << std::endl;
  double lRes = (get_p() * get_p() * no_comp_2(2, aVal, 50) * freq_tot(2));
  // std::cout << "est_2_eq(" << aVal << ") = " << lRes << std::endl;
  return lRes;
};

double 
TpchEst::est_k_eq(const uint aK, const uint aVal) const {
  if((min_B_k(aK) > aVal) || (max_B_k(aK) < aVal)) {
    return 0;
  }
  const double v = aVal;
  double lRes = 0;
  const double lCdfDiff = cdf_k(aK, v + 0.5) - cdf_k(aK, v - 0.5);
  lRes = lCdfDiff * freq_tot(aK);
  return lRes;
};

double 
TpchEst::est_1_between(const uint aLb, const uint aUb) const {
  uint a = aLb;
  uint b = aUb;
  cut_k(1, a, b);
  if(a > b) { return 0; }
  if(a == b) { return est_1_eq(a); }
  double lRes = ((double) (b - a + 1) / ((double) 50)) * freq_tot(1);
  return lRes;
};

double 
TpchEst::est_2_between(const uint aLb, const uint aUb) const {
  uint a = aLb;
  uint b = aUb;
  cut_k(2, a, b);
  double lRes = 0;
  for(uint x = a; x <= b; ++x) {
    lRes += est_2_eq(x);
  }
  return lRes;
};

double 
TpchEst::est_k_between(const uint aK, const uint aLb, const uint aUb) const {
  uint a = aLb;
  uint b = aUb;
  cut_k(aK, a, b);
  if(b < a) { return 0; }
  if(a == b) { return est_k_eq(aK, a); }
  const double lLb = a;
  const double lUb = b;
  const double lCdfDiff = cdf_k(aK, lUb + 0.5) - cdf_k(aK, lLb - 0.5);
  return (lCdfDiff * freq_tot(aK));
};

void   
TpchEst::cut_k(const uint k, uint& a, uint& b) const {
  a = std::max<uint>(a, min_B_k(k));
  b = std::min<uint>(b, max_B_k(k));
};

uint 
TpchEst::min_B_k(const uint aK) const {
  return aK;
};

uint 
TpchEst::max_B_k(const uint aK) const {
  return (aK * 51);
};

uint
TpchEst::get_true_card_cnt_eq(const uint aCnt) const {
  if((0 == aCnt) || (7 < aCnt)) {
    return 0;
  }
  return freq_tot(aCnt);
}

uint
TpchEst::get_true_card_cnt_between(const uint aLbCnt, const uint aUbCnt) const {
  const uint lLb = std::max<uint>(1, aLbCnt);
  const uint lUb = std::min<uint>(7, aUbCnt);
  uint lRes = 0;
  for(uint i = lLb; i <= lUb; ++i) {
    lRes += freq_tot(i);
  }
  return lRes;
}

uint
TpchEst::get_true_card_sum_eq(const uint aVal) const {
  uint lRes = 0;
  for(uint k = 1; k <= 7; ++k) {
    lRes += freq_sum(k, aVal);
  }
  return lRes;
}

uint
TpchEst::get_true_card_sum_eq_and_cnt_between(const uint aVal, const uint aCntLb, const uint aCntUb) const {
  uint lRes = 0;
  for(uint k = std::max<uint>(1, aCntLb); k <= std::min<uint>(7, aCntUb); ++k) {
    lRes += freq_sum(k, aVal);
  }
  return lRes;
}

uint
TpchEst::get_true_card_sum_between(const uint aLb, const uint aUb) const {
  uint lRes = 0;
  for(uint k = 1; k <= 7; ++k) {
    for(uint c = aLb; c <= aUb; ++c) {
      lRes += freq_sum(k, c);
    }
  }
  return lRes;
}

uint
TpchEst::get_true_card_sum_between_and_cnt_between(const uint aLb, const uint aUb, 
                                                   const uint aLbCnt, const uint aUbCnt) const {
  uint lRes = 0;
  for(uint k = std::max<uint>(1, aLbCnt); k <= std::min<uint>(7, aUbCnt); ++k) {
    for(uint s = aLb; s <= aUb; ++s) {
      lRes += freq_sum(k, s);
    }
  }
  return lRes;
}

uint
TpchEst::get_true_card_sum_between_or_cnt_between(const uint aLb, const uint aUb,
                                                  const uint aLbCnt, const uint aUbCnt) const {
  uint lRes = 0;
  for(uint k = 1; k <= 7; ++k) {
    if((aLbCnt <= k) && (k <= aUbCnt)) {
      lRes += freq_tot(k);
    } else {
      for(uint s = aLb; s <= aUb; ++s) {
        lRes += freq_sum(k, s);
      }
    }
  }
  return lRes;
}



uint
TpchEst::get_true_card_avg_between(const uint aLb, const uint aUb) const {
  uint lRes = 0;
  for(uint k = 1; k <= 7; ++k) {
    lRes += get_true_card_avg_between_one_count(aLb, aUb, k);
  }
  return lRes;
}

uint
TpchEst::get_true_card_avg_between_cnt(const uint aLb, const uint aUb, const uint aLbCnt, const uint aUbCnt) const {
  uint lRes = 0;
  for(uint k = std::max<uint>(1, aLbCnt); k <= std::min<uint>(7, aUbCnt); ++k) {
      lRes += get_true_card_avg_between_one_count(aLb, aUb, k);
  }
  return lRes;
}

uint
TpchEst::get_true_card_avg_between_one_count(const uint aLb, const uint aUb, const uint aCnt) const {
  const uint_vt& lFreqSum = _freq_sum[aCnt];
  const double lCnt = aCnt;
  double lAvg = 0;
  uint   lRes = 0;
  for(uint i = 0; i < 350; ++i) {
    if(0 < lFreqSum[i]) {
      lAvg = ((double) i / lCnt);
      if((aLb <= lAvg) && (lAvg <= aUb)) {
        lRes += lFreqSum[i];
      }
    }
  }
  return lRes;
}




