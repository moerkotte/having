#ifndef VARIANCE_HH
#define VARIANCE_HH

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <cinttypes>
#include <numbers>

/*
 *  class Variance
 *  calculates variance, standard deviation
 *  k-te Momente m_k for k = 1, 2, 3
 *  where m_k = E[X^k]
 *  Nicht berechnet werden k-te absolute Momente M_k
 *  M_k = E[|X|^k]
 *  oder die zentralen Momente
 *  mu_k = E[(X - mu)^k]
 *  oder die zentralen standardisierten Momente
 * mu_{Sk} = E[((X - mu) / s)^k], mit s ist Standardabweichung
 */


template<class Tdom>
class Variance {
  public:
    Variance() : _sum(0), _sum_sq(0), _sum_qu(0), _count(0) {}
    ~Variance() {}
  public:
    inline
    void init() { _sum    = 0;
                  _sum_sq = 0;
                  _sum_qu = 0;
                  _count  = 0;
                };
    inline
    void step(Tdom x) {
      _sum    += x;
      _sum_sq += (x*x);
      _sum_qu += (x*x*x);
      ++_count;
    }
    inline
    void step(const Tdom aVal, const uint64_t aFreq) {
      _sum    += (aVal * aFreq);
      _sum_sq += aFreq * (aVal * aVal);
      _sum_qu += aFreq * (aVal * aVal * aVal);
      _count  += aFreq;
    }
    inline
    void reverseStep(const Tdom x) {
      _sum    -= x;
      _sum_sq -= (x*x);
      _sum_qu -= (x*x*x);
      --_count;
    }

    inline void fin() {}
  public:
     inline uint   count() const { return _count; }
     inline double avg() const { return ((double) sum() / (double) count()); }
     inline double mean() const { return avg(); }
     inline Tdom   sum() const { return _sum; }
     inline Tdom   sum_sq() const { return _sum_sq; }
     inline Tdom   sum_qu() const { return _sum_qu; }
     // k-th moments m_k for k = 1, 2, 3
     inline double mean_1() const { return mean(); }
     inline double mean_2() const { return ((double) sum_sq() / ((double) count())); }
     inline double mean_3() const { return ((double) sum_qu() / ((double) count())); }
     double variance() const;
     double sse() const;
     double standardDeviation() const;
     double skewness() const;
     bool   get_param_skew_norm_a(double& aXi, double& aEta, double& aLambda) const;
     bool   get_param_skew_norm_b(double& aXi, double& aEta, double& aLambda) const;
     bool   get_param_beta(double& aAlpha, double& aBeta, const double aLb, const double aUb) const;
  public:
     // normalized to range [0,1] and denormalize to [lb,ub]
     inline double mean_01_1(const double aMin, const double aMax) const {
                     return ((mean_1() - aMin) / (aMax - aMin)); 
                   }
     inline double var_01_1(const double aMin, const double aMax) const {
                     const double d = aMax - aMin;
                     const double f = double{1} / d;
                     return (f * f * variance()); 
                   }
  public:
    void set(const Tdom aSum, const Tdom aSumSq, const Tdom aSumQu, const uint64_t aCnt);
  public:
    std::ostream& print(std::ostream& os, int aKindSet, int aFieldWidth = 8) const {
      os << std::setw(aFieldWidth) << count() << ' ' 
         << std::setw(aFieldWidth) << sum() << ' ' 
         << std::setw(aFieldWidth) << sum_sq() << ' ' 
         << std::setw(aFieldWidth) << variance() << ' '
         << std::endl;
      return os;
    }
  private:
    Tdom     _sum;
    Tdom     _sum_sq; // sum x^2
    Tdom     _sum_qu; // sum x^3
    uint64_t _count;
  private:
    static int _kindset;
  public:
    static int  setToPrint() { return _kindset; }
    static void setToPrint(int aKindset) { _kindset = aKindset; }
};

template<class Tdom>
double
Variance<Tdom>::variance() const {
  const double n_inv = (double{1.0} / (double) count());
  return ( ( n_inv * (double) sum_sq()) - ((n_inv * sum()) * (n_inv * sum())) );
}

template<class Tdom>
double
Variance<Tdom>::sse() const {
  return ( ((double) sum_sq()) - (sum() * sum()) );
}

template<class Tdom>
double
Variance<Tdom>::standardDeviation() const {
  return sqrt(variance());
}

template<class Tdom>
double
Variance<Tdom>::skewness() const {
  const double n_inv = double{1} / ((double) count());
  const double m1 = ((double) sum())    * n_inv;
  const double m2 = ((double) sum_sq()) * n_inv;
  const double m3 = ((double) sum_qu()) * n_inv;
  const double s = standardDeviation();
  return ((m3 - 3 * m1 * m2 + 2 * (m1 * m1 * m1)) / (s * s * s));
}

// params:
// zeta:   position
// eta:    scale
// lambda: skewness
// returns false if delta needed to be cutted
// calculation after Pewsey
// tuts nicht
template<class Tdom>
bool
Variance<Tdom>::get_param_skew_norm_a(double& aXi, double& aEta, double& aLambda) const {
  // constexpr double lPi = M_PI;
  constexpr double l_1_3 = double{1}/double{3};
  constexpr double lPi = std::numbers::pi_v<double>;
  #ifdef MACOS
  const double b   = std::sqrt(2 / lPi);
  const double c   = std::pow(2 / (4 - lPi), l_1_3);
  #else
  constexpr double b   = std::sqrt(2 / lPi);
  constexpr double c   = std::pow(2 / (4 - lPi), l_1_3);
  #endif
  // some other abbreviations
  // const double n = count();
  const double s = standardDeviation();
  const double m3 = mean_3();
  // estimated values for standardized sample (nach Pewsey)
  const double l_xi_S  = -c * std::pow(m3, l_1_3) / s;
  const double l_eta_S = std::sqrt(1 + l_xi_S * l_xi_S);
        double l_delta = - l_xi_S / (b * l_eta_S);
  // check/cut l_delta
  bool lRes = true;
  if(l_delta >= 0.995) {
    l_delta = 0.995;
    lRes = false;
  } else
  if(l_delta <= -0.995) {
    l_delta = -0.995;
    lRes = false;
  }
  // final results
  aXi     = avg() + s * l_xi_S;
  aEta    = s * l_eta_S;
  aLambda = l_delta / std::sqrt(1 - l_delta * l_delta);

  return lRes;
}

// formeln aus Diss Fent
template<class Tdom>
bool
Variance<Tdom>::get_param_skew_norm_b(double& aXi, double& aEta, double& aLambda) const {
  constexpr double l_1_3 = double{1}/double{3};
  constexpr double l_pi = std::numbers::pi_v<double>;
  // some other abbreviations
  const double s2 = variance();
  // const double m3 = mean_3();
  // rhs formulas
  const double l_gamma = skewness();
  const double n = std::pow( (2 * std::abs(l_gamma)) / (4 - l_pi) , l_1_3);
  const double m = n / std::sqrt(1 + n * n);
        double l_delta = std::sqrt(l_pi / 2) * m;
  // check/cut l_delta
  bool lRes = true;
  if(l_delta >= 0.995) {
    l_delta = 0.995;
    lRes = false;
  } else
  if(l_delta <= -0.995) {
    l_delta = -0.995;
    lRes = false;
  }
  const double l_alpha = l_delta / std::sqrt(1 - l_delta * l_delta);
  const double l_omega = std::sqrt(s2 / (1 - m * m));
  const double l_xi    = avg() - l_omega * m;
  // final results
  aXi     = l_xi;
  aEta    = l_omega;
  aLambda = l_alpha;

  return lRes;
}

template<class Tdom>
bool
Variance<Tdom>::get_param_beta(double& aAlpha, 
                               double& aBeta, 
                               const double aLb, 
                               const double aUb) const {
  const double a = aLb; // lower bound of values
  const double c = aUb; // upper bound of values
  const double mu  = mean();
  const double mu2 = mu * mu;
  const double s2  = variance();
  const double lZaehler = (a * c - a * mu - c * mu + mu2 + s2);
  const double lNenner  = s2 * (c - a);
  const double lBruch   = lZaehler/ lNenner;
  aAlpha = (a - mu) * lBruch;
  aBeta  = (mu - c) * lBruch;  // original: ((c - mu) * (...))
  // [lb,ub] = [0,1]:
  // alpha = mean * (lZaehler / s2)
  // beta  = (1 - mean) * (lZaehler / s2)
  return true;
}


template<class Tdom>
void 
Variance<Tdom>::set(const Tdom aSum, const Tdom aSumSq, const Tdom aSumQu, const uint64_t aCnt) {
  _sum = aSum;
  _sum_sq = aSumSq;
  _sum_qu = aSumQu;
  _count  = aCnt;
}

template<class Tdom>
std::ostream& 
operator<<(std::ostream& os, const Variance<Tdom>& aAggr) {
  if(0 == Variance<Tdom>::setToPrint()) {
    os << "min: " << aAggr.min() << ", ";
    os << "max: " << aAggr.max() << ", ";
    os << "avg: " << aAggr.avg() << ", ";
    os << "qmid: " << aAggr.qmiddle();
  } else {
    aAggr.print(os, Variance<Tdom>::setToPrint());
  }
  return os;
}


#endif

