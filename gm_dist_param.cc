#include "gm_dist_param.hh"

namespace gm {

double
dist_param_t::moment_3() const {
  const double lMom2 = moment_2();
  return ((skew() * std_dev() * std_dev() * std_dev())
           - 2 * mean() * mean() * mean() + 3 * mean() * lMom2);
}

dist_param_t
operator*(const dist_param_t& x,
          const dist_param_t& y) {
   // Mean is is 1st moment, just multiply
   double lMean = x.mean() * y.mean();

   // Derive 2nd moments, multiply and get resulting std_dev
   double lM2x    = x.std_dev() * x.std_dev() + x.mean() * x.mean();
   double lM2y    = y.std_dev() * y.std_dev() + y.mean() * y.mean();
   double lM2res  = lM2x * lM2y;
   double lStdDev = std::sqrt(std::max(lM2res - lMean * lMean, 0.0));

   // Derive 3rd moments, multiply and get resulting skew
   double x3M  =   (x.skew() * x.std_dev() * x.std_dev() * x.std_dev())
               - 2 * x.mean() * x.mean() * x.mean() + 3 * x.mean() * lM2x;
   double y3M  =   (y.skew() * y.std_dev() * y.std_dev() * y.std_dev())
               - 2 * y.mean() * y.mean() * y.mean() + 3 * y.mean() * lM2y;
   double lSkew =   lStdDev > .001 ?
                      (x3M * y3M - 3 * lMean * lM2res + 2 * lMean * lMean * lMean)
                    / (lStdDev * lStdDev * lStdDev)
                 : 0;

   return dist_param_t(lMean, lStdDev, lSkew);

}

bool
dist_param_t::get_param_skew_norm_b(double& aXi,
                                    double& aEta,
                                    double& aLambda) const {
  constexpr double l_1_3 = double{1}/double{3};
  #ifdef MACOS
  const double l_pi = M_PI;
  #else
  constexpr double l_pi = std::numbers::pi_v<double>;
  #endif
  // some other abbreviations
  const double s2 = variance();
  // const double m3 = mean_3();
  // rhs formulas
  const double l_gamma = skew();     
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
  const double l_xi    = mean() - l_omega * m;
  // final results
  aXi     = l_xi;
  aEta    = l_omega;
  aLambda = l_alpha;

  return lRes;
}

} // end namespace

