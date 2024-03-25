#ifndef INFRA_FLOAT_OPS_HH
#define INFRA_FLOAT_OPS_HH

// #include <ostream>
// #include <limits>
// #include <inttypes.h>
// #include <cmath>
// 

template<typename Tfloat>
inline
bool
is_equal(const Tfloat x, const Tfloat y, const Tfloat eps) {
  return ((x <= (y + eps)) && (y <= (x + eps)));
}

/*
struct print_double_t {
  double _d;
  print_double_t(const double x) : _d(x) {}
};

inline std::ostream&
operator<<(std::ostream& os, const print_double_t& x) {
  constexpr uint64_t lMax = std::numeric_limits<uint64_t>::max();
  if((lMax > ((uint64_t) x._d)) && (std::floor(x._d) == x._d)) {
    os << ((uint64_t) x._d);
  } else {
    os << x._d;
  }
  return os;
}
*/


#endif
