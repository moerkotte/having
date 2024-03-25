#include "int_comp.hh"

uint64_t
no_comp_1(const uint64_t k, const uint64_t n, const uint64_t m) {
  if((n < k) || (n > k * m)) return 0;
  if(1 == k) { return 1; }
  uint64_t lRes = 0;
  uint64_t lLim = std::min<uint64_t>(m, n - (k - 1));
  // std::cout << "lLim = " << lLim << std::endl;
  for(uint64_t i = 1; i <= lLim; ++i) {
     const uint64_t r = no_comp_1(k - 1, n - i, m);
     lRes += r;
     // std::cout << i << ' ' << r << ' ' << lRes << std::endl;
  }
  return lRes;
}

uint64_t
no_comp_2(const uint k, const uint64_t n, const uint64_t m) {
  return f_int_comp_rec(k, 1, m, n);
}


