#ifndef OPT_CARD_HAVING_LOOP_HH
#define OPT_CARD_HAVING_LOOP_HH

#include "util.hh"
#include <fstream>

/*
 *  struct loop_spec
 *  contains the specification of 
 *  - a loop 
 *    for(i = _lo_a; i <= _hi_a; i += _step_a)
 *  - and optionally a second loop on some width (for range queries)
 *    for(w = _lo_w; i <= _hi_w; i += _step_w)
 *  the parameter _is_eq distinguishes between these two cases
 *  after it has been filled, the loops can be executed by two routines:
 *  - generate_queries 
 *    (to generate queries with = or between for those cases TpchEst fails to produce the true cardinality)
 *  - check_estimates
 *    (for those cases where TpchEst can produce the true cardinality)
 */

typedef void (*genfun_t)(std::ostream&, queryinstance_t&, const querytemplate_t&);
typedef void (*checkfun_t)(const class SimpleProfileGbH&, 
                           const class TpchEst&,  
                                       queryinstance_t&, 
                           const       querytemplate_t&,
                           const class Cb&);

struct loopspec_t {
  bool _is_eq;
  uint _lo_a;
  uint _hi_a;
  uint _step_a;
  uint _lo_w;
  uint _hi_w;
  uint _step_w;
  void set(const uint aLoA, const uint aHiA, const uint aStepA);
  void set(const uint aLoA, const uint aHiA, const uint aStepA, 
           const uint aLoW, const uint aHiW, const uint aStepW);
  void generate_queries(genfun_t, uint& a, uint& b, std::ostream&, queryinstance_t&, const querytemplate_t&);
  void eval_estimates(checkfun_t, uint& a, uint& b, 
                      const class SimpleProfileGbH&, 
                      const class TpchEst&, 
                                  queryinstance_t&, 
                      const       querytemplate_t&,
                      const class Cb&);
};

#endif

