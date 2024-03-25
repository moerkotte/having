#include "loop.hh"

void
loopspec_t::set(const uint aLoA, const uint aHiA, const uint aStepA) {
  _is_eq  = true;
  _lo_a   = aLoA;
  _hi_a   = aHiA;
  _step_a = aStepA;
}

void
loopspec_t::set(const uint aLoA, const uint aHiA, const uint aStepA, 
                const uint aLoW, const uint aHiW, const uint aStepW) {
  _is_eq  = false;
  _lo_a   = aLoA;
  _hi_a   = aHiA;
  _step_a = aStepA;
  _lo_w   = aLoW;
  _hi_w   = aHiW;
  _step_w = aStepW;
}

void
loopspec_t::generate_queries(genfun_t aGenFun, uint& a, uint& b, 
                                   std::ostream&    aOsQuery, 
                                   queryinstance_t& aQueryInstance, 
                             const querytemplate_t& aQueryTemplate) {
  if(_is_eq) {
    for(uint lValA = _lo_a; lValA <= _hi_a; lValA += _step_a) {
      a = lValA;
      b = lValA;
      (*aGenFun)(aOsQuery, aQueryInstance, aQueryTemplate);
    }
  } else {
    uint lValB;
    for(uint lValA = _lo_a; lValA <= _hi_a; lValA += _step_a) {
       for(uint lWidth = _lo_w; lWidth <= _hi_w; lWidth += _step_w) {
         lValB = lValA + lWidth;
         if(_hi_a < lValB) { break; }
         a = lValA;
         b = lValB;
         (*aGenFun)(aOsQuery, aQueryInstance, aQueryTemplate);
       }
    }
  }
}

void
loopspec_t::eval_estimates(checkfun_t aCheckFun, uint& a, uint& b, 
                           const class SimpleProfileGbH& aSimpleProfile,
                           const class TpchEst&          aTpchEst,
                                       queryinstance_t&  aQueryInstance, 
                           const       querytemplate_t&  aQueryTemplate,
                           const class Cb&               aCb) {
  if(_is_eq) {
    for(uint lValA = _lo_a; lValA <= _hi_a; lValA += _step_a) {
      a = lValA;
      b = lValA;
      (*aCheckFun)(aSimpleProfile, aTpchEst, aQueryInstance, aQueryTemplate, aCb);
    }
  } else {
    uint lValB;
    for(uint lValA = _lo_a; lValA <= _hi_a; lValA += _step_a) {
       for(uint lWidth = _lo_w; lWidth <= _hi_w; lWidth += _step_w) {
         lValB = lValA + lWidth;
         if(_hi_a < lValB) { break; }
         a = lValA;
         b = lValB;
         (*aCheckFun)(aSimpleProfile, aTpchEst, aQueryInstance, aQueryTemplate, aCb);
       }
    }
  }
}

