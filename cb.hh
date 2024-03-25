#ifndef OPT_CARD_HAVING_CB_HH
#define OPT_CARD_HAVING_CB_HH

#include "util.hh"

class Cb {
  private:
    Cb(const Cb&) = delete;
    Cb& operator=(const Cb&) = delete;
  public:
   Cb();
   ~Cb();
  public:
    inline bool   help() const { return _help; }
    inline bool   trace() const { return _trace; }
    inline bool   printBad() const { return _printBad; }
    inline bool   isUiDense() const { return _isUiDense; }
    inline bool   calcMeanVar() const { return _calcMeanVar; }
    inline uint   caseNo() const { return _caseNo; }
    inline double lim_qerr() const { return _lim_qerr; }
    inline double theta() const { return _theta; }
    inline const std::string& pat_or_aggr() const { return _pat_or_aggr; }
  public:
    void help(const bool&);
    void trace(const bool&);
    void printBad(const bool&);
    void isUiDense(const bool&);
    void calcMeanVar(const bool&);
    void caseNo(const uint&);
    void lim_qerr(const double&);
    void theta(const double&);
    void pat_or_aggr(const std::string&);
  private:
    std::string _pat_or_aggr;
    double _lim_qerr;
    double _theta;
    uint   _caseNo;
    bool   _isUiDense;
    bool   _calcMeanVar;
    bool   _printBad;
    bool   _trace;
    bool   _help;
  public:
    mutable double _qerr_max;
    mutable double _theta_qerr_max;
    mutable uint   _count;
};


#endif
