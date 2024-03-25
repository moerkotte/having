#include "cb.hh"

Cb::Cb() : _pat_or_aggr(), 
           _lim_qerr(2.0), _theta(777),
           _caseNo(0), _isUiDense(false), _calcMeanVar(false), 
           _printBad(false), _trace(false), _help(false),
           _qerr_max(0), _theta_qerr_max(0), _count(0) {
}

Cb::~Cb() {
}

void Cb::help(const bool& x) { _help = x; }
void Cb::trace(const bool& x) { _trace = x; }
void Cb::printBad(const bool& x) { _printBad = x; }
void Cb::isUiDense(const bool& x) { _isUiDense = x; }
void Cb::calcMeanVar(const bool& x) { _calcMeanVar = x; }
void Cb::pat_or_aggr(const std::string& x) { _pat_or_aggr = x; }
void Cb::lim_qerr(const double& x) { _lim_qerr = x; }
void Cb::theta(const double& x) { _theta = x; }
void Cb::caseNo(const uint& x) { _caseNo = x; }
