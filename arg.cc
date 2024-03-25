#include "arg.hh"

void
construct_arg_desc(argdesc_vt& x) {
  // typedef argdesc_t<Cb, int>         iarg_t;
  typedef argdesc_t<Cb, uint>        uarg_t;
  // typedef argdesc_t<Cb, float>       farg_t;
  typedef argdesc_t<Cb, double>      darg_t;
  typedef argdesc_t<Cb, bool>        barg_t;
  // typedef argdesc_t<Cb, std::string> sarg_t;
  // typedef argdesc_t<Cb, char>        carg_t;

  x.push_back( new barg_t("--help",  false, &Cb::help, "print this message" ) );
  x.push_back( new barg_t("-h",      false, &Cb::help, "print this message" ) );
  x.push_back( new barg_t("--trace", false, &Cb::trace, "trace estimation, print single results") );
  x.push_back( new barg_t("-t",      false, &Cb::trace, "trace estimation, print single results") );
  x.push_back( new barg_t("-b",      false, &Cb::printBad, "print bad estimates with q-error > limit (see -q)") );
  x.push_back( new barg_t("-i",          0, &Cb::isUiDense, "assume attribute B is unsigned integer and dense") );
  x.push_back( new barg_t("-c",          0, &Cb::calcMeanVar, "calculate mean/variance for uniform dist") );
  x.push_back( new uarg_t("-n",          0, &Cb::caseNo,   "case number") );
  x.push_back( new darg_t("-q",        2.0, &Cb::lim_qerr, "max. allowed q-error") );
  x.push_back( new darg_t("-theta",  777.0, &Cb::theta,    "theta") );
}

