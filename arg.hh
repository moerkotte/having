#ifndef OPT_CARD_HAVING_ARG_HH
#define OPT_CARD_HAVING_ARG_HH

#include "argbase.hh"
#include "cb.hh"

typedef std::vector<argdescbase_t<Cb>*> argdesc_vt;

void construct_arg_desc(argdesc_vt& aArgDesc);


#endif

