#ifndef OPT_CARD_HAVING_INT_COMP_HH
#define OPT_CARD_HAVING_INT_COMP_HH
#pragma once

#include "util.hh"
#include "integer_composition.hh"

/*
 *  calculate the number of possibilities to 
 *  compose a positive integer n into a sum of k integers a_i in [1,m], such that
 *  sum_{i=1}^k a_i = n
 */


// calculate the number of $k$-compositions of n with a_i in [1,m]
// using a direct implementation of the recurrence
uint64_t no_comp_1(const uint64_t k, const uint64_t n, const uint64_t m);
// using explicit formulas, only possible if 
uint64_t no_comp_2(const uint k, const uint64_t n, const uint64_t m);

#endif
