#ifndef INFRA_INTEGER_DECOMPOSITION_HH
#define INFRA_INTEGER_DECOMPOSITION_HH
#pragma once

#include "util.hh"

/*
 * functions to calculate integer compositions
 * a_1 + ... + a_k = n
 * where a_i in [lb,ub]
 */

// general recurrence, slow for k > 5
uint64_t f_int_comp_rec(const uint k, const uint lb, const uint ub, const uint n);

// general form, relying on recurrence using special cases as much as possible
// TODO
uint64_t f_int_comp_csd_rec(const uint k, const uint lb, const uint ub, const uint n);


// cloased form formulas
// a_i in [1,m], n <= m; [Stanley Vol1 p15]
uint64_t f_int_comp_1_m_geq_n(const uint k, const uint lb, const uint ub, const uint n);

// a_i in [0,m], n <= m; [Stanley Vol1 p15]
uint64_t f_int_comp_0_m_geq_n(const uint k, const uint lb, const uint ub, const uint n);


// a_i in [1,m], general
uint64_t f_int_comp_1_m(const uint k, const uint lb, const uint ub, const uint n);


/*
 *  special cases for k=2 and k=3
 *  generally applicable (i.e, ub < n or ub >= n does not matter)
 */

// k = 2
uint64_t f_int_comp_k_eq_2_lb_eq_1(const uint ub, const uint n); // closed form

// k = 3
uint64_t f_int_comp_k_eq_3_lb_eq_1_sum(const int ub, const int n); // summation
uint64_t f_int_comp_k_eq_3_lb_eq_1_csd(const int ub, const int n); // closed form




/*
 * functions to calculate integer decompositions less than n
 * a_1 + ... + a_k <= n
 */

// a_i in [0,m], n <= m; [Stanley Vol1 p15]
// untested
uint64_t f_int_comp_0_m_geq_n_leq(const uint k, const uint lb, const uint ub, const uint n);




#endif

