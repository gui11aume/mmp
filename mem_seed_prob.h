#ifndef _MEM_SEED_PROB_DECLARED_
#define _MEM_SEED_PROB_DECLARED_

// Visible functions.

// Initialization and clean up.
int    set_params_mem_prob (size_t, size_t, double);
void   clean_mem_prob (void);

// Options.
void set_mem_prob_method (int);
void set_mem_prob_max_prcsn (void);
void unset_mem_prob_max_prcsn (void);

// MEM seeding probabilities.
double mem_false_pos(size_t, double, size_t);
#endif
