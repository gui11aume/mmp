#ifndef _SESAME_DECLARED_ 
#define _SESAME_DECLARED_ 

// -- BASIC FUNCTIONS -- //

// Initialization and clean up.
int sesame_set_static_params(size_t, size_t, double);
void sesame_clean(void);

// Automatic functions.
double auto_exact_seed_nullp (int, double, int);
double auto_exact_seed_offp (int, double, int);
double auto_mem_seed_nullp (int, double, int);
double auto_mem_seed_offp (int, double, int);
double auto_skip_seed_nullp (int, int, double, int);
double auto_skip_seed_offp (int, int, double, int);


// -- ADVANCED FUNCTIONS -- //

// Set precision.
void sesame_set_epsilon(double);

// Compute probabilities.
double * exact_seed_null (double, int);
double * exact_seed_offp (double, int);
double * mem_seed_null (double, int);
double * mem_seed_offp (double, int);
double * mem_seed_offp_mcmc (double, int);
double * skip_seed_null (int, double, int);
double * skip_seed_offp (int, double, int);


// Memoization, hashing, storing, loading.
double * fetch_prob (int, double, int);
int store_prob (int, double, int, double *);
void dump_prob_to_file (FILE *);
void load_prob_from_file (FILE *);
void clean_prob_storage (void);
#endif
