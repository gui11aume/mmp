#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

// SECTION 1. MACROS //

#define LIBNAME "sesame"
#define VERSION "0.9 21-05-2019"

// Overwrite at compile time for other values.

#ifndef DEFAULT_EPSILON
// Tolerance on the estimates (execpt MCMC).
#define DEFAULT_EPSILON 0.01
#endif

#ifndef DEFAULT_SAMPLINGS
// Number of MCMC samples.
#define DEFAULT_SAMPLINGS 10000000
#endif

#ifndef HSIZE
// Number of MCMC samples.
#define HSIZE 4096
#endif

#define SUCCESS 1
#define FAILURE 0

// Modulo operation that returns a positive number.
#define modulo(x,y) ((((int)x) % ((int)y) + ((int)y)) % ((int)y))

// Check if a polynomial is null.
#define iszero(poly) \
   ((poly) == NULL || ((poly)->monodeg == 0 && (poly)->coeff[0] == 0))


// NOTE 1.1. Error handling. //
//
// All the functions that may fail at runtime -- mostly the function that
// directly or indirectly call 'malloc()' -- have a 'in_case_of_failure'
// section after the final 'return' statement. This section contains the
// instructions for handling errors and exceptions without crashing. In
// general, this means freeing the allocated memory and returning a
// special value that will inform the caller that something went wrong,
// so that the error can be further propagated.
//
// The macro below simplifies error handling for memory errors. See
// function 'warning()' defined in section 4.1.

#define handle_memory_error(x) do { \
   if ((x) == NULL) { \
      warning("memory_error", __func__, __LINE__); \
      goto in_case_of_failure; \
}} while (0)



// SECTION 2. DECLARATIONS  //

// SECTION 2.1 TYPE DECLARATIONS  //
 
// NOTE 2.1.1 'trunc_pol_t' and monomials //
//
// The most important type declared here is the truncated polynomial
// 'trunc_pol_t'. It consists of a 'monodeg' member (of type 'size_t')
// and a 'coeff' member (variable-size array of type 'double').
//
// The member 'monodeg' encodes whether the polynomial is a monomial and
// if so, what its degree is. This is important to gain speed while
// multiplying truncated polynomials. If 'monodeg' is a number from 0 to
// 'K' (the maximum degree of every truncated polynomial), the polynomial
// is a monomial of degree 'monodeg'. If the value is greater than 'K',
// the polynomial is not a monomial.
//
// The member 'coeff' contains the coefficients of the truncated
// polynomial. The encoding is natural in the sense that 'coeff[i]' is
// the coefficient of degree 'i'. If the polynomial is a monomial (i.e.
// 'monodeg' is less than or equal to 'K'), then only the entry
// 'coeff[monodeg]' may be non-zero. If this is not the case (which
// should not happen for the sake of consistency), the terms of degree
// different from 'monodeg' are ignored.
//
// It will not cause problems to have a truncated polynomial with only
// one non-zero coefficient and where 'monodeg' is set to a value higher
// than 'K'. The multiplications are somewhat slower, but the results
// will be consistent.
//
// The product of two monomials is a monomial. The addition of two
// monomials is a monomial only if the monomials have the same degree.
//
// The null struct of type 'trunc_poly_t' (i.e. with all members set to
// 0) is a monomial of degree 0 whose coefficient is 0 (it is thus a
// proper definition of the null polynomial). For addition and
// multiplication of polynomials, a pointer of 'trunc_pol_t' that is set
// to 'NULL' is considered to be a pointer to the null polynomial.

typedef struct matrix_t     matrix_t;      // Transfer matrices.
typedef struct rec_t        rec_t;         // Hash table records.
typedef struct trunc_pol_t  trunc_pol_t;   // Truncated polynomials.

struct rec_t {
   int      n;      // Skipping.
	double   u;      // Divergence rate.
	int      N;      // Number of duplicates (log2).
	double * prob;   // False positive probabilities.
   rec_t  * next;   // Next record if any.
};

struct matrix_t {
   const size_t    dim;  // Column / row number.
   trunc_pol_t * term[]; // Terms of the matrix.
};

struct trunc_pol_t {
   size_t monodeg;       // Monomial (see note 2.1.1).
   double coeff[];       // Terms of the polynomial.
};


// SECTION 2.2 FUNCTION DECLARATIONS  //

// Mersenne twister functions (see license below) //
void   seedMT(unsigned long);
double runifMT(void);
// Function from the R code (see license below) //
size_t rbinom(size_t, double);



// SECTION 3. GLOBALS //

// Precision options.
static double EPSILON = DEFAULT_EPSILON;
static size_t MCMC_SAMPLINGS = DEFAULT_SAMPLINGS;

// NOTE 3.1 Static parameters. //
//
// The static parameters of the seeding problem are 'G' (a.k.a gamma, the
// minimum size of MEM seeds), 'K' (the maximum degree of all truncated
// polynomial, standing for the longest possible sequencing read) and 'P'
// (the probability of a sequencing error). They define the properties of
// the sequencing instrument and the choice of a minimum size by the
// user. The other parameters (actual read size, number of duplicates and
// rate of divergence) can depend on the read under consideration and are
// called "dynamic". Static parameters are set by initializing the
// library with 'sesame_set_static_params()'.

static int G = 0;        // Minimum size of MEM seeds.
static int K = 0;        // Max degree of 'trunc_pol_t' (read size).
static double P = 0.0;   // Probability of a read error.

static size_t  KSZ = 0;     // Memory size of the 'trunc_pol_t' struct.

// Define 6 generic hash tables to store the probabilities. Hashes 'H'
// are private; hash 'Y1' is public.
static rec_t * H1N[HSIZE] = {0}; // Exact seeds, null.
static rec_t * H1O[HSIZE] = {0}; // Exact seeds, null.
static rec_t * H2N[HSIZE] = {0}; // Exact seeds, off target.
static rec_t * H2O[HSIZE] = {0}; // Exact seeds, off target.
static rec_t * H3N[HSIZE] = {0}; // Skip seeds, null.
static rec_t * H3O[HSIZE] = {0}; // Skip seeds, null.
static rec_t * Y1[HSIZE] = {0}; // Skip seeds, off target.

static trunc_pol_t * TEMP = NULL;  // For matrix multipliciation.

static double * XI;     // Precomputed values of 'xi'.
static double * XIc;    // Precomputed values of '1-xi'.
static double * ETA;    // Precomputed values of 'eta'.
static double * ETAc;   // Precomputed values of '1-eta'.

static int PARAMS_INITIALIZED = 0;  // Parameters bookkeeping

// The following should never appear. This is only to give the end users
// a means to contact us and hopefully fix the problems with them.
const char internal_error[] =
   "internal error (please contact guillaume.filion@gmail.com)";



// SECTION 4. FUNCTION DEFINITIONS //

// SECTION 4.1. GET AND SET FUNCTIONS  (VISIBLE) //

void sesame_set_epsilon (double p) { EPSILON = p < 0 ? 0.0 : p; }

double *
to_array_of_double
(
   trunc_pol_t * x
)
// Transforms a 'trunc_pol_t' to an array of double.
{
   return (double *) memmove(x, x->coeff, K * sizeof(double));
}

// SECTION 4.2. ERROR HANDLING FUNCTIONS //

void
warning
(
   const char * msg,
   const char * function,
   const int    line
)
// SYNOPSIS:
//   Print warning message to stderr. This function is used to report
//   errors using the error macros (see macros defined in section 1 and
//   see note 1.1 about error handling).
{
   fprintf(stderr, "[%s] error in function `%s' (line %d): %s\n",
         LIBNAME, function, line, msg);
}


// SECTION 4.3. MATHEMATICAL FUNCTIONS //

double
xi
(
   const size_t i,
   const double u
)
{

   return 1 - pow(1-u, i);

}


double
omega
(
   const size_t m,
   const double u,
   const int N
)
{

   if (m > N) {
      return 0.0;
   }

   double log_N_choose_m = lgamma(N+1)-lgamma(m+1)-lgamma(N-m+1);
   return exp(log_N_choose_m + (N-m)*log(1-u/3) + m*log(u/3));

}


double
psi
(
   const size_t i,
   const size_t m,
   const size_t n,
   const size_t r,
   const double u,
   const int N
)
{

   if (N < (m+r)) return 0.0;

   // Take out the terms of the sum that do not depend on 'q'.
   const double lC = lgamma(m+1) + lgamma(N-m-r+1);
   const double lx = log(xi(i,u));
   double val = 0.0;
   size_t topq = n-r < m ? n-r : m;
   for (int q = 0 ; q <= topq ; q++) {
      val += exp((n-r-q) * lx - lgamma(q+1) - lgamma(m-q+1)
                     - lgamma(n-r-q+1) - lgamma(N-m-n+q+1) + lC);
   }

   return val;

}


double
zeta
(
   const size_t i,
   const size_t m,
   const size_t n,
   const double u,
   const int N
)
{

   if (N < m) return 0.0;

   // Take out the terms of the sum that do not depend on 'r'.
   const double lC = lgamma(N-m+1)+lgamma(n+1)+lgamma(N-n+1)-lgamma(N+1);
   const double lu = i * log(1-u);
   double val = 0.0;
   size_t topr = N-m < n ? N-m : n;
   for (int r = 1 ; r <= topr ; r++) {
      val += psi(i,m,n,r,u,N) * \
             exp(r*lu - lgamma(r+1) - lgamma(N-m-r+1) + lC);
   }

   return val;

}



// SECTION 4.4. INITIALIZATION, BOOKKEEPING AND CLEAN UP FUNCTIONS //

// SECTION 4.4.1. DESTRUCTORS //

void
destroy_mat
(
   matrix_t * mat
)
// SYNOPSIS:
//   Free the memory for all entries of the matrix passed as argument.
{

   // Do nothing if 'mat' is the NULL pointer.
   if (mat == NULL) return;

   size_t nterms = (mat->dim) * (mat->dim);
   for (size_t i = 0 ; i < nterms ; i++)
      free(mat->term[i]);
   free(mat);

}


void
clean_hash
(
   rec_t ** hash
)
// SYNOPSIS:
//   Free the memory for all entries of the hash.
{
   for (int i = 0 ; i < HSIZE ; i++) {
      for (rec_t *rec = hash[i] ; rec != NULL ; ) {
         rec_t *next = rec->next;
         free(rec->prob);
         free(rec);
         rec = next;
      }    
      hash[i] = NULL;
   }
}



// SECTION 4.4.2. CONSTRUCTORS //

trunc_pol_t *
new_zero_trunc_pol
(void)
// SYNOPSIS:
//   Allocate memory for a new struct of type 'trunc_pol_t', initialize
//   it to 0 (i.e. null polynomial).
//
// RETURN:
//   A pointer to the new struct of type 'trunc_pol_t' or 'NULL' in case
//   of failure.
//
// FAILURE:
//   Fails if parameters are not initialized (if 'K' is unknown, we
//   cannot call 'malloc()' with the proper size) or if there is a
//   memory error.
{

   // Cannot be called before 'sesame_set_static_params()'.
   if (!PARAMS_INITIALIZED) {
      warning("parameters unset: call `sesame_set_static_params'",
            __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = calloc(1, KSZ);
   handle_memory_error(new);

   return new;

in_case_of_failure:
   return NULL;

}


matrix_t *
new_null_matrix
(
   const size_t dim
)
// SYNOPSIS:
//   Allocate memory for a new struct of type 'matrix_t', initialize all
//   entries to 0 (i.e. all are NULL pointers). Since the entries of
//   matrix are pointers to 'trunc_pol_t', each of them corresponds to a
//   null polynomial. Still, the entries are not writable because no
//   memory has been allocatd for them (see 'new_zero_matrix()').
//
// RETURN:
//   A pointer to the new struct of type 'matrix_t' or 'NULL' in case
//   of failure.
//
// FAILURE:
//   Fails if there is a memory error.
{

   // Initialize to zero.
   size_t sz = sizeof(matrix_t) + dim*dim * sizeof(trunc_pol_t *);
   matrix_t *new = calloc(1, sz);
   handle_memory_error(new);

   // The dimension is set upon creation and never changes afterwards.
   // In the declaration of the struct, the member 'dim' is a constant,
   // the syntax below is way to go around the 'const' qualifier.
   *(size_t *)&new->dim = dim;

   return new;

in_case_of_failure:
   return NULL;

}

matrix_t *
new_zero_matrix
(
   const size_t dim
)
// SYNOPSIS:
//   Allocate memory for a new struct of type 'matrix_t', initialize all
//   entries to different null polynomials (i.e. all are non NULL
//   pointers to struct of type 'trunc_pol_t' with all their coefficients
//   set to 0). Unlike 'new_null_matrix()', the entries are writable.
//
// RETURN:
//   A pointer to the new struct of type 'matrix_t' or 'NULL' in case
//   of failure.
//
// FAILURE:
//   Fails if parameters are not initialized (if 'K' is unknown, we
//   cannot call 'malloc()' with the proper size for the truncated
//   polynomials) or if there is a memory error. Checking that parameters
//   are set is done indirectly in the call to 'new_zero_trunc_pol()'.
{

   matrix_t *new = new_null_matrix(dim);
   handle_memory_error(new);

   // NOTE: cannot be called before 'sesame_set_static_params()'.
   for (int i = 0 ; i < dim*dim ; i++) {
      handle_memory_error(new->term[i] = new_zero_trunc_pol());
   }

   return new;

in_case_of_failure:
   destroy_mat(new);
   return NULL;

}


// SECTION 4.4.3 HASH MANIPULATION FUNCTIONS //

size_t
squish
(
   int N
)
// SYNOPSIS:
//   Assign 'N' to predefined buckets.
//
// RETURN:
//   The reference value for 'N'.
{

   // NOTE 4.4.3.1. Shared buckets for the values of N.
   //
   // The value of log2(N+1) changes at N = 0, 1, 3, 7, 15, 31, 63, 127,
   // 255, 511, 1023, 2047, 4095, 8191, 16383. In other words, 0 has its
   // own bucket, 1 and 2 have the same bucket, so do 3, 4, 5, 6 etc.
   // The normalized values of 'N' are the mid-points of the intervals.

   const size_t normalized_values[] =
      {0,1,4,10,22,46,94,190,382,766,1534,3070,6142,12286,24574};

   size_t lgN = N < 16383 ? floor(log2(N+1)) : 14;
   return normalized_values[lgN];

}


rec_t *
lookup
(
          rec_t  ** hash,
    const int       n,
    const double    u,
    const int       N
)
// SYNOPSIS:
//   Look for the entry of the hash associated with the key (N,u).
//   'u' is "coarse-grained" to 0.01.
//
// RETURN:
//   A pointer to the struct of type 'rec_t' that corresponds to the key
//   (or its coarse-grained variants) if present, or NULL otherwise.
{

   size_t coarse_u = (100 * u);
   size_t addr = (739*n + 37*N + coarse_u) % HSIZE;
   
   for (rec_t *rec = hash[addr] ; rec != NULL ; rec = rec->next) {
      if (rec->n == n && rec->u == coarse_u && rec->N == N) return rec;
   }

   // Entry not found.
   return NULL;
   
}


rec_t *
insert
(
          rec_t  ** hash,
    const int       n,
    const double    u,
    const int       N,
          double  *  prob
)
// SYNOPSIS:
//   Create an entry in the hash associated with the key (N,u).
//   'u' is "coarse-grained" to 0.01.
//
// RETURN:
//   A pointer to the struct of type 'rec_t' that corresponds to the
//   inserted key, or NULL in case of failure.
//
// FAILURE:
//   Fails if 'malloc()' fails.
{

   // See note 4.4.3.2.
   size_t coarse_u = (100 * u);
   size_t addr = (739*n + 37*N + coarse_u) % HSIZE;

   rec_t *new = malloc(sizeof(rec_t));
   handle_memory_error(new);

   new->n = n;
   new->u = coarse_u;
   new->N = N;
   new->prob = prob;
   new->next = hash[addr];

   hash[addr] = new;

   return new;

in_case_of_failure:
   return NULL;

}


double *
fetch_prob   // VISIBLE // 
(
    const int    n,
    const double u,
    const int    N
)
// SYNOPSIS:
//   Look for the entry of the hash associated with the key (N,u).
//   'u' is "coarse-grained" to 0.01.
//
// RETURN:
//   An array of numbers that corresponds to the key  (or its
//   coarse-grained variants) if present, or NULL otherwise.
{

   rec_t *record = lookup(Y1, n, u, N);
   return record == NULL ? NULL : record->prob;

}



int
store_prob   // VISIBLE // 
(
    const int      n,
    const double   u,
    const int      N,
          double * prob
)
// SYNOPSIS:
//   Create an entry in the hash associated with the key
//   'u' is "coarse-grained" to 0.01.
//
// RETURN:
//   'SUCCESS' (1) in case of success, 'FAILURE' (0) in case of
//   failure.
//
// FAILURE:
//   Fails if 'insert()' fails.
{

   rec_t *record = insert(Y1, n, u, N, prob);
   return record == NULL ? FAILURE : SUCCESS;

}



// SECTION 4.4.4. LIBRARY INITIALIATION AND CLEAN UP //


int
dynamic_params_OK
(
   const int k,    // Segment or read size.
   const double u,    // Divergence rate.
   const int N     // Number of duplicates.
)
// SYNOPSIS:
//   Check if dynamic parameters (see note 3.1) are conform to
//   specifications.
//
// RETURN:
//   SUCCESS (1) upon success or FAILURE (0) upon failure.
//
// FAILURE:
//   Fails if dynamic parameters do not conform specifications or if
//   static parameters are not initialized.
{

   // Check if sequencing parameters were set. Otherwise
   // warn user and fail gracefully (return nan).
   if (!PARAMS_INITIALIZED) {
      warning("parameters unset: call `sesame_set_static_params'",
            __func__, __LINE__);
      return FAILURE;
   }

   if (k > K) {
      char msg[128];
      snprintf(msg, 128,
            "argument k (%d) greater than set value (%d)", k, K);
      warning(msg, __func__, __LINE__);
      return FAILURE;
   }

   if (u <= 0.0 || u >= 1.0) {
      warning("parameter u must be between 0 and 1",
            __func__, __LINE__); 
      return FAILURE;
   }

   return SUCCESS;

}

void
sesame_clean // VISIBLE //
(void)
// SYNOPSIS:
//   Rest global parameters to their uninitialized state and free
//   allocated memory. The function must be called after using the
//   library to prevent memory leaks.
//
// SIDE-EFFETS:
//   Change the values of the global constants 'G', 'K', 'P', 'KSZ'
//   and 'PARAMS_INITIALIZED'. Free global variables 'TEMP', 'XI',
//   'XIc', 'ETA' and 'ETAc'. Destroy global hashes.
{

   PARAMS_INITIALIZED = 0;

   // Set global variables to "fail" values.
   G = 0;
   K = 0;
   P = 0.0;
   KSZ = 0;

   PARAMS_INITIALIZED = 0;

   // Free everything.
   free(TEMP);  TEMP = NULL;
   free(XI);    XI   = NULL;
   free(XIc);   XIc  = NULL;
   free(ETA);   ETA  = NULL;
   free(ETAc);  ETAc = NULL;

   // Clean hash tables.
   clean_hash(H1N);
   clean_hash(H1O);
   clean_hash(H2N);
   clean_hash(H2O);
   clean_hash(H3N);
   clean_hash(H3O);
   clean_hash(Y1);
   
   return;

}


int
sesame_set_static_params // VISIBLE //
(
   int g,
   int k,
   double p
)
// SYNOPSIS:
//   Initialize static parameters from user-defined values.
//
// RETURN:
//   SUCCESS (1) upon success or FAILURE (0) upon failure.
//
// SIDE-EFFETS:
//   Change the values of the global constants 'G', 'K', 'P', 'KSZ'
//   and 'PARAMS_INITIALIZED'. Destroy global hashes. Reset the
//   global 'TEMP'.
//
// FAILURE:
//   Fails if static parameters do not conform specifications or if
//   'malloc()' fails (the truncated polynomial 'TEMP' must be
//   allocated).
{

   // Check input
   if (g == 0 || k == 0) {
      warning("parameters g and k must greater than 0",
            __func__, __LINE__); 
      goto in_case_of_failure;
   }

   if (p <= 0.0 || p >= 1.0) {
      warning("parameter p must be between 0 and 1",
            __func__, __LINE__); 
      goto in_case_of_failure;
   }

   G = g;  // MEM size.
   K = k;  // Read size.
   P = p;  // Sequencing error.

   // All 'trunc_pol_t' must be congruent.
   KSZ = sizeof(trunc_pol_t) + (K+1) * sizeof(double);

   PARAMS_INITIALIZED = 1;

   // Clean previous values (if any).
   clean_hash(H1N);
   clean_hash(H1O);
   clean_hash(H2N);
   clean_hash(H2O);
   clean_hash(H3N);
   clean_hash(H3O);
   clean_hash(Y1);

   free(TEMP); // Nothing happens if 'TEMP' is NULL.
   handle_memory_error(TEMP = new_zero_trunc_pol());

   return SUCCESS;

in_case_of_failure:
   sesame_clean();
   return FAILURE;

}



// SECTION 4.5. WEIGHTED GENERATING FUNCTIONS //

// NOTE 4.5.1. Internal errors in weighted generating functions //
// The functions below do not have detailed synopsis because they all
// follow the same pattern of allocating a particular polynomial defined
// in the theory of computing MEM seeding probabilities by the symbolic
// method. Some functions check the value of the input paramters and can
// trigger an "internal error" (see note 1.1), whereas others do not. The
// logic is that the polynomials can be mathematically undefined (when
// the expression involves a division by zero for instance) or simply not
// used in the present theory. The former case triggers an internal
// error, whereas the second does not. This is to allow anyone to reuse
// the code for their own purpose, without "hurting themselves" by
// triggering undefined, unpredictable or unwanted behaviors. Calling the
// visible functions of this library should never trigger an internal
// error.

trunc_pol_t *
new_trunc_pol_A
(
   const size_t m,    // Initial state (down).
   const size_t n,    // Final state (down).
   const double u,    // Divergence rate.
   const int N     // Number of duplicates.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.
   
   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // When N is 0, these polynomials are non-null
   // only if 'm' and 'n' are zero.
   if (N == 0) {
      if (m == 0 && n == 0) {
         new->monodeg = 1;
         new->coeff[1] = P;
      }
      return new;
   }

   // This is not a monomial.
   new->monodeg = K+1;

   double omega_p_pow_of_q = omega(n,u,N) * P;
   for (int i = 1 ; i <= G ; i++) {
      new->coeff[i] = (1 - pow(xi(i-1,u),N)) * omega_p_pow_of_q;
      omega_p_pow_of_q *= (1.0-P);
   }
   for (int i = G+1 ; i <= K ; i++) {
      new->coeff[i] = (1 - pow(xi(i-1,u),m) *
            (1- zeta(i-1,m,n,u,N))) * omega_p_pow_of_q;
      omega_p_pow_of_q *= (1.0-P);
   }
   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_B
(
   const size_t i,    // Final state (up).
   const double u,    // Divergence rate.
   const int N     // Number of duplicates.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   if (i < 1) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }   

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   if (N == 0) {
      // Zero if i is not 1.
      if (i != 1) return new;
      // Monomial of degree 1 otherwise.
      new->monodeg = 1;
      new->coeff[1] = 1-P;
      return new;
   }

   // This is a monomial.
   new->monodeg = i;
   new->coeff[i] = (pow(xi(i,u),N) - pow(xi(i-1,u),N)) * pow(1-P,i);
   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_C
(
   const size_t m,    // Initial state (down).
   const double u,    // Divergence rate.
   const int N     // Number of duplicates.
)
// Note: the polynomials are defined in the case m > N, but
// they have no meaning in the present context. Here they are
// simply not forbidden, but also not used (and not tested).
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // Special case that 'N' is 0. Assuming that 'm' is 0 (see
   // remark above), this is a monomial equal to 1.
   if (N == 0) {
      new->monodeg = 0;
      new->coeff[0] = 1.0;
      return new;
   }

   // This is not a monomial (except in the case that
   // 'G' is 1 and 'm' is 1, but this is too rare to justify
   // the extra code (it won't trigger a mistake anyway).
   new->monodeg = K+1;

   double pow_of_q = 1.0;
   for (int i = 0 ; i <= G-1 ; i++) {
      new->coeff[i] = (1 - pow(xi(i,u),N)) * pow_of_q;
      pow_of_q *= (1.0-P);
   }
   for (int i = G ; i <= K ; i++) {
      new->coeff[i] = (1 - pow(xi(i,u),m)) * pow_of_q;
      pow_of_q *= (1.0-P);
   }
   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_D
(
   const size_t j,    // Initial state (up).
   const size_t m,    // Final state (down).
   const double u,    // Divergence rate.
   const int N     // Number of duplicates.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   if (j > G-1) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }   

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // In the case that 'N' is 0, the only polynomial that is
   // defined is D(1,0,0) -- the others are set to 0. In the
   // case that 'm' is greater than 'N' the polynomial is 0.
   if ((N == 0 && !(j == 1 && m == 0)) || m > N) {
      return new;
   }

   // This is a monomial only if j is equal to G-1.
   new->monodeg = j == G-1 ? 1 : K+1;

   double omega_p_pow_of_q = omega(m,u,N) * P;
   for (int i = 1 ; i <= G-j ; i++) {
      new->coeff[i] = omega_p_pow_of_q;
      omega_p_pow_of_q *= (1.0-P);
   }

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_E
(
   const size_t j    // Initial state (up).
)
{

   if (j > G-1) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // This is a monomial only if j is equal to G-1
   new->monodeg = j == G-1 ? 0 : K+1;

   double pow_of_q = 1.0;
   for (int i = 0 ; i <= G-j-1 ; i++) {
      new->coeff[i] = pow_of_q;
      pow_of_q *= (1.0-P);
   }

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_F
(
   const size_t j,
   const double u    // Divergence rate.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   if (j > G-1) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // This is a monomial only when 'i' is 0.
   new->monodeg = j == 0 ? 0 : K+1;

   const double a = (1-P)*(1-u);
   double pow_of_a = 1.0;
   for (int i = 0 ; i <= j ; i++) {
      new->coeff[i] = pow_of_a;
      pow_of_a *= a;
   }

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_H
(
   const size_t r,   // Source phase.
   const size_t s,   // Target phase.
   const size_t n    // Skipping. XXX Make a global constant???
)
{

   if (s > n || r > n) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   const int x = modulo((int)r-(int)s-1, (int)n+1);
   const int m = (G-1+r-x) / (n+1);

   // This is a monomial when 'm' is zero.
   new->monodeg = m == 0 ? x+1 : K+1;

   double pqx_times_the_rest = P*pow(1-P,x);
   for (int i = 0 ; i <= m ; i++) {
      if (x+1+i*(n+1) > K) break;
      new->coeff[x+1+i*(n+1)] = pqx_times_the_rest;
      pqx_times_the_rest *= pow(1-P,n+1);
   }

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_J
(
   const size_t r,  // Source phase.
   const size_t n   // Skipping.
)
{

   if (r > n) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   new->monodeg = K+1;

   double pow_of_q = 1.0;
   const int imax = K > G-1+r ? G-1+r : K;
   for (int i = 0 ; i <= imax ; i++) {
      new->coeff[i] = pow_of_q;
      pow_of_q *= 1-P;
   }

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_N
(
   size_t m
)
{

   if (m > K) m = K;

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   if (m > 0) new->monodeg = K+1;

   for (int i = 0 ; i <= m ; i++) {
      new->coeff[i] = 1.0;
   }

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_R
(
   const size_t j,
   const double u    // Divergence rate.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   if (j > G-1) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // This is a monomial only when 'i' is 0.
   new->monodeg = j == 0 ? 1 : K+1;

   const double a = (1-P)*(1-u);
   const double d = P*(1-u/3.0);
   double pow_of_a = 1.0;
   for (int i = 0 ; i <= j ; i++) {
      new->coeff[i+1] = pow_of_a * d;
      pow_of_a *= a;
   }

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_r_plus
(
   const size_t i,
   const double u    // Divergence rate.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   if (i > G-2) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // This is a monomial.
   new->monodeg = i+1;
   new->coeff[i+1] = i == 0 ?
      P * u/3.0 : pow((1-P)*(1-u), i) * P * u/3.0;

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_r_minus
(
   const size_t i,
   const double u    // Divergence rate.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   if (i > G-2) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // This is a monomial.
   new->monodeg = i+1;
   new->coeff[i+1] = i == 0 ?
      (1-P)*u : pow((1-P)*(1-u), i) * (1-P)*u;

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_ss
(
   const int i,  // Match length.
   const int j,  // Destination phase.
   const int n,  // Skipping.
   const double u   // Divergence rate.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   if (i > G-1 || j > G-1 || i < 1 || j < 1) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   const int x = modulo(n-i+1, n+1);

   // Polynomial is null in the following case.
   if (i+x+j > G) return new;

   const double a = (1-P)*(1-u);
   const double b = (1-P)*u;

   // This is a monomial when 'm' is zero.
   new->monodeg = x+j;
   new->coeff[x+j] = b*pow(a,x+j-1);

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_tt
(
   const int i,  // Match length.
   const int j,  // Destination phase.
   const int n,  // Skipping.
   const double u   // Divergence rate.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   if (i > G-1 || j > G-1 || i < 1 || j < 1) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   const int x = modulo(n-i+1, n+1);

   // Polynomial is null in the following case.
   if (i+x+j > G) return new;

   const double a = (1-P)*(1-u);
   const double c = P*u/3;

   // This is a monomial when 'm' is zero.
   new->monodeg = x+j;
   new->coeff[x+j] = c*pow(a,x+j-1);

   return new;

in_case_of_failure:
   return NULL;

}



trunc_pol_t *
new_trunc_pol_U
(
   const int i,  // Match length.
   const int j,  // Destination phase.
   const int n,  // Skipping.
   const double u   // Divergence rate.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   if (j > n || i > G-1 || i <= 0 || j < 0) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   const int x = modulo(-(i+j+1), n+1);
   const int y = modulo(n-i+1, n+1);
   const int m = (G-1-i-x) / (n+1);

   // Polynomial is null in the following case.
   if (i+x > G-1) return new;

   const double a = (1-P)*(1-u);
   const double b = (1-P)*u;
   const double d = P*(1-u/3); 

   // This is a monomial when 'm' is zero.
   new->monodeg = m == 0 ? x+1 : K+1;

   double dax_times_the_rest = d*pow(a,x);
   for (int r = 0 ; r <= m ; r++) {
      new->coeff[x+1+r*(n+1)] = dax_times_the_rest;
      dax_times_the_rest *= pow(a,n+1);
   }

   if (j < y) new->coeff[x+1] += b*pow(a,x);

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_V
(
   const int i,  // Source phase.
   const int j,  // Destination phase.
   const int n,  // Skipping.
   const double u   // Divergence rate.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   if (j > n || i > G-1 || i == 0 || j < 0) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   const int x = modulo(-(i+j+1), n+1);
   const int y = modulo(n-i+1, n+1);
   const int m = (G-1-i-x) / (n+1);

   // Polynomial is null in the following case.
   if (i+x > G-1) return new;

   const double a = (1-P)*(1-u);
   const double c = P*u/3;
   const double d = P*(1-u/3); 

   // This is a monomial when 'm' is zero.
   new->monodeg = m == 0 ? x+1 : K+1;

   double dax_times_the_rest = d*pow(a,x);
   for (int r = 0 ; r <= m ; r++) {
      new->coeff[x+1+r*(n+1)] = dax_times_the_rest;
      dax_times_the_rest *= pow(a,n+1);
   }

   if (j < y) new->coeff[x+1] += c*pow(a,x);

   return new;

in_case_of_failure:
   return NULL;

}


trunc_pol_t *
new_trunc_pol_W
(
   const int j,  // Destination phase.
   const int n,  // Skipping.
   const double u   // Divergence rate.
)
{

   // NOTE: 'u' must be in (0,1), we assume the caller checked.

   if (j > n) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   trunc_pol_t *new = new_zero_trunc_pol();
   handle_memory_error(new);

   // In this case the polynomial is null.
   if (j < n-G+1) return new;

   const int x = modulo(n-j, n+1);
   const int m = (G-1-x) / (n+1);

   // This is a monomial when 'm' is zero.
   new->monodeg = m == 0 ? x+1 : K+1;

   const double a = (1-P)*(1-u);
   const double d = P*(1-u/3.0);

   double dax_times_the_rest = d*pow(a,x);
   for (int i = 0 ; i <= m ; i++) {
      new->coeff[x+1+i*(n+1)] = dax_times_the_rest;
      dax_times_the_rest *= pow(a,n+1);
   }

   return new;

in_case_of_failure:
   return NULL;

}



// SECTION 4.6. TRANSFER MATRICES //

matrix_t *
new_matrix_M
(
   const double u,   // Divergence rate.
   const int N    // Number of duplicates.
)
// SYNOPSIS:
//   Allocate memory for a new struct of type 'matrix_t' and initialize
//   the entries (with a mix of NULL and non NULL pointers) to obtain the
//   transfer matrix M(z) of reads without good seed for the
//   given static and dynamic parameters.
//
// RETURN:
//   A pointer to the new struct of type 'matrix_t' or 'NULL' in case
//   of failure.
//
// FAILURE:
//   Fails if there is a memory error.
{

   matrix_t *M = NULL;

   // NOTE 4.6.1. Checking dynamic parameters //
   //
   // In the library, this function is called only by
   // 'wgf_mem()', itself called only by
   // sesame_false_positive_or_negative(),
   // which already checks that dynamic parameters 'k', 'u' and 'N' are
   // conform. If this code is reused somewhere else, it may be a good
   // thing to check them on every call by uncommenting the following
   // code block.
   //
   // -- BEGIN BLOCK -- //
   // Check dynamic parameters. Only 'u' and 'N' must be checked, so
   // we pass 'K' as a first parameter, which cannot trigger an error.
   //if (!dynamic_params_OK(K,u,N)) {
   //   goto in_case_of_failure;
   //}
   // -- END BLOCK -- //
   
   const size_t dim = G+N+1;
   M = new_null_matrix(dim);
   handle_memory_error(M);

   // First N+1 series of rows.
   for (int j = 0 ; j <= N ; j++) {
      for (int i = 0 ; i <= N ; i++) {
         handle_memory_error(
            M->term[j*dim+i] = new_trunc_pol_A(j,i,u,N)
         );
      }
      for (int i = N+1 ; i <= N+G-1 ; i++) {
         handle_memory_error(
            M->term[j*dim+i] = new_trunc_pol_B(i-N,u,N)
         );
      }
      handle_memory_error(
         M->term[j*dim+(N+G)] = new_trunc_pol_C(j,u,N)
      );
   }

   // Next G-1 rows.
   for (int j = N+1 ; j <= N+G-1 ; j++) {
      for (int i = 0 ; i <= N ; i++) {
         handle_memory_error(
            M->term[j*dim+i] = new_trunc_pol_D(j-N,i,u,N)
         );
      }
      handle_memory_error(
         M->term[j*dim+(N+G)] = new_trunc_pol_E(j-N)
      );
   }

   // Last row is null (nothing to do).

   return M;

in_case_of_failure:
   destroy_mat(M);
   return NULL;

}


matrix_t *
new_matrix_L
(
   const double u        // Divergence rate.
)
// SYNOPSIS:
//   Allocate memory for a new struct of type 'matrix_t' and initialize
//   the entries (with a mix of NULL and non NULL pointers) to obtain the
//   transfer matrix L(z) of reads with neither on-target nor off-target
//   exact gamma seed when there is N = 1 duplicate.
//
// RETURN:
//   A pointer to the new struct of type 'matrix_t' or 'NULL' in case
//   of failure.
//
// FAILURE:
//   Fails if there is a memory error.
{

   matrix_t *L = NULL;

   // Assume parameters were checked by the caller.  See note 4.6.1.
   // about checking dynamic parameters if this code is reused.
   // -- BEGIN BLOCK -- //
   // Check dynamic parameters. Only 'u' needs to be checked, so
   // we pass 'K' as a first parameter and 0 as a third parameter,
   // which never generate errors.
   //if (!dynamic_params_OK(K,u,0)) {
   //   goto in_case_of_failure;
   //}
   // -- END BLOCK -- //

   const size_t dim = 2*G;
   L = new_null_matrix(dim);
   handle_memory_error(L);

   // First row.
   handle_memory_error(
      L->term[0] = new_trunc_pol_R(G-1,u)
   );
   for (int j = 0 ; j <= G-2 ; j++) {
      handle_memory_error(
         L->term[j+1] = new_trunc_pol_r_plus(j,u)
      );
   }
   for (int j = 0 ; j <= G-2 ; j++) {
      handle_memory_error(
         L->term[j+G] = new_trunc_pol_r_minus(j,u)
      );
   }
   handle_memory_error(
      L->term[dim-1] = new_trunc_pol_F(G-1,u)
   );

   // Next 'G-1' rows -- matrices A(z) and B(z).
   for (int i = 1 ; i <= G-1 ; i++) {
      handle_memory_error(
         L->term[i*dim] = new_trunc_pol_R(G-1-i,u)
      );
      for (int j = 0 ; j <= ((int)G)-2-i ; j++) {
         handle_memory_error(
            L->term[i*dim+i+j+1] = new_trunc_pol_r_plus(j,u)
         );
      }
      for (int j = 0 ; j <= G-1-i ; j++) {
         handle_memory_error(
            L->term[i*dim+j+G] = new_trunc_pol_r_minus(j,u)
         );
      }
      handle_memory_error(
         L->term[i*dim+dim-1] = new_trunc_pol_F(G-1-i,u)
      );
   }

   // Next G-1 rows -- matrices C(z) and D(z).
   for (int i = G ; i <= 2*G-2 ; i++) {
      handle_memory_error(
         L->term[i*dim] = new_trunc_pol_R(2*G-2-i,u)
      );
      for (int j = 0 ; j <= 2*G-2-i ; j++) {
         handle_memory_error(
            L->term[i*dim+j+1] = new_trunc_pol_r_plus(j,u)
         );
      }
      for (int j = 0 ; j <= (2*(int)G)-3-i ; j++) {
         handle_memory_error(
            L->term[i*dim+j+G+(i-G+1)] = new_trunc_pol_r_minus(j,u)
         );
      }
      handle_memory_error(
         L->term[i*dim+dim-1] = new_trunc_pol_F(2*G-2-i,u)
      );
   }

   // Last row is null (nothing to do).

   return L;

in_case_of_failure:
   destroy_mat(L);
   return NULL;

}


matrix_t *
new_matrix_S
(
   const size_t n       // Skipping.
)
// SYNOPSIS:
//   Allocate memory for a new struct of type 'matrix_t' and initialize
//   the entries (with a mix of NULL and non NULL pointers) to obtain the
//   transfer matrix S(z) of reads without on-target skip-n seed.
//
// RETURN:
//   A pointer to the new struct of type 'matrix_t' or 'NULL' in case
//   of failure.
//
// FAILURE:
//   Fails if static parameters are not initialized or if there is a
//   memory error.
{

   matrix_t *S = NULL;

   const size_t dim = n+2;
   S = new_null_matrix(dim);
   handle_memory_error(S);

   // First n+1 rows.
   for (int i = 0 ; i <= n ; i++) {
      for (int j = 0 ; j <= n ; j++) {
         handle_memory_error(
            S->term[i*dim+j] = new_trunc_pol_H(i,j,n)
         );
      }
      handle_memory_error(
         S->term[i*dim+dim-1] = new_trunc_pol_J(i,n)
      );
   }

   // Last row is null (nothing to do).

   return S;

in_case_of_failure:
   destroy_mat(S);
   return NULL;

}


matrix_t *
new_matrix_T
(
   const size_t n,  // Skipping.
   const double u   // Divergence rate.
)
// SYNOPSIS:
//   Allocate memory for a new struct of type 'matrix_t' and initialize
//   the entries (with a mix of NULL and non NULL pointers) to obtain the
//   transfer matrix T(z) of reads without skip-n seed for either the
//   target or a specific duplicate.
//
// RETURN:
//   A pointer to the new struct of type 'matrix_t' or 'NULL' in case
//   of failure.
//
// FAILURE:
//   Fails if there is a memory error or if the static and dynamic
//   parameters are not properly set.
{

   matrix_t *T = NULL;

   // Check dynamic parameters. Only 'u' must be checked, so
   // we pass 'K' as a first parameter, and 0 as third paramater,
   // which cannot trigger an error.
   if (!dynamic_params_OK(K,u,0)) {
      goto in_case_of_failure;
   }
   
   const size_t dim = n+2*G;
   T = new_null_matrix(dim);
   handle_memory_error(T);

   // First row.
   for (int i = 0 ; i < n+1 ; i++) {
      handle_memory_error(
         T->term[i] = new_trunc_pol_W(i,n,u)
      );
   }
   for (int i = n+1 ; i < n+G ; i++) {
      handle_memory_error(
         T->term[i] = new_trunc_pol_r_plus(i-n-1,u)
      );
   }
   for (int i = n+G ; i < dim-1 ; i++) {
      handle_memory_error(
         T->term[i] = new_trunc_pol_r_minus(i-n-G,u)
      );
   }
   // Last term (tail).
   handle_memory_error(
      T->term[dim-1] = new_trunc_pol_F(G-1, u)
   );

   // Next 'n' rows.
   for (int j = 1 ; j <= n ; j++) {
      // First term (monomial equal to a power of z).
      handle_memory_error(
         T->term[j*dim+0] = new_zero_trunc_pol()
      );
      if (j <= K) {
         T->term[j*dim+0]->coeff[j] = 1.0;
         T->term[j*dim+0]->monodeg = j;
      }
      // Tail term (tail).
      handle_memory_error(
         T->term[j*dim+dim-1] = new_trunc_pol_N(j-1)
      );
   }

   // Next 'G-1' rows.
   for (int j = n+1 ; j < n+G ; j++) {
      for (int i = 0 ; i <= n ; i++) {
         handle_memory_error(
            T->term[j*dim+i] = new_trunc_pol_U(j-n,i,n,u)
         );
      }
      // Matrix A(z).
      for (int i = n+1 ; i <= n+G-1 ; i++) {
         int deg = i-j-1;
         if (deg < 0) continue;
         handle_memory_error(
            T->term[j*dim+i] = new_trunc_pol_r_plus(deg,u)
         );
      }
      // Matrix ~B(z).
      for (int i = n+G ; i <= n+2*G-2 ; i++) {
         const int from = j-n;
         const int to = i-n-G+1;
         handle_memory_error(
            T->term[j*dim+i] = new_trunc_pol_ss(from,to,n,u)
         );
      }
      // Tail term (F polynomial).
      handle_memory_error(
         T->term[j*dim+dim-1] = new_trunc_pol_F(G-1+n-j, u)
      );
   }

   // Next 'G-1' rows.
   for (int j = n+G ; j < dim-1 ; j++) {
      for (int i = 0 ; i <= n ; i++) {
         handle_memory_error(
            T->term[j*dim+i] = new_trunc_pol_V(j-n-G+1,i,n,u)
         );
      }
      // Matrix ~C(z).
      for (int i = n+1 ; i <= n+G-1 ; i++) {
         const int from = j-n-G+1;
         const int to = i-n;
         handle_memory_error(
            T->term[j*dim+i] = new_trunc_pol_tt(from,to,n,u)
         );
      }
      // Matrix D(z).
      for (int i = n+G ; i <= n+2*G-2 ; i++) {
         int deg = i-j-1;
         if (deg < 0) continue;
         handle_memory_error(
            T->term[j*dim+i] = new_trunc_pol_r_minus(deg,u)
         );
      }
      // Tail term (polynomial F).
      handle_memory_error(
         T->term[j*dim+dim-1] = new_trunc_pol_F(2*G-2+n-j, u)
      );
   }

   // Last row is null (nothing to do).

   return T;

in_case_of_failure:
   destroy_mat(T);
   return NULL;

}



// SECTION 4.7. POLYNOMIAL MANIPULATION FUNCTIONS //

trunc_pol_t *
trunc_pol_mult
(
         trunc_pol_t * dest,
   const trunc_pol_t * a,
   const trunc_pol_t * b
)
// SYNOPSIS:
//   Multiply two truncated polynomials 'a' and 'b', and store the result
//   in 'dest'.
//
// RETURN:
//   A pointer to 'dest'.
{

   // Erase 'dest' (set it to zero).
   bzero(dest, KSZ);

   if (iszero(a) || iszero(b)) return dest;

   if (a->monodeg < K+1 && b->monodeg < K+1) {
      // Both are monomials, just do one multiplication.
      // If degree is too high, all coefficients are zero.
      if (a->monodeg + b->monodeg > K) return NULL;
      // Otherwise do the multiplication.
      dest->monodeg = a->monodeg + b->monodeg;
      dest->coeff[dest->monodeg] =
         a->coeff[a->monodeg] * b->coeff[b->monodeg];
      return dest;
   }

   // The result cannot be a monomial.
   dest->monodeg = K+1;
   if (a->monodeg < K+1) {
      // 'a' is a monomial, do one "row" of multiplications.
      for (int i = a->monodeg ; i <= K ; i++)
         dest->coeff[i] = a->coeff[a->monodeg] * b->coeff[i-a->monodeg];
   }
   else if (b->monodeg < K+1) {
      // 'b' is a monomial, do one "row" of multiplications.
      for (int i = b->monodeg ; i <= K ; i++)
         dest->coeff[i] = b->coeff[b->monodeg] * a->coeff[i-b->monodeg];
   }
   else {
      // Standard convolution product.
      for (int i = 0 ; i <= K ; i++) {
         dest->coeff[i] = a->coeff[0] * b->coeff[i];
         for (int j = 1 ; j <= i ; j++) {
            dest->coeff[i] += a->coeff[j] * b->coeff[i-j];
         }
      }
   }

   return dest;

}


void
trunc_pol_update_add
(
         trunc_pol_t * dest,
   const trunc_pol_t * a
)
// SYNOPSIS:
//   Add polynomial 'a' to 'dest' and store the result in 'dest'.
{

   // No update when adding zero.
   if (iszero(a)) return;

   // Only adding two monomials of same degree returns
   // a monomial (the case of adding a zero polynomial
   // was taken care of above.
   if (a->monodeg != dest->monodeg) {
      dest->monodeg = K+1;
   }

   for (int i = 0 ; i <= K ; i++) {
      dest->coeff[i] += a->coeff[i];
   }

}



// SECTION 4.8. MATRIX MANIPULATION FUNCTIONS //

matrix_t *
matrix_mult
(
         matrix_t * dest,
   const matrix_t * a,
   const matrix_t * b
)
// SYNOPSIS:
//   Multiply matrices 'a' and 'b', and store the result in 'dest'.
//
// RETURN:
//   A pointer to 'dest'.
//
// FAILURE:
//   Fails if matrices are not congruent (i.e. square matrices with
//   identical dimensions).
{

   if (a->dim != dest->dim || b->dim != dest->dim) {
      warning(internal_error, __func__, __LINE__);
      goto in_case_of_failure;
   }

   size_t dim = dest->dim;

   for (int i = 0 ; i < dim ; i++) {
   for (int j = 0 ; j < dim ; j++) {
      // Erase destination entry.
      bzero(dest->term[i*dim+j], KSZ);
      for (int m = 0 ; m < dim ; m++) {
         trunc_pol_update_add(dest->term[i*dim+j],
            trunc_pol_mult(TEMP, a->term[i*dim+m], b->term[m*dim+j]));
      }
   }
   }

   return dest;

in_case_of_failure:
   return NULL;

}


// SECTION 4.9. HIGH LEVEL WGF COMPUTATION FUNCTIONS //

trunc_pol_t *
wgf_seed
(void)
// SYNOPSIS:
//   Compute the probabilities that reads do not contain an on-target
//   exact gamma seed for specified static parameters using recurrence
//   formula (12) from doi:10.3390/a11010003 "Analytic Combinatorics for
//   Computing Seeding Probabilities".
//
// RETURN:
//   A pointer to a struct of type 'trunc_pol_t' containing the
//   probabilitie of interest, or NULL in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized or if 'malloc()'
//   fails. Initialization is checked indirectly through the call to
//   'new_zero_trunc_pol()'.
{

   trunc_pol_t *pol = new_zero_trunc_pol();
   handle_memory_error(pol);

   // Not a monomial.
   pol->monodeg = K+1;

   const double q_pow_gamma = pow(1-P,G);
   const double pq_pow_gamma = P * q_pow_gamma;

   for (int i = 0 ; i < G ; i++) pol->coeff[i] = 1.0;
   pol->coeff[G] = 1.0 - q_pow_gamma;
   for (int i = G+1 ; i <= K ; i++) {
      pol->coeff[i] = pol->coeff[i-1] - pq_pow_gamma * pol->coeff[i-G-1];
   }

   return pol;

in_case_of_failure:
   return NULL;

}


// Need this snippet to compute a bound of the numerical imprecision.
double HH(double x, double y)
         {  return x*log(x/y)+(1-x)*log((1-x)/(1-y));  }

trunc_pol_t *
wgf_mem
(
   const double u, // Divergence rate.
   const int N     // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probabilities that reads do not contain an on-target
//   MEM seed for specified static and dynamic parameters.
//
// RETURN:
//   A pointer to a struct of type 'trunc_pol_t' containing the
//   probabilitie of interest, or NULL in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized or if 'malloc()'
//   fails. Initialization is checked indirectly through the call to
//   'new_zero_trunc_pol()'.
{

   // Assume parameters were checked by the caller.  See note 4.6.1.
   // about checking dynamic parameters if this code is reused.
   // -- BEGIN BLOCK -- //
   // Check dynamic parameters. Only 'u' and 'N' must be checked, so
   // we pass 'K' as a first parameter, which cannot trigger an error.
   //if (!dynamic_params_OK(K,u,N)) {
   //   goto in_case_of_failure;
   //}
   // -- END BLOCK -- //

   // The truncated polynomial 'w' stands the weighted generating
   // function of the reads without good seed for set parameters
   // and 'N' duplicates.
   trunc_pol_t *w = new_zero_trunc_pol();

   // The matrix 'M' is the transfer matrix of reads without MEM seed.
   // The row and column of the head state have been deleted because the
   // reads can be assumed to start in state "down/0".  The entries of
   // 'M' are truncated polynomials that stand for weighted generating
   // functions of segments joining different states. The m-th power of
   // 'M' contains the weighted generating functions of sequences of m
   // segments joining the dfferent states. We thus update 'w' with the
   // entry of M^m that joins the initial "down/0" state to the tail
   // state.
   matrix_t *M = new_matrix_M(u, N);

   // Create two temporary matrices 'powM1' and 'powM2' to compute the
   // powers of M. On first iteration, powM1 = M * M, and later perform
   // the operations powM2 = M * powM1 and powM1 = M * powM2, so that
   // powM1 is M^2m at iteration m. Using two matrices allows the
   // operations to be performed without requesting additional memory.
   // The temporary matrices are implicitly erased (not freed) in the
   // course of the multiplication by 'matrix_mult()'.
   matrix_t *powM1 = new_zero_matrix(G+N+1);
   matrix_t *powM2 = new_zero_matrix(G+N+1);

   handle_memory_error(w);
   handle_memory_error(M);
   handle_memory_error(powM1);
   handle_memory_error(powM2);

   // Update weighted generating function with one-segment reads (i.e.
   // tail only). The weighte generating function of tail segments
   // following the state "down/0" is in the top-right entry of 'M', i.e.
   // at 'M[dim]', where dim = G+N is the dimension of the matrix.
   trunc_pol_update_add(w, M->term[G+N]);

   matrix_mult(powM1, M, M);

   // Update weighted generating function with two-segment reads. Same
   // comment as above, 'M[dim]' contains the weighted generating
   // function of tail segments in two-segment reads starting in state
   // "down/0".
   trunc_pol_update_add(w, powM1->term[G+N]);

   // A read with m segments has at least x = floor((m-1)/2) errors.
   // Indeed, a non-error terminator (states "up/i") cannot follow
   // another, so the read must contain at least so many error
   // terminators (the tail segment has no terminator). We bound the
   // probability that a read of size k has more than 'm' segments -- and
   // thus more than x errors -- by a formula for the binomial
   // distribution, where m is the number of segments, i.e. the power of
   // matirx M. For detail see link below.
   // https://en.wikipedia.org/wiki/Binomial_distribution#Tail_Bounds
   for (int n = 2 ; n < K ; n += 2) {
      // Increase the number of segments and update
      // the weighted generating function accordingly.
      matrix_mult(powM2, M, powM1);
      trunc_pol_update_add(w, powM2->term[G+N]);
      matrix_mult(powM1, M, powM2);
      trunc_pol_update_add(w, powM1->term[G+N]);
      // In max precision mode, get all possible digits by computing all
      // the powers of 'M' up to 'K'.
      // Otherwise, stop when reaching 1% precision.
      double x = floor((n-1)/2) / ((double) K);
      double bound_on_imprecision = exp(-HH(x, P)*K);
      if (bound_on_imprecision / w->coeff[K] < EPSILON) break;
   }

   // Clean temporary variables.
   destroy_mat(powM1);
   destroy_mat(powM2);
   destroy_mat(M);

   return w;

in_case_of_failure:
   // Clean everything.
   destroy_mat(powM1);
   destroy_mat(powM2);
   destroy_mat(M);
   free(w);
   return NULL;

}


trunc_pol_t *
wgf_dual
(
   const double u     // Divergence rate.
)
// SYNOPSIS:
//   Compute the probabilities that reads contain neither an on-target
//   exact gamma seed nor an off-target exact gamma seed for the
//   specified static and dynamic parameters, where there is N = 1
//   duplicate of the target.
//
// RETURN:
//   A pointer to a struct of type 'trunc_pol_t' containing the
//   probabilitie of interest, or NULL in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized or if 'malloc()'
//   fails. Initialization is checked indirectly through the call to
//   'new_zero_trunc_pol()'.
{

   // Assume parameters were checked by the caller.  See note 4.6.1.
   // about checking dynamic parameters if this code is reused.
   // -- BEGIN BLOCK -- //
   // Check dynamic parameters. Only 'u' and 'N' must be checked, so
   // we pass 'K' as a first parameter, which cannot trigger an error.
   //if (!dynamic_params_OK(K,u,N)) {
   //   goto in_case_of_failure;
   //}
   // -- END BLOCK -- //

   // The logic of this function is the same as 'wgf_mem()'.
   // See the comments there for more detail on what is going on.
   trunc_pol_t *w = new_zero_trunc_pol();
   matrix_t *L = new_matrix_L(u);

   matrix_t *powL1 = new_zero_matrix(2*G);
   matrix_t *powL2 = new_zero_matrix(2*G);

   handle_memory_error(w);
   handle_memory_error(L);
   handle_memory_error(powL1);
   handle_memory_error(powL2);

   // Update weighted generating function with
   // one-segment reads (i.e. tail only).
   trunc_pol_update_add(w, L->term[2*G-1]);

   matrix_mult(powL1, L, L);

   // Update weighted generating function with two-segment reads.
   trunc_pol_update_add(w, powL1->term[2*G-1]);

   // Every non-tail segment contains a terminator that is a mismatch for
   // at least one of the two sequences. The probabilities of occurrence
   // of these symbols are 'b', 'c' and 'd' (below), so the probability
   // of occurrence of a terminator is bounded by 'prob' = max(b,c,d).
   // With an argument similar to that expalained in the function
   // 'wgf_mem()', we can bound the probabibility that a
   // read contains more than m segments with the Binomial distribution.
   // https://en.wikipedia.org/wiki/Binomial_distribution#Tail_Bounds
   const double b = P * u/3.0;
   const double c = (1-P) * u;
   const double d = P * (1-u/3.0);
   const double prob = b > c ? (b > d ? b : d) : (c > d ? c : d);
   for (int n = 2 ; n <= K+1 ; n += 2) {
      // Increase the number of segments and update
      // the weighted generating function accordingly.
      matrix_mult(powL2, L, powL1);
      trunc_pol_update_add(w, powL2->term[2*G-1]);
      matrix_mult(powL1, L, powL2);
      trunc_pol_update_add(w, powL1->term[2*G-1]);
      // In max precision mode, get all possible digits.
      // Otherwise, stop when reaching 1% precision.
      double x = floor(n-1) / ((double) K);
      double bound_on_imprecision = exp(-HH(x, prob)*K);
      if (bound_on_imprecision / w->coeff[K] < EPSILON) break;
   }

   // Clean temporary variables.
   destroy_mat(powL1);
   destroy_mat(powL2);
   destroy_mat(L);

   return w;

in_case_of_failure:
   // Clean everything.
   destroy_mat(powL1);
   destroy_mat(powL2);
   destroy_mat(L);
   free(w);
   return NULL;

}


trunc_pol_t *
wgf_skip
(
   const size_t n     // Skipping.
)
// SYNOPSIS:
//   Compute the probabilities that reads do not contain an on-target
//   skip-n seed of size gamma for the specified static and dynamic
//   parameters.
//
// RETURN:
//   A pointer to a struct of type 'trunc_pol_t' containing the
//   probabilitie of interest, or NULL in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized or if 'malloc()'
//   fails. Initialization is checked indirectly through the call to
//   'new_zero_trunc_pol()'.
{

   // The logic of this function is the same as 'wgf_mem()'.
   // See the comments there for more detail on what is going on.
   trunc_pol_t *w = new_zero_trunc_pol();
   matrix_t *S = new_matrix_S(n);

   matrix_t *powS1 = new_zero_matrix(n+2);
   matrix_t *powS2 = new_zero_matrix(n+2);

   handle_memory_error(w);
   handle_memory_error(S);
   handle_memory_error(powS1);
   handle_memory_error(powS2);

   // Update weighted generating function with
   // one-segment reads (i.e. tail only).
   trunc_pol_update_add(w, S->term[n+1]);

   matrix_mult(powS1, S, S);

   // Update weighted generating function with two-segment reads.
   trunc_pol_update_add(w, powS1->term[n+1]);

   // Every non-tail segment contains an error. We can bound the
   // probabibility that a read contains more than m segments with
   // the Binomial distribution.
   // https://en.wikipedia.org/wiki/Binomial_distribution#Tail_Bounds
   for (int s = 2 ; s <= K+1 ; s += 2) {
      // Increase the number of segments and update
      // the weighted generating function accordingly.
      matrix_mult(powS2, S, powS1);
      trunc_pol_update_add(w, powS2->term[n+1]);
      matrix_mult(powS1, S, powS2);
      trunc_pol_update_add(w, powS1->term[n+1]);
      // In max precision mode, get all possible digits.
      // Otherwise, stop when reaching 1% precision.
      double x = floor(s-1) / ((double) K);
      double bound_on_imprecision = exp(-HH(x, P)*K);
      if (bound_on_imprecision / w->coeff[K] < EPSILON) break;
   }

   // Clean temporary variables.
   destroy_mat(powS1);
   destroy_mat(powS2);
   destroy_mat(S);

   return w;

in_case_of_failure:
   // Clean everything.
   destroy_mat(powS1);
   destroy_mat(powS2);
   destroy_mat(S);
   free(S);
   return NULL;

}


trunc_pol_t *
wgf_skip_dual
(
   const size_t n,    // Skipping.
   const double u     // Divergence rate.
)
// SYNOPSIS:
//   Compute the probabilities that reads contain neither an on-target
//   skip-n gamma seed nor an off-target skip-n gamma seed for the
//   specified static and dynamic parameters, where there is N = 1
//   duplicate of the target.
//
// RETURN:
//   A pointer to a struct of type 'trunc_pol_t' containing the
//   probabilitie of interest, or NULL in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized or if 'malloc()'
//   fails. Initialization is checked indirectly through the call to
//   'new_zero_trunc_pol()'.
{

   // Assume parameters were checked by the caller.  See note 4.6.1.
   // about checking dynamic parameters if this code is reused.
   // -- BEGIN BLOCK -- //
   // Check dynamic parameters. Only 'u' and 'N' must be checked, so
   // we pass 'K' as a first parameter, which cannot trigger an error.
   //if (!dynamic_params_OK(K,u,N)) {
   //   goto in_case_of_failure;
   //}
   // -- END BLOCK -- //

   // The logic of this function is the same as 'wgf_mem()'.
   // See the comments there for more detail on what is going on.
   trunc_pol_t *w = new_zero_trunc_pol();
   matrix_t *T = new_matrix_T(n,u);

   matrix_t *powT1 = new_zero_matrix(n+2*G);
   matrix_t *powT2 = new_zero_matrix(n+2*G);

   handle_memory_error(w);
   handle_memory_error(T);
   handle_memory_error(powT1);
   handle_memory_error(powT2);

   // Update weighted generating function with
   // one-segment reads (i.e. tail only).
   trunc_pol_update_add(w, T->term[n+2*G-1]);

   matrix_mult(powT1, T, T);

   // Update weighted generating function with two-segment reads.
   trunc_pol_update_add(w, powT1->term[n+2*G-1]);

   // https://en.wikipedia.org/wiki/Binomial_distribution#Tail_Bounds
   const double b = P * u/3.0;
   const double c = (1-P) * u;
   const double d = P * (1-u/3.0);
   const double prob = b > c ? (b > d ? b : d) : (c > d ? c : d);
   for (int s = 2 ; s <= K+1 ; s += 2) {
      // Increase the number of segments and update
      // the weighted generating function accordingly.
      matrix_mult(powT2, T, powT1);
      trunc_pol_update_add(w, powT2->term[n+2*G-1]);
      matrix_mult(powT1, T, powT2);
      trunc_pol_update_add(w, powT1->term[n+2*G-1]);
      // In max precision mode, get all possible digits.
      // Otherwise, stop when reaching 1% precision.
      double x = floor(s/2) / ((double) K);
      double bound_on_imprecision = exp(-HH(x, prob)*K);
      if (bound_on_imprecision / w->coeff[K] < EPSILON) break;
   }

   // Clean temporary variables.
   destroy_mat(powT1);
   destroy_mat(powT2);
   destroy_mat(T);

   return w;

in_case_of_failure:
   // Clean everything.
   destroy_mat(powT1);
   destroy_mat(powT2);
   destroy_mat(T);
   free(w);
   return NULL;

}


// SECTION 4.10. LOW-LEVEL MONTE CARLO MARKOV CHAIN FUNCTIONS //

size_t
rgeom
(
   const double prob
)
{
   return ceil( log(runifMT()) / log(1-prob) );
}


size_t
rpos
(
   const size_t m,    // The number of hard masks.
   const size_t i     // Segment size (max value).
)
// SYNOPSIS:
//   Simulate the random position where the last hard mask fails.
//
// RETURN:
//   A number greater than or equal to 1 indicating the position where
//   the last hard mask failed.
{
   if (m == 0) return 1;

   // The principle is to invert the probability distribution. The
   // probability that all m hard masks have failed at position j is
   // (1-(1-u)^j)^m. This is the cumulative distribution of the random
   // variable j. The probability that they have all failed at position
   // j, given that they have all failed at position i >= j is the ratio
   // (1-(1-u)^j)^m / (1-(1-u)^i)^m, which is now the cumulative
   // distribution of the random variable j, given that it is less than
   // or equal to i.
   //
   // If X is the m-th root of a U(0,1) random variable, the event
   //
   //          j  <  log(1-(1-(1-u)^i)*X) / log(1-u)  <  j+1 
   //
   // is equivalent to the event
   //
   //  (1-(1-u)^j) / (1-(1-u)^i) <  X  < (1-(1-u)^(j+1)) / (1-(1-u)^i),
   //
   // which corresponds to the probability distribution of j given that
   // it is less than i. The formula below directly generates a random
   // variable with the requested probability distribution.
   const double mth_root_of_unif = pow(runifMT(), 1.0/m);
   return ceil( log(1-XI[i]*mth_root_of_unif) / log(XIc[1]) );
}


void
one_mcmc_mem
(
   const int      N,
         double * pos
)
// SYNOPSIS:
//   Simulate one read with MEM seeds using Monte Carlo Markov chains
//   with the specified static parameters. Update the vector 'pos' with
//   the positions of the read where a MEM seed is present, i.e., the
//   positions such that a MEM seed would be present if the read
//   terminated at this position.
{

   size_t i;     // Size of the error-free segment.
   size_t m = 0; // Number of hard masks.

   // Note: stop if last nucleotide is an error.
   for (int sz = K ; sz > 0 ; sz -= (i+1)) {

      // Get size of the error-free segment.
      i = rgeom(P) - 1;

      if (i >= sz) {
         // If the error-free segment is longer than the read, resize
         // the error-free segment to 'sz' and ignore soft masks.
         
         // If the size of the read is shorter than the minimum
         // seed size 'g', then there is no good seed.
         int long_enough = sz >= G;
         // Otherwise, there is a good seed if there is
         // no survivng hard mask.
         int unmasked = rbinom(m, XIc[sz]) == 0;
         if (unmasked && long_enough) {
            // The segment is unmasked and long enough: there is a global
            // MEM seed. We need to find its position.
            size_t failpos = rpos(m, sz);
            size_t from = failpos < G ? K-sz + G : K-sz + failpos;
            for (int j = from ; j <= K ; j++) pos[j]++;
         }
         return;
      }
      else {
         // The next error is before the end of the read.
         // Get the number of hard masks that survive i
         // nucleotides -- each survives with probability (1-u)^i.
         size_t hsurv = rbinom(m, XIc[i]);

         // Get the number of soft masks that survive i+1
         // nucleotides -- each survives with probability (1-u)^i*u/3.
         size_t ssurv = rbinom(N-m, ETAc[i]);

         // If the segment has size greater than the
         // minimum seed length and there are no surviving
         // threads then there is a good seed.
         int unmasked = hsurv == 0;
         int long_enough = i >= G;
         if (unmasked && long_enough) {
            int all_soft_masks_failed = ssurv == 0;
            size_t failpos = rpos(m, i);
            size_t from = failpos < G ? K-sz + G : K-sz + failpos;
            if (all_soft_masks_failed) {
               for (int j = from ; j <= K ; j++) pos[j]++;
               return;
            }
            else {
               size_t to = K-sz+i;
               for (int j = from ; j <= to ; j++) pos[j]++;
            }
         }

         // Otherwise, get the number of strictly masking threads (m).
         // All the surviving soft masks become strictly masking. The
         // hard masks have a probability u/3 of matching the error.
         // The soft masks that did not survive have a probability
         // pp = u/3 * xi / eta of matching the error. 
         double pp = 1 - ETA[0] / ETA[i];
         // The new state is down/m.
         m = ssurv + rbinom(m, ETAc[0]) + rbinom(N-m-ssurv, pp);

      }

   }

}


void
one_mcmc_skip
(
   const size_t   n,    // Skipping.
   const double   u,    // Divergence.
         double * pos
)
// SYNOPSIS:
//   Simulate one read with skip seeds using Monte Carlo Markov chains
//   with the specified static parameters. Update the vector 'pos' with
//   the positions of the read where a skip seed is present, i.e., the
//   positions such that a skip seed would be present if the read
//   terminated at this position.
//
//   This function is not actuall used in the library because it is
//   quite slow. It is here only for testing purpose, to make sure
//   that the computations with weighted generating functions are
//   correct.
{

   const double a = (1-P)*(1-u);
   const double b = (1-P)*u;
   const double c = P*u/3;

   int spls = 0;
   int smin = 0;

   size_t i;     // Size of the double-match segment.

   // Note: stop if last nucleotide is an error.
   for (int sz = 0 ; sz < K ; ) {
      // Get size of the error-free segment.
      i = rgeom(1-a) - 1;
      if (i >= K-sz) {
         spls += K-sz;
         smin += K-sz;
         break;
      }
      spls += i;
      smin += i;
      if (spls >= G || smin >= G) {
         int pospls = spls >= G ? sz+i - (spls-G) : K+1;
         int posmin = smin >= G ? sz+i - (smin-G) : K+1;
         int from = pospls < posmin ? pospls : posmin;
         for (int j = from ; j <= K ; j++) pos[j]++;
         return;
      }
      sz += i+1;
      double r = runifMT() * (1-a);
      if (r < b) {
         spls++;
         smin = -modulo(-sz,n+1);
      }
      else if (r < b+c) {
         spls = -modulo(-sz,n+1);
         smin++;
      }
      else {
         spls = smin = -modulo(-sz,n+1);
      }
   }
   if (spls >= G || smin >= G) {
      int pospls = spls >= G ? K - (spls-G) : K+1;
      int posmin = smin >= G ? K - (smin-G) : K+1;
      int from = pospls < posmin ? pospls : posmin;
      for (int j = from ; j <= K ; j++) pos[j]++;
   }
   return;
}


// SECTION 4.11. HIGH-LEVEL MONTE CARLO MARKOV CHAIN FUNCTIONS //

trunc_pol_t *
mem_mcmc
(
   const double u, // Divergence rate.
   const int N     // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probabilities that reads do not contain an on-target
//   MEM seed for specified static and dynamic parameters, using the
//   Monte Carlo Markov chain approach.
//
// RETURN:
//   A pointer to a struct of type 'trunc_pol_t' containing the
//   probabilitie of interest, or NULL in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized or if 'malloc()'
//   fains. Initialization is checked indirectly through the call to
//   'new_zero_trunc_pol()'.
{

   // Assume parameters were checked by the caller.  See note 4.6.1.
   // about checking dynamic parameters if this code is reused.
   // -- BEGIN BLOCK -- //
   // Check dynamic parameters. Only 'u' and 'N' must be checked, so
   // we pass 'K' as a first parameter, which cannot trigger an error.
   //if (!dynamic_params_OK(K,u,N)) {
   //   goto in_case_of_failure;
   //}
   // -- END BLOCK -- //

   
   trunc_pol_t *w = new_zero_trunc_pol();
   handle_memory_error(w);

   // Precompute intermediate quantities.
   handle_memory_error(XI = malloc((K+1) * sizeof(double)));
   handle_memory_error(XIc = malloc((K+1) * sizeof(double)));
   handle_memory_error(ETA = malloc((K+1) * sizeof(double)));
   handle_memory_error(ETAc = malloc((K+1) * sizeof(double)));

   for (int i = 0 ; i <= K ; i++) {
      XI[i]   = 1 - pow(1-u,i);
      XIc[i]  = 1 - XI[i];
      ETA[i]  = 1 - pow(1-u,i) * u/3.0;
      ETAc[i] = 1 - ETA[i];
   }

   w->monodeg = K+1;
   for (int i = 0 ; i < MCMC_SAMPLINGS ; i++) {
      one_mcmc_mem(N, w->coeff);
   }

   for (int i = 0 ; i <= K ; i++) {
      w->coeff[i] = 1.0 - w->coeff[i] / MCMC_SAMPLINGS;
   }

   free(XI);    XI   = NULL;
   free(XIc);   XIc  = NULL;
   free(ETA);   ETA  = NULL;
   free(ETAc);  ETAc = NULL;

   return w;

in_case_of_failure:
   free(XI);    XI   = NULL;
   free(XIc);   XIc  = NULL;
   free(ETA);   ETA  = NULL;
   free(ETAc);  ETAc = NULL;
   free(w);
   return NULL;

}


trunc_pol_t *
compute_skipseedp_mcmc
(
   const size_t n,   // Skipping.
   const double u    // Divergence rate.
)
// SYNOPSIS:
//   Compute the probabilities that reads do not contain a skip-n
//   seed for the target and a single given off-target sequence, for
//   specified static and dynamic parameters, using the Monte Carlo
//   Markov chain approach.
//
// RETURN:
//   A pointer to a struct of type 'trunc_pol_t' containing the
//   probabilitie of interest, or NULL in case of failure.
//
// FAILURE:
//   Fails if parameters are unininitialized or if 'malloc()'
//   fains. Initialization is checked indirectly through the call to
//   'new_zero_trunc_pol()'.
{

   // Assume parameters were checked by the caller.  See note 4.6.1.
   // about checking dynamic parameters if this code is reused.
   // -- BEGIN BLOCK -- //
   // Check dynamic parameters. Only 'u' and 'N' must be checked, so
   // we pass 'K' as a first parameter, which cannot trigger an error.
   //if (!dynamic_params_OK(K,u,N)) {
   //   goto in_case_of_failure;
   //}
   // -- END BLOCK -- //

   trunc_pol_t *w = new_zero_trunc_pol();
   handle_memory_error(w);

   w->monodeg = K+1;
   for (int i = 0 ; i < MCMC_SAMPLINGS ; i++) {
      one_mcmc_skip(n, u, w->coeff);
   }

   for (int i = 0 ; i <= K ; i++) {
      w->coeff[i] = 1.0 - w->coeff[i] / MCMC_SAMPLINGS;
   }

   return w;

in_case_of_failure:
   free(w);
   return NULL;

}
   

// SECTION 4.12. HIGH-LEVEL LIBRARY FUNCTIONS //


double *
exact_seed_nullp        // VISIBLE //
(
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that exact seeding is null.
//
// RETURN:
//   A pointer to newly allocated 'trunc_pol_t', or NULL in case
//   of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform, if 'n' is negative or if 'malloc()'
//   fails.
{

   trunc_pol_t *PS0 = NULL;
   trunc_pol_t *PS0S1 = NULL;
   
   // Check dynamic parameters 'u' and 'N'.
   if (!dynamic_params_OK(K,u,N)) {
      goto in_case_of_failure;
   }

   PS0   = wgf_seed();     // Prob no exact seed.
   PS0S1 = wgf_dual(u);    // Prob no hit (two seq model).

   if (PS0 == NULL || PS0S1 == NULL)
      goto in_case_of_failure;

   // This is the probability that there is no exact seed
   // (on-target or off-target) given that there is no on-target
   // exact seed, multiplied by the probability that there is no
   // on-target exact seed. See equation (8).
   for (int i = 0 ; i <= K ; i++) {
      // Store everything in 'PS0'.
      PS0->coeff[i] = 
         PS0->coeff[i] * pow(PS0S1->coeff[i] / PS0->coeff[i], N);
      // Suppress fluctuations in very small numbers.
      if (PS0->coeff[i] < 0) PS0->coeff[i] = 0.0;
   }

   free(PS0S1);

   return to_array_of_double(PS0);

in_case_of_failure:
   free(PS0);
   free(PS0S1);
   return NULL;

}


double *
exact_seed_offp        // VISIBLE //
(
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that exact seeding is off target.
//
// RETURN:
//   A pointer to newly allocated 'trunc_pol_t', or NULL in case
//   of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform or if 'malloc()' fails.
{

   trunc_pol_t *PS0 = NULL;
   trunc_pol_t *PS0S1 = NULL;
   
   // Check dynamic parameters 'u' and 'N'.
   if (!dynamic_params_OK(K,u,N)) {
      goto in_case_of_failure;
   }

   PS0   = wgf_seed();     // Prob no exact seed.
   PS0S1 = wgf_dual(u);    // Prob no hit (two seq model).

   if (PS0 == NULL || PS0S1 == NULL)
      goto in_case_of_failure;

   // This is the probability that there is no on-target exact
   // seed minus the probability that there is no positive at
   // all. The latter is computed as the probability that
   // there is no exact seed (on-target or off-target) given
   // that there is no on-target exact seed, multiplied by the
   // probability that there is no on-target exact seed. See
   // equation (9).
   for (int i = 0 ; i <= K ; i++) {
      // Store everything in 'PS0'.
      PS0->coeff[i] = PS0->coeff[i] -
         PS0->coeff[i] * pow(PS0S1->coeff[i] / PS0->coeff[i], N);
      // Suppress fluctuations in very small numbers.
      if (PS0->coeff[i] < 0) PS0->coeff[i] = 0.0;
   }

   free(PS0S1);

   return to_array_of_double(PS0);

in_case_of_failure:
   free(PS0);
   free(PS0S1);
   return NULL;

}


double *
skip_seed_nullp        // VISIBLE //
(
   const int    n,   // Skipping.
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that exact seeding is off null.
//
// RETURN:
//   A pointer to newly allocated 'trunc_pol_t', or NULL in case
//   of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform or if 'malloc()' fails.
{

   trunc_pol_t *PS0 = NULL;
   trunc_pol_t *PS0S1 = NULL;
   
   // Check dynamic parameters 'u' and 'N'.
   if (!dynamic_params_OK(K,u,N)) {
      goto in_case_of_failure;
   }

   if (n < 0) {
      warning("skipping (n) must be non-negative",
            __func__, __LINE__);
      goto in_case_of_failure;
   }

   PS0   = wgf_skip(n);          // Prob no exact seed.
   PS0S1 = wgf_skip_dual(n, u);  // Prob no hit (two seq model).

   if (PS0 == NULL || PS0S1 == NULL)
      goto in_case_of_failure;

   // This is the probability that there is no skip seed
   // (on-target or off-target) given that there is no on-target
   // skip seed, multiplied by the probability that there is no
   // on-target skip seed. See equation (8).
   for (int i = 0 ; i <= K ; i++) {
      // Store everything in 'PS0'.
      PS0->coeff[i] = 
         PS0->coeff[i] * pow(PS0S1->coeff[i] / PS0->coeff[i], N);
      // Suppress fluctuations in very small numbers.
      if (PS0->coeff[i] < 0) PS0->coeff[i] = 0.0;
   }

   free(PS0S1);

   return to_array_of_double(PS0);

in_case_of_failure:
   free(PS0);
   free(PS0S1);
   return NULL;

}


double *
skip_seed_offp        // VISIBLE //
(
   const int    n,   // Skipping.
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that exact seeding is off target.
//
// RETURN:
//   A pointer to newly allocated 'trunc_pol_t', or NULL in case
//   of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform, if 'n' is negative or if 'malloc()'
//   fails.
{

   trunc_pol_t *PS0 = NULL;
   trunc_pol_t *PS0S1 = NULL;
   
   // Check dynamic parameters 'u' and 'N'.
   if (!dynamic_params_OK(K,u,N)) {
      goto in_case_of_failure;
   }

   if (n < 0) {
      warning("skipping (n) must be non-negative",
            __func__, __LINE__);
      goto in_case_of_failure;
   }

   PS0   = wgf_skip(n);          // Prob no exact seed.
   PS0S1 = wgf_skip_dual(n, u);  // Prob no hit (two seq model).

   if (PS0 == NULL || PS0S1 == NULL)
      goto in_case_of_failure;

   // This is the probability that there is no on-target skip
   // seed minus the probability that there is no positive at
   // all. The latter is computed as the probability that
   // there is no skip seed (on-target or off-target) given
   // that there is no on-target skip seed, multiplied by the
   // probability that there is no on-target skip seed. See
   // equation (9).
   for (int i = 0 ; i <= K ; i++) {
      // Store everything in 'PS0'.
      PS0->coeff[i] = PS0->coeff[i] -
         PS0->coeff[i] * pow(PS0S1->coeff[i] / PS0->coeff[i], N);
      // Suppress fluctuations in very small numbers.
      if (PS0->coeff[i] < 0) PS0->coeff[i] = 0.0;
   }

   free(PS0S1);

   return to_array_of_double(PS0);

in_case_of_failure:
   free(PS0);
   free(PS0S1);
   return NULL;

}


double *
mem_seed_nullp        // VISIBLE //
(
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that MEM seeding is null.
//
// RETURN:
//   A pointer to newly allocated 'trunc_pol_t', or NULL in case
//   of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform or if 'malloc()' fails.
{

   // This is exactly the same probability.
   return exact_seed_nullp(u, N);

}


double *
mem_seed_offp        // VISIBLE //
(
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that MEM seeding is off target.
//
// RETURN:
//   A pointer to newly allocated 'trunc_pol_t', or NULL in case
//   of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform or if 'malloc()' fails.
{

   trunc_pol_t *PS0 = NULL;
   trunc_pol_t *PS0S1 = NULL;
   trunc_pol_t *PM0 = NULL;
   
   // Check dynamic parameters 'u' and 'N'.
   if (!dynamic_params_OK(K,u,N)) {
      goto in_case_of_failure;
   }

   PS0   = wgf_seed();     // Prob no exact seed.
   PS0S1 = wgf_dual(u);    // Prob no hit (two seq model).
   PM0   = wgf_mem(u, N);  // Prob no on-target MEM seed.

   if (PS0 == NULL || PS0S1 == NULL || PM0 == NULL)
      goto in_case_of_failure;

   // This is the probability that there is no on-target MEM
   // seed minus the probability that there is no positive at
   // all. The latter is computed as the probability that
   // there is no exact seed (on-target or off-target) given
   // that there is no on-target exact seed, multiplied by the
   // probability that there is no on-target exact seed. See
   // equation (21).
   for (int i = 0 ; i <= K ; i++) {
      // Store everything in 'M0'.
      PM0->coeff[i] = PM0->coeff[i] -
         PS0->coeff[i] * pow(PS0S1->coeff[i] / PS0->coeff[i], N);
      // Suppress fluctuations in very small numbers.
      if (PM0->coeff[i] < 0) PM0->coeff[i] = 0.0;
   }

   free(PS0);
   free(PS0S1);

   return to_array_of_double(PM0);

in_case_of_failure:
   free(PS0);
   free(PS0S1);
   free(PM0);
   return NULL;

}


double *
mem_seed_offp_mcmc   // VISIBLE //
(
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that MEM seeding is off target, using
//   the MCMC algorithm (for speed).
//
// RETURN:
//   A pointer to newly allocated 'trunc_pol_t', or NULL in case
//   of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform or if 'malloc()' fails.
{

   trunc_pol_t *PS0 = NULL;
   trunc_pol_t *PS0S1 = NULL;
   trunc_pol_t *PM0 = NULL;
   
   // Check dynamic parameters 'u' and 'N'.
   if (!dynamic_params_OK(K,u,N)) {
      goto in_case_of_failure;
   }

   PS0   = wgf_seed();      // Prob no exact seed.
   PS0S1 = wgf_dual(u);     // Prob no hit (two seq model).
   PM0   = mem_mcmc(u, N);  // Prob no on-target MEM seed.

   if (PS0 == NULL || PS0S1 == NULL || PM0 == NULL)
      goto in_case_of_failure;

   // This is the probability that there is no on-target MEM
   // seed minus the probability that there is no positive at
   // all. The latter is computed as the probability that
   // there is no exact seed (on-target or off-target) given
   // that there is no on-target exact seed, multiplied by the
   // probability that there is no on-target exact seed. See
   // equation (21).
   for (int i = 0 ; i <= K ; i++) {
      // Store everything in 'M0'.
      PM0->coeff[i] = PM0->coeff[i] -
         PS0->coeff[i] * pow(PS0S1->coeff[i] / PS0->coeff[i], N);
      // Suppress fluctuations in very small numbers.
      if (PM0->coeff[i] < 0) PM0->coeff[i] = 0.0;
   }

   free(PS0);
   free(PS0S1);

   return to_array_of_double(PM0);

in_case_of_failure:
   free(PS0);
   free(PS0S1);
   free(PM0);
   return NULL;

}


double
auto_exact_seed_nullp   // VISIBLE //
(
   const int    k,   // Segment or read size.
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that the exact seeding process is null.
//   Take care of the detail regarding the caching so that
//   computations are fast.
//
// RETURN:
//   A double-precision number with the probability of interest,
//   or 'nan' in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform or if 'malloc()' fails.
{

   double *prob = NULL;

   // Check dynamic parameters.
   if (!dynamic_params_OK(k,u,N)) {
      goto in_case_of_failure;
   }

   // Squish 'N' for "coarseness".
   size_t sqN = squish(N);

   // Retrieve probabilities from H1N.
   rec_t *record = lookup(H1N, 0, u, sqN);

   // Otherwise compute them.
   if (record == NULL) {
      prob = exact_seed_nullp(u, sqN);
      // Insert in H1N.
      record = insert(H1N, 0, u, sqN, prob);
      if (record == NULL) {
         goto in_case_of_failure;
      }
   }

   return record->prob[k];

in_case_of_failure:
   free(prob);
   return 0.0/0.0;

}


double
auto_exact_seed_offp   // VISIBLE //
(
   const int    k,   // Segment or read size.
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that the exact seeding process is off
//   target. Take care of the detail regarding caching so that
//   computations are fast.
//
// RETURN:
//   A double-precision number with the probability of interest,
//   or 'nan' in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform or if 'malloc()' fails.
{

   double *prob = NULL;

   // Check dynamic parameters.
   if (!dynamic_params_OK(k,u,N)) {
      goto in_case_of_failure;
   }

   // Squish 'N' for "coarseness".
   size_t sqN = squish(N);

   // Retrieve probabilities from H1O.
   rec_t *record = lookup(H1O, 0, u, sqN);

   // Otherwise compute them.
   if (record == NULL) {
      prob = exact_seed_offp(u, sqN);
      // Insert in H1O.
      record = insert(H1O, 0, u, sqN, prob);
      if (record == NULL) {
         goto in_case_of_failure;
      }
   }

   return record->prob[k];

in_case_of_failure:
   free(prob);
   return 0.0/0.0;

}


double
auto_skip_seed_nullp   // VISIBLE //
(
   const int    k,   // Segment or read size.
   const int    n,   // Skipping.
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that the skip seeding process is null.
//   Take care of the detail regarding caching so that computations
//   are fast.
//
// RETURN:
//   A double-precision number with the probability of interest,
//   or 'nan' in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform or if 'malloc()' fails.
{

   double *prob = NULL;

   // Check dynamic parameters.
   if (!dynamic_params_OK(k,u,N)) {
      goto in_case_of_failure;
   }

   if (n < 0) {
      warning("skipping (n) must be non-negative",
            __func__, __LINE__);
      goto in_case_of_failure;
   }

   // Squish 'N' for "coarseness".
   size_t sqN = squish(N);

   // Retrieve probabilities from H2N.
   rec_t *record = lookup(H2N, n, u, sqN);

   // Otherwise compute them.
   if (record == NULL) {
      prob = skip_seed_nullp(n, u, sqN);
      // Insert in H2N.
      record = insert(H2N, n, u, sqN, prob);
      if (record == NULL) {
         goto in_case_of_failure;
      }
   }

   return record->prob[k];

in_case_of_failure:
   free(prob);
   return 0.0/0.0;

}


double
auto_skip_seed_offp   // VISIBLE //
(
   const int    k,   // Segment or read size.
   const int    n,   // Skipping.
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that the skip seeding process is off
//   target. Take care of the detail regarding caching so that
//   computations are fast.
//
// RETURN:
//   A double-precision number with the probability of interest,
//   or 'nan' in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform or if 'malloc()' fails.
{

   double *prob = NULL;

   // Check dynamic parameters.
   if (!dynamic_params_OK(k,u,N)) {
      goto in_case_of_failure;
   }

   if (n < 0) {
      warning("skipping (n) must be non-negative",
            __func__, __LINE__);
      goto in_case_of_failure;
   }

   // Squish 'N' for "coarseness".
   size_t sqN = squish(N);

   // Retrieve probabilities from H2O.
   rec_t *record = lookup(H2O, n, u, sqN);

   // Otherwise compute them.
   if (record == NULL) {
      prob = skip_seed_offp(n, u, sqN);
      // Insert in H2O.
      record = insert(H2O, n, u, sqN, prob);
      if (record == NULL) {
         goto in_case_of_failure;
      }
   }

   return record->prob[k];

in_case_of_failure:
   free(prob);
   return 0.0/0.0;

}


double
auto_mem_seed_nullp   // VISIBLE //
(
   const int    k,   // Segment or read size.
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that the MEM seeding process is null.
//   Take care of the detail regarding caching so that computations
//   are fast.
//
// RETURN:
//   A double-precision number with the probability of interest,
//   or 'nan' in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform or if 'malloc()' fails.
{

   double *prob = NULL;

   // Check dynamic parameters.
   if (!dynamic_params_OK(k,u,N)) {
      goto in_case_of_failure;
   }

   // Squish 'N' for "coarseness".
   size_t sqN = squish(N);

   // Retrieve probabilities from H3N.
   rec_t *record = lookup(H3N, 0, u, sqN);

   // Otherwise compute them.
   if (record == NULL) {
      prob = mem_seed_nullp(u, sqN);
      // Insert in H3N.
      record = insert(H3N, 0, u, sqN, prob);
      if (record == NULL) {
         goto in_case_of_failure;
      }

   }

   return record->prob[k];

in_case_of_failure:
   free(prob);
   return 0.0/0.0;

}


double
auto_mem_seed_offp   // VISIBLE //
(
   const int    k,   // Segment or read size.
   const double u,   // Divergence rate.
   const int    N    // Number of duplicates.
)
// SYNOPSIS:
//   Compute the probability that the MEM seeding process is off
//   target. Take care of the detail regarding the choice of the
//   algorithm and the caching so that computations are fast.
//
// RETURN:
//   A double-precision number with the probability of interest,
//   or 'nan' in case of failure.
//
// FAILURE:
//   Fails if static parameters are unininitialized, if dynamic
//   parameters are not conform or if 'malloc()' fails.
{

   double *prob = NULL;

   // Check dynamic parameters.
   if (!dynamic_params_OK(k,u,N)) {
      goto in_case_of_failure;
   }

   // Squish 'N' for "coarseness".
   size_t sqN = squish(N);

   // Retrieve probabilities from H3O.
   rec_t *record = lookup(H3O, 0, u, sqN);

   // Otherwise compute them.
   if (record == NULL) {
      // If 'N' > 20 and 'P' < 0.05 use MCMC.
      if (N < 21 && P < .05) {
         prob = mem_seed_offp(u, sqN);
      }
      else {
         prob = mem_seed_offp_mcmc(u, sqN);
      }

      // Insert in H3O.
      record = insert(H3O, 0, u, sqN, prob);
      if (record == NULL) {
         goto in_case_of_failure;
      }

   }

   return record->prob[k];

in_case_of_failure:
   free(prob);
   return 0.0/0.0;

}


void clean_prob_storage (void) { clean_hash(Y1); } // VISIBLE //

void
dump_prob_to_file   // VISIBLE //
(
    FILE   * f
)
// SYNOPSIS:
{

   // Write the static parameters.
   fprintf(f, "# gamma:%d, k:%d, p:%.3f\n", G, K, P);

   // Iterate over the values of the hash table.
   for (int i = 0 ; i < HSIZE ; i++) {
      for (rec_t *rec = Y1[i] ; rec != NULL ; rec = rec->next) {
         fprintf(f, "%d\t%.2f\t%d", rec->n, rec->u / 100.0, rec->N);
         for (int j = 0 ; j <= K ; j++) {
            fprintf(f, "\t%.8f", rec->prob[j]);
         }
         fprintf(f, "\n");
      }
   }

}


int
load_prob_from_file   // VISIBLE //
(
    FILE   * f
)
// SYNOPSIS:
{

   // For maximum portability, we are not using any external code
   // to read the file. Unfortunately, this makes it less readable.
   char line[4096] = {0};
   char *s = NULL;
   
   // Read the header.
   if ((s = fgets(line, 4096, f)) == NULL) {
      warning("input file is empty", __func__, __LINE__);
      goto in_case_of_failure;
   }

   int g = -1;
   int k = -1;
   double p = -1.;

   // Decode the header.
   if (sscanf(line, "# gamma:%d, k:%d, p:%lf", &g, &k, &p) != 3) {
      warning("could not read header", __func__, __LINE__);
      goto in_case_of_failure;
   }

   sesame_set_static_params(g, k, p);
   fprintf(stderr, "sesame static parameters: %d, %d, %.3f\n", g, k, p);

   while((s = fgets(line, 4096, f)) != NULL) {

      int n = -1;
      double u = -1.;
      int N = -1;
      double *prob = malloc((K+1) * sizeof(double));
      handle_memory_error(prob);

      int counter = 0;
      char *token = strtok(line, "\t");
      while (token != NULL) {
         if      (counter == 0) n = strtol(token, NULL, 10);
         else if (counter == 1) u = strtod(token, NULL);
         else if (counter == 2) N = strtol(token, NULL, 10);
         else                   prob[counter-3] = strtod(token, NULL);
         counter++;
         token = strtok(NULL, "\t");
      }

      if (counter != K+4) {
         warning("could not parse input file", __func__, __LINE__);
         goto in_case_of_failure;
      }

      if (!store_prob(n, u, N, prob)) {
         warning("memory error", __func__, __LINE__);
         goto in_case_of_failure;
      }

   }

   return SUCCESS;

in_case_of_failure:
   clean_prob_storage();
   return FAILURE;

}


#if 0
double
average_errors
(
   const int g,
   const int k,
   const double p
)
// Use recurrence.
{

   if (k <= g) return p*k;

   // Allocate at least '2*g+2' numbers for convenience.
   size_t sz = (k+1 < 2*g+2 ? 2*g+2 : k+1) * sizeof(double);
   double *s = malloc(sz);
   handle_memory_error(s);

   const double two_pq_pow_gamma = 2 * p * pow(1-p,g);
   const double pq_pow_gamma_square = pow(p * pow(1-p,g),2);

   for (int i = 0 ; i <= g ; i++) s[i] = p*i;

   s[g+1] = p*(g+1) - two_pq_pow_gamma;

   for (int i = g+2 ; i < 2*g+1 ; i++) {
      s[i] = 2*s[i-1] - s[i-2] -
         two_pq_pow_gamma * (s[i-g-1] - s[i-g-2]);
   }

   s[2*g+1] = 2*s[2*g] - s[2*g-1] -
      two_pq_pow_gamma * (s[g] - s[g-1]) + p*pow(1-p,2*g);

   for (int i = 2*g+2 ; i <= k ; i++) {
      s[i] = 2*s[i-1] - s[i-2] -
         two_pq_pow_gamma * (s[i-g-1] - s[i-g-2]) -
         pq_pow_gamma_square * s[i-2*g-2];
   }

   double value = s[k];
   free(s);

   return value;

in_case_of_failure:
   return 0.0 / 0.0; // nan

}

double
special_average
(
   const int g,
   const int k,
   const double p
)
// TODO: explain the logic of this in the LaTex document.
{

   const double num = p*k - average_errors(g,k,p);
   const double denom = 1.0 - exact_seed_prob(g,k,0) - pow(1-p,k);
   return num / denom; 

}
#endif


// http://www.math.keio.ac.jp/~matumoto/ver980409.html

// This is the ``Mersenne Twister'' random number generator MT19937,
// which generates pseudorandom integers uniformly distributed in
// 0..(2^32 - 1) starting from any odd seed in 0..(2^32 - 1).  This
// version is a recode by Shawn Cokus (Cokus@math.washington.edu) on
// March 8, 1998 of a version by Takuji Nishimura (who had suggestions
// from Topher Cooper and Marc Rieffel in July-August 1997).
//
// Effectiveness of the recoding (on Goedel2.math.washington.edu, a DEC
// Alpha running OSF/1) using GCC -O3 as a compiler: before recoding:
// 51.6 sec. to generate 300 million random numbers; after recoding: 24.0
// sec. for the same (i.e., 46.5% of original time), so speed is now
// about 12.5 million random number generations per second on this
// machine.
//
// According to the URL <http://www.math.keio.ac.jp/~matumoto/emt.html>
// (and paraphrasing a bit in places), the Mersenne Twister is ``designed
// with consideration of the flaws of various existing generators,'' has
// a period of 2^19937 - 1, gives a sequence that is 623-dimensionally
// equidistributed, and ``has passed many stringent tests, including the
// die-hard test of G. Marsaglia and the load test of P. Hellekalek and
// S.  Wegenkittl.''  It is efficient in memory usage (typically using
// 2506 to 5012 bytes of static data, depending on data type sizes, and
// the code is quite short as well).  It generates random numbers in
// batches of 624 at a time, so the caching and pipelining of modern
// systems is exploited.  It is also divide- and mod-free.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as
// published by the Free Software Foundation (either version 2 of the
// License or, at your option, any later version).  This library is
// distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY, without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
// License for more details.  You should have received a copy of the GNU
// Library General Public License along with this library; if not, write
// to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
// Boston, MA 02111-1307, USA.
//
// The code as Shawn received it included the following notice:
//
//   Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.  When you
//   use this, send an e-mail to <matumoto@math.keio.ac.jp> with an
//   appropriate reference to your work.
//
// It would be nice to CC: <Cokus@math.washington.edu> when you write.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//
// uint32 must be an unsigned integer type capable of holding at least
// 32 bits; exactly 32 should be fastest, but 64 is better on an Alpha
// with GCC at -O3 optimization so try your options and see what's best
// for you
//

typedef unsigned long uint32;

// NOTE: The macros 'N', 'M' and 'K' were renamed 'NN', 'MM' and 'KK
// to avoid possible confusions (and they are undefined below).

#define NN             (624)                // length of state vector
#define MM             (397)                // a period parameter
#define KK             (0x9908B0DFU)        // a magic constant
#define hiBit(u)       ((u) & 0x80000000U)  // mask all but highest bit
#define loBit(u)       ((u) & 0x00000001U)  // mask all but lowest bit
#define loBits(u)      ((u) & 0x7FFFFFFFU)  // mask the highest bit
#define mixBits(u, v)  (hiBit(u)|loBits(v)) // move hi bit u to hi bit v

static uint32   state[NN+1]; // state vector + 1 extra to respect ANSI C
static uint32   *next;      // next random value is computed from here
static int      left = -1;  // can *next++ so many times before reloading


void
seedMT
(
   uint32 seed
)
//
// We initialize state[0..(NN-1)] via the generator
//
//   x_new = (69069 * x_old) mod 2^32
//
// from Line 15 of Table 1, p. 106, Sec. 3.3.4 of Knuth's
// _The Art of Computer Programming_, Volume 2, 3rd ed.
//
// Notes (SJC): I do not know what the initial state requirements
// of the Mersenne Twister are, but it seems this seeding generator
// could be better.  It achieves the maximum period for its modulus
// (2^30) iff x_initial is odd (p. 20-21, Sec. 3.2.1.2, Knuth); if
// x_initial can be even, you have sequences like 0, 0, 0, ...;
// 2^31, 2^31, 2^31, ...; 2^30, 2^30, 2^30, ...; 2^29, 2^29 + 2^31,
// 2^29, 2^29 + 2^31, ..., etc. so I force seed to be odd below.
//
// Even if x_initial is odd, if x_initial is 1 mod 4 then
//
//   the          lowest bit of x is always 1,
//   the  next-to-lowest bit of x is always 0,
//   the 2nd-from-lowest bit of x alternates  ... 0 1 0 1 0 1 0 1 ... ,
//   the 3rd-from-lowest bit of x 4-cycles    ... 0 1 1 0 0 1 1 0 ... ,
//   the 4th-from-lowest bit of x has 8-cycle ... 0 0 0 1 1 1 1 0 ... ,
//    ...
//
// and if x_initial is 3 mod 4 then
//
//   the          lowest bit of x is always 1,
//   the  next-to-lowest bit of x is always 1,
//   the 2nd-from-lowest bit of x alternates  ... 0 1 0 1 0 1 0 1 ... ,
//   the 3rd-from-lowest bit of x 4-cycles    ... 0 0 1 1 0 0 1 1 ... ,
//   the 4th-from-lowest bit of x has 8-cycle ... 0 0 1 1 1 1 0 0 ... ,
//    ...
//
// The generator's potency (min. s>=0 with (69069-1)^s = 0 mod 2^32) is
// 16, which seems to be alright by p. 25, Sec. 3.2.1.3 of Knuth.  It
// also does well in the dimension 2..5 spectral tests, but it could be
// better in dimension 6 (Line 15, Table 1, p. 106, Sec. 3.3.4, Knuth).
//
// Note that the random number user does not see the values generated
// here directly since reloadMT() will always munge them first, so maybe
// none of all of this matters.  In fact, the seed values made here could
// even be extra-special desirable if the Mersenne Twister theory says
// so-- that's why the only change I made is to restrict to odd seeds.
//
{

    register uint32 x = (seed | 1U) & 0xFFFFFFFFU, *s = state;
    register int    j;

    for(left=0, *s++=x, j=NN; --j;
          *s++ = (x*=69069U) & 0xFFFFFFFFU);
 }


uint32 reloadMT(void)
{
   register uint32 *p0=state, *p2=state+2, *pM=state+MM, s0, s1;
   register int    j;

   if(left < -1)
      seedMT(4357U);

   left=NN-1, next=state+1;

   for(s0=state[0], s1=state[1], j=NN-MM+1; --j; s0=s1, s1=*p2++)
      *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? KK : 0U);

   for(pM=state, j=MM; --j; s0=s1, s1=*p2++)
      *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? KK : 0U);

   s1=state[0], *p0 = *pM ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? KK : 0U);
   s1 ^= (s1 >> 11);
   s1 ^= (s1 <<  7) & 0x9D2C5680U;
   s1 ^= (s1 << 15) & 0xEFC60000U;
   return(s1 ^ (s1 >> 18));
}

#undef NN
#undef MM
#undef KK


uint32
randomMT
(void)
{
   uint32 y;

   if(--left < 0)
      return(reloadMT());

   y  = *next++;
   y ^= (y >> 11);
   y ^= (y <<  7) & 0x9D2C5680U;
   y ^= (y << 15) & 0xEFC60000U;
   return(y ^ (y >> 18));
}

// Tiny function to return a uniform pseudo random number between 0
// and 1 using the Mersenne twister algorithm above.

double runifMT (void) { return randomMT() / 4294967295.0; }




/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2014 The R Core Team
 *  Copyright (C) 2007 The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  DESCRIPTION
 *
 *	Random variates from the binomial distribution.
 *
 *  REFERENCE
 *
 *	Kachitvichyanukul, V. and Schmeiser, B. W. (1988).
 *	Binomial random variate generation.
 *	Communications of the ACM 31, 216-222.
 *	(Algorithm BTPEC).
 */


#define repeat for(;;)

size_t
rbinom
(
   size_t n,
   double pp
)
{

   static double c, fm, npq, p1, p2, p3, p4, qn;
   static double xl, xll, xlr, xm, xr;

   static double psave = -1.0;
   static int nsave = -1;
   static int m;

   double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
   double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
   int i, ix, k;

   if (n == 0 || pp == 0.) return 0;
   if (pp == 1.) return n;

   p = pp < .5 ? pp : 1. - pp;
   q = 1. - p;
   np = n * p;
   r = p / q;
   g = r * (n + 1);

// Setup, perform only when parameters change using static (globals).

   if (pp != psave || n != nsave) {
      psave = pp;
      nsave = n;
      if (np < 30.0) {
         /* inverse cdf logic for mean less than 30 */
         qn = pow(q, n);
         goto L_np_small;
      } else {
         ffm = np + p;
         m = (int) ffm;
         fm = m;
         npq = np * q;
         p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
         xm = fm + 0.5;
         xl = xm - p1;
         xr = xm + p1;
         c = 0.134 + 20.5 / (15.3 + fm);
         al = (ffm - xl) / (ffm - xl * p);
         xll = al * (1.0 + 0.5 * al);
         al = (xr - ffm) / (xr * q);
         xlr = al * (1.0 + 0.5 * al);
         p2 = p1 * (1.0 + c + c);
         p3 = p2 + c / xll;
         p4 = p3 + c / xlr;
      }
   } else if (n == nsave) {
      if (np < 30.0)
         goto L_np_small;
   }

   //-------------------------- np = n*p >= 30 : ------------------- //
   repeat {
      u = runifMT() * p4;
      v = runifMT();
      /* triangular region */
      if (u <= p1) {
         ix = (int)(xm - p1 * v + u);
         goto finis;
      }
      /* parallelogram region */
      if (u <= p2) {
         x = xl + (u - p1) / c;
         v = v * c + 1.0 - fabs(xm - x) / p1;
         if (v > 1.0 || v <= 0.)
            continue;
         ix = (int) x;
      } else {
         if (u > p3) {	/* right tail */
            ix = (int)(xr - log(v) / xlr);
            if (ix > n)
               continue;
            v = v * (u - p3) * xlr;
         } else {/* left tail */
            ix = (int)(xl + log(v) / xll);
            if (ix < 0)
               continue;
            v = v * (u - p2) * xll;
         }
      }
      // determine appropriate way to perform accept/reject test //
      k = abs(ix - m);
      if (k <= 20 || k >= npq / 2 - 1) {
         /* explicit evaluation */
         f = 1.0;
         if (m < ix) {
            for (i = m + 1; i <= ix; i++)
               f *= (g / i - r);
         } else if (m != ix) {
            for (i = ix + 1; i <= m; i++)
               f /= (g / i - r);
         }
         if (v <= f)
            goto finis;
      } else {
         // squeezing using upper and lower bounds on log(f(x)) //
         amaxp = (k / npq) *
            ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
         ynorm = -k * k / (2.0 * npq);
         alv = log(v);
         if (alv < ynorm - amaxp)
            goto finis;
         if (alv <= ynorm + amaxp) {
            // stirling's formula to machine accuracy //
            // for the final acceptance/rejection test //
            x1 = ix + 1;
            f1 = fm + 1.0;
            z = n + 1 - fm;
            w = n - ix + 1.0;
            z2 = z * z;
            x2 = x1 * x1;
            f2 = f1 * f1;
            w2 = w * w;
            if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
               goto finis;
         }
      }
   }

L_np_small:
   //---------------------- np = n*p < 30 : ------------------------- //

   repeat {
      ix = 0;
      f = qn;
      u = runifMT();
      repeat {
         if (u < f)
            goto finis;
         if (ix > 110)
            break;
         u -= f;
         ix++;
         f *= (g / ix - r);
      }
   }
finis:
   if (psave > 0.5)
      ix = n - ix;
   return ix;
}
