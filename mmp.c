#define VERSION "1.0"
#define _GNU_SOURCE

#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>

#include "bwt.h"
#include "map.h"
#include "sesame.h"

#define GAMMA 19
#define PROBDEFAULT 0.01
#define SKIPQUALDEFAULT 10
#define QUICK_DUPLICATES 20
#define MAXTHREADSDEFAULT 1

#define min(a,b) (((a) < (b)) ? (a) : (b))

typedef struct uN0_t uN0_t;
typedef struct seedp_t seedp_t;
typedef struct read_t read_t;
typedef struct writerarg_t writerarg_t;
typedef struct batch_t batch_t;

struct uN0_t {
   double u;
   size_t N0;
   double p0;
};

struct read_t {
   char * name;
   char * seq;
   char * phred;
};

struct seedp_t {
   double off;
   double nul;
};

struct writerarg_t {
   batch_t         * first_batch;
   pthread_mutex_t * mutex;
   pthread_cond_t  * cond_writer;
   pthread_cond_t  * cond_reader;
};

enum status_t {idle = 0, map = 1, output = 2, die = 3};

struct batch_t {
   // Data
   index_t           idx;
   size_t            lineid;
   wstack_t        * reads;
   wstack_t        * output;
   // Thread signaling
   enum status_t     status;
   int             * act_threads;
   pthread_mutex_t * mutex;
   pthread_mutex_t * mutex_sesame;
   pthread_cond_t  * cond_reader;
   pthread_cond_t  * cond_writer;
   pthread_cond_t  * cond_reads;
   pthread_cond_t  * cond_threads;
   batch_t         * next_batch;
};

enum fmt_t {unset = 0, fasta = 1, fastq = 2};

static double PROB = PROBDEFAULT;
static double SKIPQUAL = SKIPQUALDEFAULT;
static int    MAXTHREADS = MAXTHREADSDEFAULT;

char* HELP_MSG =
  "Usage:\n"
  "   index:  mmp  --index  index.fasta\n"
  "   map:    mmp [-e 0.01 | -t 1] index.fasta reads.fasta\n"
  "\n"
  "Options:\n"
  "  -e: sequencing error rate (default: 0.01)\n"
  "  -t: number of threads (default: 1)\n"
  "\n";


double digamma(double);
double trigamma(double);


void
build_index
(
   const char * fname
)
{
   // Open fasta file.
   FILE * fasta = fzopen(fname, "r");
   if (fasta == NULL) exit_cannot_open(fname);

   // Aux variables for file writing.
   char * data;
   ssize_t ws;
   size_t sz;

   char * fn = malloc(strlen(fname)+10);
   char * fn2 = malloc(strlen(fname)+10);

   // Read and pack FASTA
   fprintf(stderr, "packing sequence... ");
   pack_fasta(fname);
   fprintf(stderr, "done.\n");

   // Compute bwt using bwt_gen
   fprintf(stderr, "computing bwt...\n");
   sprintf(fn, "%s.dna", fname);
   sprintf(fn2, "%s.bwt", fname);
   size_t bwt_blocksize = 10000000;
   bwt_bwtgen2(fn, fn2, bwt_blocksize);
   fprintf(stderr, "done.\n");

   // Compute Occ table from BWT
   fprintf(stderr, "computing occ from bwt...\n");
   bwt2occ(fname);
   fprintf(stderr, "\rdone.     \n");

   // Compute sampled SA using backward search
   fprintf(stderr, "computing sampled SA with bw search...\n");
   bwt2sa(fname);
   fprintf(stderr, "\rdone.     \n");

   // Load OCC
   sprintf(fn, "%s.occ", fname);
   int focc = open(fn, O_RDONLY);
   if (focc < 0) exit_cannot_open(fn);

   ssize_t mmsz = lseek(focc, 0, SEEK_END);
   occ_t *occ = (occ_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, focc, 0);
   exit_error(occ == NULL);
   close(focc);

   // Create LUT
   fprintf(stderr, "filling lookup table... ");
   lut_t * lut = malloc(sizeof(lut_t));
   fill_lut(lut, occ, (range_t) {.bot=1, .top=occ->txtlen}, 0, 0);
   fprintf(stderr, "done.\n");

   // Unmap occ
   munmap(lut, mmsz);

   // Write the lookup table
   sprintf(fn, "%s.lut", fname);
   int flut = creat(fn, 0644);
   if (flut < 0) exit_cannot_open(fn);

   ws = 0;
   sz = sizeof(lut_t);
   data = (char *) lut;
   while (ws < sz) ws += write(flut, data + ws, sz - ws);
   close(flut);

   // Free memory
   free(lut);
}


index_t
load_index
(
   const char  * fname
)
{

   size_t mmsz;
   char buff[256];

   chr_t * chr = index_load_chr(fname);
   exit_error(chr == NULL);

   sprintf(buff, "%s.dna", fname);
   int fdna = open(buff, O_RDONLY);
   if (fdna < 0) exit_cannot_open(buff);

   mmsz = lseek(fdna, 0, SEEK_END);
   char *dna = (char *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, fdna, 0);
   exit_error(dna == NULL);
   close(fdna);


   sprintf(buff, "%s.sa", fname);
   int fsar = open(buff, O_RDONLY);
   if (fsar < 0) exit_cannot_open(buff);

   mmsz = lseek(fsar, 0, SEEK_END);
   csa_t *csa = (csa_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, fsar, 0);
   exit_error(csa == NULL);
   close(fsar);


   sprintf(buff, "%s.bwt", fname);
   int fbwt = open(buff, O_RDONLY);
   if (fbwt < 0) exit_cannot_open(buff);

   mmsz = lseek(fbwt, 0, SEEK_END);
   bwt_t *bwt = (bwt_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, fbwt, 0);
   exit_error(bwt == NULL);
   close(fbwt);


   sprintf(buff, "%s.occ", fname);
   int focc = open(buff, O_RDONLY);
   if (focc < 0) exit_cannot_open(buff);

   mmsz = lseek(focc, 0, SEEK_END);
   occ_t *occ = (occ_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, focc, 0);
   exit_error(occ == NULL);
   close(focc);


   sprintf(buff, "%s.lut", fname);
   int flut = open(buff, O_RDONLY);
   if (flut < 0) exit_cannot_open(buff);

   mmsz = lseek(flut, 0, SEEK_END);
   lut_t *lut = (lut_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, focc, 0);
   exit_error(lut == NULL);
   close(flut);

   return (index_t) {.chr = chr, .csa = csa, .bwt = bwt,
	 .occ = occ, .lut = lut, .dna = dna};

}


uN0_t
estimate_N0
(
 seed_t         L1,
 seed_t         L2,
 const index_t  idx,
 const double   mu
)
{

   int fwd = L1.end - L1.beg + 1;

   double a = 1. - pow(1-mu,fwd+1);
   double b = 1. - pow(1-mu,fwd);

   int N1 = log(log(a)/log(b)) / (log(b) - log(a));

   if (L1.beg == 0) {
      // No need to go reverse: the answer will be the same.
      // We switch gear and use Newton-Raphson iterations.
      long int m = L1.range.top - L1.range.bot + 1;
      double p = pow(1-mu, fwd);
      double N = N1 > m ? N1 : m;
      for (int j = 0 ; j < 8 ; j++) {
         double fN = digamma(N+1)-digamma(N-m+1) + log(1.-p);
         double dfN = trigamma(N+1)-trigamma(N-m+1);
         N = N - fN / dfN;
      }
      return (uN0_t) {mu, N, 1.};
   }

   int bwd = L2.end - L2.beg + 1;

   a = 1. - pow(1-mu,bwd+1);
   b = 1. - pow(1-mu,bwd);

   int N2 = log(log(a)/log(b)) / (log(b) - log(a));

   // Estimate N0.
   int N0 = (N1 + N2) / 2;

   // Compute p0 (prob of the min).
   int m = fwd < bwd ? fwd : bwd;
   size_t G = idx.chr->gsize * 2;
   // The term "(pow(1-pow(.25,m+1), G) - pow(1-pow(.25,m), G))" often
   // underflows. We divide every probability below by this term. This
   // simplifies "prob_2" and "prob_3", and in "prob_1", we can turn
   // "pow(1-pow(.25, m), G)" into "(exp(3G/4^(m+1))-1)^-1".

   // Probability of the event if both ends are duplicated.
   const double approx = exp(3*G * pow(.25, m+1)) - 1.;
   const double denom = approx < 1e-6 ? 1e-6 : approx;
   double prob_1 = mu * pow(1-mu, 2*m) / denom +
      (1-pow(1-mu,m)) * pow(1-mu,m) + mu * pow(1-mu, 2*m);
   // Probability of the event if one end is duplicated.
   double prob_2 = pow(1-mu,m);
   // Probability of the event if no end is duplicated.
   double prob_3 = (1.-pow(1-pow(.25,m+1), G));
   // Bayes formula with 1:9 prior for duplication. We assume
   // this because there is ~1:9 chance that the sequence
   // is duplicated with exactly 1 duplicate.
   double p0 = prob_1 / (prob_1 + prob_2 + 8*prob_3);

   // No nonsense. If just one error can explain the difference
   // between the estimates, assume that this error exists.
   if (L2.end + 2 >= L1.beg) p0 = 1.;

   return (uN0_t) {.06, N0, p0};


}


int
test_20mer_uniqueness
(
   const char    * seq,
   const index_t   idx
)
{

   size_t merid = 0;
   int mlen;

   // Look up the beginning (reverse) of the query in lookup table.
   for (mlen = 0 ; mlen < LUTK ; mlen++) {
      // Note: every "N" is considered a "A".
      uint8_t c = ENCODE[(uint8_t) seq[20-mlen-1]];
      merid = c + (merid << 2);
   }
   range_t range = idx.lut->kmer[merid];

   for ( ; mlen < 20 ; mlen++) {
      // When there are "N" in the reference, the estimation
      // must bail out because we cannot find the answer.
      if (NONALPHABET[(uint8_t) seq[20-mlen-1]])
         goto in_case_of_failure;
      int c = ENCODE[(uint8_t) seq[20-mlen-1]];
      range.bot = get_rank(idx.occ, c, range.bot - 1);
      range.top = get_rank(idx.occ, c, range.top) - 1;
      // If only 1 hit remains, we can bail out.
      if (range.bot >= range.top) return 1;
   }

   // Target is not unique.
   return 0;

in_case_of_failure:
   return -1;

}


int
cmpN0
(
   const void * a,
   const void * b
)
{
   uN0_t A = *(uN0_t *) a;
   uN0_t B = *(uN0_t *) b;

   return (A.N0 > B.N0) - (A.N0 < B.N0);

}

double
quality_low
(
   int      slen,
   seed_t * mem,
   uN0_t    uN0
)
{

  int len = mem->end - mem->beg + 1;
  double l = uN0.u*(1-PROB) + PROB*(1-uN0.u/3); // lambda
  double a = 0;
  int nends = (mem->beg == 0) + (mem->end == slen-1);
  switch (nends) {
    case 0:
      a = len * (1-PROB) * (1 - pow(1 - l*l*l/3*pow(1-l, len-1), uN0.N0));
      return a / (PROB + a);
    case 1:
      a = len * (1 - pow(1 - l*l/3*pow(1-l, len-1), uN0.N0));
      return a / (1 + a);
    case 2:
      a = len * PROB * (1 - pow(1 - l/3*pow(1-l, len-1), uN0.N0));
      return a / (1-PROB + a);
  }
  return 0./0.; // Oops...
}


//seedp_t
double
quality
(
         aln_t             aln,
   const char *            seq,  // Read.
         index_t           idx,
         uN0_t             uN0_read,
         pthread_mutex_t * mutex
)
{

   double slen = strlen(seq);
   // FIXME: assert is very weak (here only to declare 'uN0' on stack).
   assert(slen < 250);
   assert(slen >= GAMMA);

   int tot = (slen/10) - 2;
   int yes_max_evidence_N_is_0 = 1;
   double prob_p0 = .5;

   for (int s = 0 ; s <= slen-20 ; s += 10) {
      int unique = test_20mer_uniqueness(aln.refseq + s, idx);
      if (!unique) {
        yes_max_evidence_N_is_0 = 0;
        break;
      }
      prob_p0 *= .29; // = .94^20
   }

   if (yes_max_evidence_N_is_0 && slen >= 30) {
      // NB: for the super reads, we assume a frequency of 10%
      // in the genome. For Drosophila this is much more, for
      // human this is approximately half the value of Drosophila
      // and for pine this is 25x less. So this value is actually
      // genome-dependent. One way to get to it would just be to
      // count the proportion of reads that go to the super category
      // and plug the value in the formula. But we would need to
      // buffer the first ~10,000 reads to get to this estimate.
      // Those are the "super reads".
      const int mm = 1 + tot/2;
      if (aln.score == 0) {
         // Odd number of 20 nt block: we place a mm in the
         // event triplets. There are 11 positions on the first
         // and the last triplets, 10 positions on the internal
         // ones. The number of ways to choose the positions is
         // 11^2 * 10^{mm-2}. But we just say that every odd
         // block has 11 positions. Those are mm errors (occurrence
         // PROB), the other nucleotides are correct (occurrence
         // 1-PROB). The duplicate has compensating mutations
         // at those positions (occurrence u/3), the other
         // positions are not mutations (occurrence 1-u). We
         // divide by the probability that the read has score 0,
         // approximately equal to the probability that there is
         // no mutation. We also divide by the probability that
         // p0 is small, approximately 2/3 per segment.
         // Note that we divide 'tot' by 3 to approximate
         // the dependence between consecutive segments.
         // We also count a 1/10 probability that the read is
         // repeated with the specified level of u and with
         // only one extra copy.
         double u = mm / (double) slen; // Worst value of 'u'.
         return .1 * pow(10*PROB*u/3 / (1-PROB), mm) *
                     pow(1-u,slen-mm);
         // 0.1 * (11*pu/3)^mm * (1-u)^k-mm * (1-p)^k-mm / (1-p)^k
      }
      else if (aln.score == 1) {
         // Here we also have to consider the probability that there
         // would be a seed for the target. That makes it complex.
         // Case 1: the mismatch is in the second (or but-to-last)
         // segment, no further than GAMMA from the border. An
         // error can occur if the mismatch is an error (frequency p)
         // combined with an uncompensated mutation (2u/3). Other
         // hidden compensated errors are in the even segments.
         // Case 2: the mismatch is somewhere else. The most likely
         // scenario for an error is that the mismatch is a mutation
         // (frequency u) and that there are hidden and
         // compensated errors in even segments.
         // We also count a 1/10 probability that the read is
         // repeated with the specified level of u and with
         // only one extra copy.
         int errpos = aln.read_beg == 0 ? aln.read_end+1 : aln.read_beg-1;
         int case_1 = (errpos >= 10 && errpos < GAMMA) ||
            (errpos < slen-10 && errpos >= slen-GAMMA);
         if (case_1) {
            double u = mm / (double) slen; // Worst value of 'u'.
            return .1 * 2*u/3 * pow(10*PROB*u/3 / (1-PROB), mm-1) *
                  pow(1-u, slen-mm);
            // .1 * 2*pu/3 * (11*pu/3)^mm-1 * (1-u)^k-mm *
            // (1-p)^k-mm / p*(1-p)^k-1
         }
         else {
            double u = (mm+1) / (double) slen; // Worst value of 'u'.
            return .1 * 2*pow(10*PROB*u/3 / (1-PROB), mm) * u*(1-PROB) *
                  pow(1-u,slen-mm-1) / PROB;
            // .1 * 2*u*(1-p) * (11*pu/3)^mm * (1-u)^k-mm-1 *
            // (1-p)^k-mm-1 / p*(1-p)^k-1
         }
      }
      else if (aln.score == 2) {
         double u = mm / (double) slen; // Worst value of 'u'.
         // There are several ways this can be wrong, but we collapse
         // it to an "average" case with compensated errors in all
         // even segments except one. This last segment contains a
         // mutation and an uncompensated error.
         if (slen < 50);
         else if (slen == 50) {
            // Special case where there are 231 possibilities in
            // the order mutation, error, compensated error, plus
            // 84 possibilities in the order error, mutation,
            // compensated error (two times by symmetry).
            return 0.2 * 2*315*u/3*u*(1-u/3)*pow(1-u,47) /
                  slen / (slen-1);
         }
         else {
            return mm * 6*pow(11*PROB*u/3 / (1-PROB), mm) *
                  pow((1-PROB)/PROB,2) * pow(1-u,slen-mm) /
                  slen / (slen-1);
            // 0.1 *mm*11*10*(1-p)*u * p*(1-u/3) * (11*pu/3)^{mm-1} *
            // (1-u)^k-mm-1 * (1-p)^k-mm-1 / k*(k-1)/2*p^2*(1-p)^k-2
         }
      }
   }

   double u = uN0_read.u;
   int N0 = uN0_read.N0;
   double p0 = uN0_read.p0;

   // Estimate N0 on the hit.
   seed_t L1, L2;
   extend_L1L2(aln.refseq, slen, idx, &L1, &L2);
   uN0_t uN0_hit = estimate_N0(L1, L2, idx, u);
   if (uN0_hit.N0 * uN0_hit.p0 > N0 * p0) {
      N0 = uN0_hit.N0;
      p0 = uN0_hit.p0;
   }

   if (yes_max_evidence_N_is_0) {
     N0 = 1;
     prob_p0 = (prob_p0 + p0) / 2;
   }
   else {
     prob_p0 = p0;
   }

   if (aln.score == 0) {
     // Special case for perfect alignment score.
     const double cond = pow(1-PROB, slen);
     double p1 = slen * PROB * pow(1-PROB, slen-1);
     double p2 = u/3. * pow(1-u, slen-1);
     return p0 * p1 * (1. - pow(1-p2, N0)) / cond;
   }

   pthread_mutex_lock(mutex);
   double poff = auto_mem_seed_offp(slen, u, N0);
   pthread_mutex_unlock(mutex);

   // Weight of the evidence if mapping is correct...
   double A = aln.score * log(PROB) + (slen-aln.score) * log(1-PROB);
   // ... if mapping is on a duplicate...
   double B = aln.score * log(u) + (slen-aln.score) * log(1-u);
   // ... and if mapping is random.
   int naln = (aln.read_beg == 0 || aln.read_end == slen-1) ?
      slen - aln.read_end + aln.read_beg - 2 :
      slen - aln.read_end + aln.read_beg - 3;
   double C = aln.score * log(.75)  + (naln-aln.score) * log(.25);

   if (A - B > 1.) A = B + 1.;

   // Here 'term2' uses non-informative prios instead of Sesame
   // priors for the sake of simplicity. In practice, 'term2' is
   // either very close to 0 or very close to 1 and the priors
   // do not matter.
   double term1 = prob_p0 * poff / ( poff + exp(A-B)*(1-poff) );
   double term2 = aln.score < slen / 5 ? 0 : 1. / (1. + exp(A-C));

   return term1 + term2 > 1. ? 1. : term1 + term2;

}

void
print_sam_header
(
   const index_t idx
)
{
   for (int i = 0 ; i < idx.chr->nchr; i++) {
     size_t sz = i >= idx.chr->nchr-1 ?
       idx.chr->gsize - idx.chr->start[i] :
       idx.chr->start[i+1] - idx.chr->start[i];
     fprintf(stdout, "@SQ\tSN:%s\tLN:%ld\n", idx.chr->name[i], sz);
   }
}


void
remove_N
(
   char * seq
)
{
   const char DNA[4] = "gatc";
   static int c = 0;
   for (char * p = strchr(seq, 'N') ; p != NULL ; p = strchr(p, 'N')) {
      *p = DNA[c++ % 4];
   }
}


int
parse_read
(
        FILE      * inputf,
   enum fmt_t       format,
	     wstack_t ** stack
)
{
  size_t n = 0;
  ssize_t len = 0;
  read_t * read = calloc(1, sizeof(read_t));
  exit_error(read == NULL);

  if (format == fasta) {
    // Read 2 lines.
    //name
    len = getline(&(read->name), &n, inputf);
    if (len == -1) {
       goto end_of_file_return;
    }
    // Replace empty space by '\0'.
    strtok(read->name, " \t\n");

    //seq
    n = 0;
    len = getline(&(read->seq), &n, inputf);
    if (len == -1) {
       fprintf(stderr, "format error in fasta file\n");
       goto file_error_return;
    }
    if (read->seq[len-1] == '\n') read->seq[len-1] = 0;
    // Replace new line by '\0'.
    strtok(read->seq, "\n");
    remove_N(read->seq);

    //no quality
    read->phred = strdup("*");

    push(read, stack);
    return 1;
  }

  if (format == fastq) {
    // name
    len = getline(&(read->name), &n, inputf);
    if (len == -1) {
       goto end_of_file_return;
    }
    // Replace empty space by '\0'.
    strtok(read->name, " \t\n");

    // seq
    n = 0;
    len = getline(&(read->seq), &n, inputf);
    if (len == -1) {
       fprintf(stderr, "format error in fastq file\n");
       goto file_error_return;
    }
    // Replace new line by '\0'.
    strtok(read->seq, "\n");
    remove_N(read->seq);

    // + (skip line)
    if (fscanf(inputf, "%*[^\n]\n") < 0) {
       fprintf(stderr, "format error in fastq file\n");
       goto file_error_return;
    }

    // phred
    n = 0;
    len = getline(&(read->phred), &n, inputf);
    if (len == -1) {
       fprintf(stderr, "format error in fastq file\n");
       goto file_error_return;
    }
    // Replace new line by '\0'.
    strtok(read->phred, "\n");

    push(read, stack);
    return 1;

  }

  // Wrong format.
  free(read);
  return -1;

end_of_file_return:
  free(read->name);
  free(read->seq);
  free(read->phred);
  free(read);
  return 0;

file_error_return:
  free(read->name);
  free(read->seq);
  free(read->phred);
  free(read);
  return -1;

}

void *
batch_map
(
 void * arg
)
{
   batch_t * batch  = (batch_t *)arg;
   index_t idx = batch->idx;
   int * act_threads = batch->act_threads;
   char * outstr;

   while (1) {
      // Two blocks:
      // 1. No reads available yet
      // 2. No threads available

      // 1. Reads available
      pthread_mutex_lock(batch->mutex);
      while (batch->status != map) {
         if (batch->status != die) {
            pthread_cond_wait(batch->cond_reads, batch->mutex);
         } else {
            pthread_mutex_unlock(batch->mutex);
            return NULL;
         }
      }

      // 2. Wait for available thread slot
      while (*act_threads >= MAXTHREADS)
         pthread_cond_wait(batch->cond_threads, batch->mutex);

      // Take thread slot
      *act_threads = *act_threads + 1;
      pthread_mutex_unlock(batch->mutex);

      // Reset output
      batch->output->pos = 0;
      for (size_t i = 0; i < batch->reads->pos; i++) {

         read_t * read = (read_t *)batch->reads->ptr[i];
         size_t rlen = strlen(read->seq);

         // DEBUG //
         fprintf(stderr, "%s\n", read->seq);

         // Compute L1, L2 and MEMs.
         seed_t L1, L2;
         extend_L1L2(read->seq, rlen, idx, &L1, &L2);

         // Compute seeds.
         wstack_t * seeds = mem_seeds(read->seq, idx, GAMMA);

         // Return if no seeds were found
         if (seeds->pos == 0) {
            // Did not find anything.
            free(seeds);
            // Output in sam format.
            int olen = snprintf(NULL, 0, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",
                  read->name+1, read->seq, read->phred);
            outstr = malloc(olen+1);
            exit_error(outstr == NULL);
            sprintf(outstr, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",
                  read->name+1, read->seq, read->phred);
            batch->output->ptr[batch->output->pos++] = outstr;
            continue;
         }

         // Compute N(L1,L2)
         const double lambda = (1-PROB)*.06 + PROB*(1-.06/3);
         uN0_t uN0 = estimate_N0(L1, L2, idx, lambda);

         // Quick mode: only align longest MEMs
         seed_t *longest_mem = NULL;
         if (uN0.N0 > QUICK_DUPLICATES) {
            longest_mem = filter_longest_mem(seeds);
         }

         alnstack_t * alst = mapread(
            seeds,
            read->seq,
            idx,
            rlen,
            batch->lineid + i
         );

         // Did not find anything.
         // TODO: use skip seeds in this case.
         if (alst->pos == 0) {
            int olen = snprintf(NULL, 0, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",
                  read->name+1, read->seq, read->phred);
            outstr = malloc(olen+1);
            exit_error(outstr == NULL);
            sprintf(outstr, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",
                  read->name+1, read->seq, read->phred);
            batch->output->ptr[batch->output->pos++] = outstr;
            free(alst);
            // Output in sam format.
            // Free seeds.
            for (size_t i = 0; i < seeds->pos; i++) {
               seed_t * s = (seed_t *) seeds->ptr[i];
               free(s->sa);
               free(s);
            }
            free(seeds);
            // Free read.
            free(read->name);
            free(read->seq);
            free(read->phred);
            free(read);
            continue;
         }

         // Pick a top alignment at "random".
         aln_t a = alst->aln[(batch->lineid + i) % alst->pos];

         if (alst->pos == 1) {
            a.qual = uN0.N0 > QUICK_DUPLICATES ?
               quality_low(rlen, longest_mem, uN0) :
               quality(a, read->seq, idx, uN0, batch->mutex_sesame);
         }
         else {
            a.qual = 1-1./alst->pos;
         }

        // If the results are so-so, use skip seeds to see if we
        // can find something better.
        int lo_mapq = a.qual > .001;
        int hi_err = a.score > 1;
        int lo_N = longest_mem == NULL;
        int no_tie = alst->pos == 1;
        if (lo_mapq && hi_err && lo_N && no_tie) {
           // Run skip-9 seeding of minimum size 19. This scheme has a
           // ~1% chance of not finding the target on reads of 50 nt.
           // In addition, it allows us to ascertain that 50 nt reads
           // have at least two errors if they have only one seed.
           wstack_t * skipseeds = skip_seeds(read->seq, idx, 19, 9);
           if (skipseeds->pos > 0) {
              alnstack_t * salst = remap_with_skip_seeds(
                  skipseeds,
                  alst,  // Contains loci aligned with MEMS.
                  read->seq,
                  idx,
                  a.score,
                  batch->lineid + i
              );
              if (salst->pos > 0) {
                 // Found a better hit.
                 for (size_t i = 0 ; i < salst->pos ; i++) {
                    aln_t skipa = salst->aln[i];
                    fprintf(stdout, ">>Better score: %d\n", skipa.score);
                    free(salst->aln[i].refseq);
                 }
              }
              free(salst);
           }
           // Free seeds.
           for (size_t i = 0 ; i < skipseeds->pos ; i++) {
              seed_t * s = (seed_t *) skipseeds->ptr[i];
              free(s->sa);
              free(s);
           }
           free(skipseeds);
        }

         // Report mapping results
         pos_t pos = get_pos(a.refpos, idx.chr);
         // Output in sam format.
         int bits = pos.strand ? 0 : 16;
         size_t leftpos = pos.strand ? pos.pos : pos.pos - rlen+1;
         int olen = snprintf(NULL, 0, "%s\t%d\t%s\t%ld\t%d\t%ldM\t*\t0\t0\t%s\t%s\tXS:i:%d\n",
               read->name+1, bits, pos.rname, leftpos, (int) (-10*log10(a.qual)),
               rlen, read->seq, read->phred, a.score);

         outstr = malloc(olen+1);
         exit_error(outstr == NULL);
         sprintf(outstr, "%s\t%d\t%s\t%ld\t%d\t%ldM\t*\t0\t0\t%s\t%s\tXS:i:%d\n",
               read->name+1, bits, pos.rname, leftpos, (int) (-10*log10(a.qual)),
               rlen, read->seq, read->phred, a.score);
         batch->output->ptr[batch->output->pos++] = outstr;

         // Free seeds
         for (size_t i = 0; i < seeds->pos; i++) {
            seed_t * s = (seed_t *) seeds->ptr[i];
            free(s->sa);
            free(s);
         }
         free(seeds);

         // Free alignments
         for(size_t i = 0; i < alst->pos; i++) free(alst->aln[i].refseq);
         free(alst);

         // Free read
         free(read->name);
         free(read->seq);
         free(read->phred);
         free(read);
      }

      // Decrease active threads
      pthread_mutex_lock(batch->mutex);
      batch->status = output;
      *act_threads = *act_threads - 1;
      pthread_cond_signal(batch->cond_threads);
      pthread_cond_signal(batch->cond_writer);
      pthread_mutex_unlock(batch->mutex);
   }
   return NULL;
}


void *
mt_write
(
void * arg
)
{
   writerarg_t * warg = (writerarg_t *)arg;

   batch_t * batch = warg->first_batch;

   while(1) {
      while (batch && batch->status == output) {
         // Write all output strings
         for (size_t i = 0; i < batch->output->pos; i++) {
            fprintf(stdout, "%s", (char *)batch->output->ptr[i]);
            free(batch->output->ptr[i]);
         }

         // Wait until next batch is assigned
         pthread_mutex_lock(warg->mutex);
         while (batch->next_batch == (void *)-1)
            pthread_cond_wait(warg->cond_writer, warg->mutex);

         // Got new batch, release old one
         batch->status = idle;
         batch = batch->next_batch;

         // Signal reader, batch is idle
         pthread_cond_signal(warg->cond_reader);
         pthread_mutex_unlock(warg->mutex);
      }

      if (!batch) break;

      pthread_mutex_lock(warg->mutex);
      pthread_cond_wait(warg->cond_writer, warg->mutex);
      pthread_mutex_unlock(warg->mutex);
   }

   return NULL;
}

void
mt_read
(
   const char * indexfname,
   const char * readsfname
)
{
   size_t maxlen = 0; // Max 'k' value for seeding probabilities.

   fprintf(stderr, "loading index... ");
   // Load index files.
   index_t idx = load_index(indexfname);

   fprintf(stderr, "done.\n");

   FILE * inputf = fzopen(readsfname, "r");
   if (inputf == NULL) exit_cannot_open(readsfname);

   // Set the input type.
   enum fmt_t format = unset;

   // Look at the first character.
   char first = getc(inputf);
   if (first == '>') format = fasta;
   if (first == '@') format = fastq;
   ungetc(first, inputf);
   print_sam_header(idx);

   int success = 0;

   size_t sz = 256, line = 0;
   char * buff = malloc(sz);
   exit_error(buff == NULL);

   // Multithreading variables
   int BATCHSIZE = 10000;
   int act_threads = 0;

   // Create mutex and monitors
   pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
   pthread_mutex_t mutex_sesame = PTHREAD_MUTEX_INITIALIZER;
   pthread_cond_t  cond_reader  = PTHREAD_COND_INITIALIZER;
   pthread_cond_t  cond_writer  = PTHREAD_COND_INITIALIZER;
   pthread_cond_t  cond_threads = PTHREAD_COND_INITIALIZER;
   pthread_cond_t * cond_reads  = malloc(2*MAXTHREADS*sizeof(pthread_cond_t));
   exit_error(!cond_reads);

   // Allocate batches
   pthread_t * worker = malloc(2*MAXTHREADS*sizeof(pthread_t));
   batch_t ** batch = malloc(2*MAXTHREADS*sizeof(batch_t));
   exit_error(batch == NULL || worker == NULL);

   // Allocate and fill static batch info
   for (int i = 0; i < 2*MAXTHREADS; i++) {
      // Alloc batch
      batch[i] = malloc(sizeof(batch_t));
      exit_error(batch[i] == NULL);
      // Initialize static info
      batch[i]->idx          = idx;
      batch[i]->reads        = stack_new(BATCHSIZE);
      batch[i]->output       = stack_new(BATCHSIZE);
      batch[i]->status       = idle;
      batch[i]->act_threads  = &act_threads;
      batch[i]->mutex        = &mutex;
      batch[i]->mutex_sesame = &mutex_sesame;
      batch[i]->cond_reader  = &cond_reader;
      batch[i]->cond_writer  = &cond_writer;
      batch[i]->cond_reads   = cond_reads+i;
      batch[i]->cond_threads = &cond_threads;

      exit_error(!batch[i]->reads || !batch[i]->output);

      // Spawn thread
      pthread_cond_init(cond_reads+i, NULL);
      pthread_create(worker+i, NULL, batch_map, batch[i]);
   }

   // Spawn writer thread
   writerarg_t warg = (writerarg_t){batch[0], &mutex, &cond_writer, &cond_reader};
   pthread_t writer_thread;
   pthread_create(&writer_thread, NULL, mt_write, &warg);

   int last_worker = -1;
   while (1) {
      // Find idle worker
      int w = -1;
      pthread_mutex_lock(&mutex);
      for (int i = (last_worker+1)%(2*MAXTHREADS), j = 0; j < 2*MAXTHREADS; i = (i+1)%(2*MAXTHREADS), j++) {
         if (batch[i]->status == idle && i != last_worker) {
            w = i;
            break;
         }
      }
      if (w == -1) {
         pthread_cond_wait(&cond_reader, &mutex);
         pthread_mutex_unlock(&mutex);
         continue;
      }
      pthread_mutex_unlock(&mutex);

      // Thread is idle, can modify without race conditions
      batch[w]->reads->pos = 0;
      while ((success = parse_read(inputf, format, &(batch[w]->reads)))) {
         // Set sesame static params
         size_t rlen = strlen(((read_t *)batch[w]->reads->ptr[batch[w]->reads->pos-1])->seq);
         if (rlen > maxlen) {
            maxlen = rlen;
            pthread_mutex_lock(&mutex_sesame);
            sesame_set_static_params(GAMMA, maxlen, PROB);
            pthread_mutex_unlock(&mutex_sesame);
         }

         line++;
         if (batch[w]->reads->pos == batch[w]->reads->max)
            break;
      }

      if (success < 0) {
         fprintf(stderr, "error while reading input file (line %ld)\n", line+1);
         exit(EXIT_FAILURE);
      }

      // Batch linked list
      if (last_worker >= 0) {
         pthread_mutex_lock(&mutex);
         batch[last_worker]->next_batch = batch[w];
         pthread_cond_signal(&cond_writer);
         pthread_mutex_unlock(&mutex);
      }

      // Flag next batch
      if (success) {
         batch[w]->next_batch = (void *) -1;
         last_worker = w;
      } else {
         // Last batch
         batch[w]->next_batch = NULL;
         last_worker = -1;
      }

      // Assign line id for output consistency
      batch[w]->lineid = line - batch[w]->reads->pos;

      // Signal worker
      pthread_mutex_lock(&mutex);
      batch[w]->status = map;
      pthread_cond_signal(batch[w]->cond_reads);
      pthread_mutex_unlock(&mutex);

      // Break on last batch
      if (last_worker == -1) break;
   }

   // Wait for writer.
   pthread_join(writer_thread, NULL);

   // Signal worker threads
   pthread_mutex_lock(&mutex);
   for (int i = 0; i < 2*MAXTHREADS; i++) {
      batch[i]->status = die;
      pthread_cond_signal(batch[i]->cond_reads);
   }
   pthread_mutex_unlock(&mutex);

   // Free memory
   for (int i = 0; i < 2*MAXTHREADS; i++) {
      pthread_join(worker[i], NULL);
      free(batch[i]->reads);
      free(batch[i]->output);
      free(batch[i]);
   }
   free(buff);
   free(worker);
   free(batch);
   free(cond_reads);
   fclose(inputf);
   free_index_chr(idx.chr);
   sesame_clean();
}


int
main
(
int argc,
char ** argv
)
{

   // Sanity checks.
   if (argc < 2) {
      fprintf(stderr, "First argument must be \"--index\" or index file.\n");
      exit(EXIT_FAILURE);
   }

   if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) {
      fprintf(stderr, "MEM Mapper Prototype version %s\n", VERSION);
      fprintf(stderr, "%s", HELP_MSG);
   }
   else if (strcmp(argv[1], "--index") == 0) {
      if (argc < 3) {
         fprintf(stderr, "Specify file to index.\n");
         exit(EXIT_FAILURE);
      }
      build_index(argv[2]);
   }
   else {
      if (argc < 3) {
         fprintf(stderr, "error: not enough arguments.\n");
         fprintf(stderr, "%s", HELP_MSG);
         exit(EXIT_FAILURE);
      }

      char * index_path = NULL;
      char * reads_path = NULL;

      for (int i = 1; i < argc; i++) {
         if (argv[i][0] == '-') {
            if (argv[i][1] == 'q') {
               SKIPQUAL = strtod(argv[++i], NULL);
            } else if (argv[i][1] == 'e') {
               PROB = strtod(argv[++i], NULL);
               if (PROB <= 0 || PROB >= 1) {
                  fprintf(stderr, "Sequencing error must be in (0,1).\n");
                  exit(EXIT_FAILURE);
               }
            } else if (argv[i][1] == 't') {
               MAXTHREADS = atoi(argv[++i]);
               if (MAXTHREADS <= 0) {
                  fprintf(stderr, "Max threads must be greater than 0.\n");
                  exit(EXIT_FAILURE);
               }
            }
         } else {
            if (!index_path) index_path = argv[i];
            else if (!reads_path) reads_path = argv[i];
            else {
               fprintf(stderr, "error: too many arguments.\n");
               fprintf(stderr, "%s", HELP_MSG);
               exit(EXIT_FAILURE);
            }
         }
      }

      if (!index_path || !reads_path) {
         fprintf(stderr, "error: not enough arguments.\n");
         fprintf(stderr, "%s", HELP_MSG);
         exit(EXIT_FAILURE);
      }

      mt_read(index_path, reads_path);

   }
}


double
digamma
(
   double x
)
{
  double neginf = -INFINITY;
  static const double c = 12,
    digamma1 = -0.57721566490153286,
    trigamma1 = 1.6449340668482264365, /* pi^2/6 */
    s = 1e-6,
    s3 = 1./12,
    s4 = 1./120,
    s5 = 1./252,
    s6 = 1./240,
    s7 = 1./132;
    //s8 = 691./32760,
    //s9 = 1./12,
    //s10 = 3617./8160;
  double result;
  /* Illegal arguments */
  if((x == neginf) || isnan(x)) {
    return NAN;
  }
  /* Singularities */
  if((x <= 0) && (floor(x) == x)) {
    return neginf;
  }
  /* Negative values */
  /* Use the reflection formula (Jeffrey 11.1.6):
   * digamma(-x) = digamma(x+1) + pi*cot(pi*x)
   *
   * This is related to the identity
   * digamma(-x) = digamma(x+1) - digamma(z) + digamma(1-z)
   * where z is the fractional part of x
   * For example:
   * digamma(-3.1) = 1/3.1 + 1/2.1 + 1/1.1 + 1/0.1 + digamma(1-0.1)
   *               = digamma(4.1) - digamma(0.1) + digamma(1-0.1)
   * Then we use
   * digamma(1-z) - digamma(z) = pi*cot(pi*z)
   */
  if(x < 0) {
    return digamma(1-x) + M_PI / tan(-M_PI*x);
  }
  /* Use Taylor series if argument <= S */
  if(x <= s) return digamma1 - 1/x + trigamma1*x;
  /* Reduce to digamma(X + N) where (X + N) >= C */
  result = 0;
  while(x < c) {
    result -= 1/x;
    x++;
  }
  /* Use de Moivre's expansion if argument >= C */
  /* This expansion can be computed in Maple via asympt(Psi(x),x) */
  if(x >= c) {
    double r = 1/x;
    result += log(x) - 0.5*r;
    r *= r;
    result -= r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
  }
  return result;
}

/* The trigamma function is the derivative of the digamma function.

   Reference:

    B Schneider,
    Trigamma Function,
    Algorithm AS 121,
    Applied Statistics,
    Volume 27, Number 1, page 97-99, 1978.

    From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
    (with modification for negative arguments and extra precision)
*/


double trigamma(double x)
{
  double neginf = -INFINITY,
    small = 1e-4,
    large = 8,
    trigamma1 = 1.6449340668482264365, /* pi^2/6 = Zeta(2) */
    tetragamma1 = -2.404113806319188570799476,  /* -2 Zeta(3) */
    b2 =  1./6,  /* B_2 */
    b4 = -1./30, /* B_4 */
    b6 =  1./42, /* B_6 */
    b8 = -1./30, /* B_8 */
    b10 = 5./66; /* B_10 */
  double result;
  /* Illegal arguments */
  if((x == neginf) || isnan(x)) {
    return NAN;
  }
  /* Singularities */
  if((x <= 0) && (floor(x) == x)) {
    return neginf;
  }
  /* Negative values */
  /* Use the derivative of the digamma reflection formula:
   * -trigamma(-x) = trigamma(x+1) - (pi*csc(pi*x))^2
   */
  if(x < 0) {
    result = M_PI/sin(-M_PI*x);
    return -trigamma(1-x) + result*result;
  }
  /* Use Taylor series if argument <= small */
  if(x <= small) {
    return 1/(x*x) + trigamma1 + tetragamma1*x;
  }
  result = 0;
  /* Reduce to trigamma(x+n) where ( X + N ) >= B */
  while(x < large) {
    result += 1/(x*x);
    x++;
  }
  /* Apply asymptotic formula when X >= B */
  /* This expansion can be computed in Maple via asympt(Psi(1,x),x) */
  if(x >= large) {
    double r = 1/(x*x);
    result += 0.5*r + (1 + r*(b2 + r*(b4 + r*(b6 + r*(b8 + r*b10)))))/x;
  }
  return result;
}
