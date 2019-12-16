#define _GNU_SOURCE

#include <math.h>

#include "bwt.h"
#include "map.h"
#include "divsufsort.h"
#include "sesame.h"

#define GAMMA 17
#define PROBDEFAULT 0.01
#define SKIPQUALDEFAULT 10
#define QUICK_DUPLICATES 20

// Index parameters
#define OCC_INTERVAL_SIZE 8

// ------- Machine-specific code ------- //
// The 'mmap()' option 'MAP_POPULATE' is available
// only on Linux and from kern version 2.6.23.
#if __linux__
  #include <linux/version.h>
  #if LINUX_VERSION_CODE > KERNEL_VERSION(2,6,22)
    #define _MAP_POPULATE_AVAILABLE
  #endif
#endif

#ifdef _MAP_POPULATE_AVAILABLE
  #define MMAP_FLAGS (MAP_PRIVATE | MAP_POPULATE)
#else
  #define MMAP_FLAGS MAP_PRIVATE
#endif

// ------- Machine-specific code ------- //
// The 'mmap()' option 'MAP_POPULATE' is available
// only on Linux and from kern version 2.6.23.
#if __linux__
  #include <linux/version.h>
  #if LINUX_VERSION_CODE > KERNEL_VERSION(2,6,22)
    #define _MAP_POPULATE_AVAILABLE
  #endif
#endif

#ifdef _MAP_POPULATE_AVAILABLE
  #define MMAP_FLAGS (MAP_PRIVATE | MAP_POPULATE)
#else
  #define MMAP_FLAGS MAP_PRIVATE
#endif


// Error-handling macro.
#define exit_cannot_open(x) \
   do { fprintf(stderr, "cannot open file '%s' %s:%d:%s()\n", (x), \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); } while(0)

#define exit_error(x) \
   do { if (x) { fprintf(stderr, "error: %s %s:%d:%s()\n", #x, \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); }} while(0)

typedef struct uN0_t uN0_t;
typedef struct seedp_t seedp_t;

struct uN0_t {
   double u;
   size_t N0;
   double lev;
};

struct seedp_t {
   double off;
   double nul;
};

static double PROB = PROBDEFAULT;
static double SKIPQUAL = SKIPQUALDEFAULT;

char* HELP_MSG =
   "usage: smmfdp ([-e 0.01] [-q 10] index-file file.fasta | --index index-file)\n"
   "\n"
   "mapping options:\n"
   "  -e  sets the expected error rate (default: 0.01)\n"
   "  -q  sets the threshold qual1ty to use skip seeds (default 10)\n"
   "\n";
   

void
build_index
(
   const char * fname
)
{

   // Open fasta file.
   FILE * fasta = fopen(fname, "r");
   if (fasta == NULL) exit_cannot_open(fname);

   // Aux variables for file writing.
   char * data;
   ssize_t ws;
   size_t sz;

   char buff[256];
   size_t gsize;

   // Read and normalize genome
   sprintf(buff, "%s.chr", fname);
   fprintf(stderr, "reading genome... ");
   char * genome = normalize_genome(fasta, buff, &gsize);
   fprintf(stderr, "done.\n");

   fclose(fasta);

   fprintf(stderr, "compressing nucleotides... ");
   char * dna = compress_genome(genome, gsize);

   // Write the compressed genome
   sprintf(buff, "%s.dna", fname);
   int fdna = creat(buff, 0644);
   if (fdna < 0) exit_cannot_open(buff);
   
   ws = 0;
   sz = gsize/4+1;
   data = (char *) dna;
   while (ws < sz) ws += write(fdna, data + ws, sz - ws);
   close(fdna);

   free(dna);
   fprintf(stderr, "done.\n");

   fprintf(stderr, "creating suffix array... ");
   int64_t * sa = compute_sa(genome);
   fprintf(stderr, "done.\n");

   fprintf(stderr, "creating bwt... ");
   bwt_t * bwt = create_bwt(genome, sa);
   fprintf(stderr, "done.\n");

   free(genome);

   fprintf(stderr, "creating Occ table... ");
   occ_t * occ = create_occ(bwt, OCC_INTERVAL_SIZE);
   fprintf(stderr, "done.\n");

   fprintf(stderr, "filling lookup table... ");
   lut_t * lut = malloc(sizeof(lut_t));
   fill_lut(lut, occ, (range_t) {.bot=1, .top=gsize}, 0, 0);
   fprintf(stderr, "done.\n");

   fprintf(stderr, "compressing suffix array... ");
   csa_t * csa = compress_sa(sa);
   fprintf(stderr, "done.\n");

   free(sa);

   // Write the compressed suffix array file.
   sprintf(buff, "%s.sa", fname);
   int fsar = creat(buff, 0644);
   if (fsar < 0) exit_cannot_open(buff);

   ws = 0;
   sz = sizeof(csa_t) + csa->nint64 * sizeof(int64_t);
   data = (char *) csa;
   while (ws < sz) ws += write(fsar, data + ws, sz - ws);
   close(fsar);

   free(csa);

   // Write the Burrows-Wheeler transform.
   sprintf(buff, "%s.bwt", fname);
   int fbwt = creat(buff, 0644);
   if (fbwt < 0) exit_cannot_open(buff);

   ws = 0;
   sz = sizeof(bwt_t) + bwt->nslots * sizeof(uint8_t);
   data = (char *) bwt;
   while (ws < sz) ws += write(fbwt, data + ws, sz - ws);
   close(fbwt);

   free(bwt);

   // Write the Occ table.
   sprintf(buff, "%s.occ", fname);
   int focc = creat(buff, 0644);
   if (focc < 0) exit_cannot_open(buff);

   ws = 0;
   sz = sizeof(occ_t) + occ->occ_size*SIGMA*sizeof(uint64_t);
   data = (char *) occ;
   while (ws < sz) ws += write(focc, data + ws, sz - ws);
   close(focc);

   free(occ);

   // Write the lookup table
   sprintf(buff, "%s.lut", fname);
   int flut = creat(buff, 0644);
   if (flut < 0) exit_cannot_open(buff);

   ws = 0;
   sz = sizeof(lut_t);
   data = (char *) lut;
   while (ws < sz) ws += write(flut, data + ws, sz - ws);
   close(flut);

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
estimate_uN0
(
   const char    * seq,
   const index_t   idx
)
{

   const size_t n = 21;
   const double MU[3] = { .06, .04, .02 };

   size_t L1, L2, R1, R2;
   L1 = L2 = R1 = R2 = 0;
   int mlen;

   size_t merid = 0;

   // Look up the beginning (reverse) of the query in lookup table.
   for (mlen = 0 ; mlen < LUTK ; mlen++) {
      // Note: every "N" is considered a "A".
      uint8_t c = ENCODE[(uint8_t) seq[30-mlen-1]];
      merid = c + (merid << 2);
   }
   range_t range = idx.lut->kmer[merid];

   for ( ; mlen < 30 ; mlen++) {
      // When there are "N" in the reference, the estimation
      // must bail out because we cannot find the answer.
      if (NONALPHABET[(uint8_t) seq[30-mlen-1]])
         goto in_case_of_failure;
      int c = ENCODE[(uint8_t) seq[30-mlen-1]];
      range.bot = get_rank(idx.occ, c, range.bot - 1);
      range.top = get_rank(idx.occ, c, range.top) - 1;
      // If only 1 hit remains, we can bail out.
      if (range.bot == range.top) break;
      if (mlen == 20) {
         L1 = range.top - range.bot;
         L2 = n * L1;
      }
      if (mlen > 20) {
         L2 += range.top - range.bot;
      }
   }

   L1 -= (range.top - range.bot);

   merid = 0;

   // Look up the beginning (forward) of the query in lookup table.
   for (mlen = 0 ; mlen < LUTK ; mlen++) {
      // Note: every "N" is considered a "T".
      uint8_t c = REVCMP[(uint8_t) seq[mlen]];
      merid = c + (merid << 2);
   }
   range = idx.lut->kmer[merid];

   for ( ; mlen < 30 ; mlen++) {
      // When there are "N" in the reference, the estimation
      // must bail out because we cannot find the answer.
      if (NONALPHABET[(uint8_t) seq[mlen]])
         goto in_case_of_failure;
      int c = REVCMP[(uint8_t) seq[mlen]];
      range.bot = get_rank(idx.occ, c, range.bot - 1);
      range.top = get_rank(idx.occ, c, range.top) - 1;
      // If only 1 hit remains, we can bail out.
      if (range.bot == range.top) break;
      if (mlen == 20) {
         R1 = range.top - range.bot;
         R2 = n * R1;
      }
      if (mlen > 20) {
         R2 += range.top - range.bot;
      }
   }

   R1 -= (range.top - range.bot);

   double loglik = -INFINITY;
   size_t best_N0 = 0;
   double best_mu = 0.0;

   double C[3] = {
       (1 - pow(1-.06,n)) / (n * pow(1-.06,n-1)),
       (1 - pow(1-.04,n)) / (n * pow(1-.04,n-1)),
       (1 - pow(1-.02,n)) / (n * pow(1-.02,n-1)),
   };

   for (int iter = 0 ; iter < 3 ; iter++) {

      double mu = MU[iter];

      double L3 = L2 / (1-mu) - L1 / mu;
      double R3 = R2 / (1-mu) - R1 / mu;

      // Compute the number of duplicates.
      double N0 = (L1+R1 + C[iter]*(L3+R3)) / 2.0;
      if (N0 < 1.0) N0 = 1.0;

      // TODO: precompute some of this?
      double tmp = 2 * lgamma(N0+1) + (L1+R1) * log(mu) +
         (L2+R2) * log(1-mu) + (2*N0-L1-R1) * log(1-pow(1-mu,n)) -
         (lgamma(N0-L1+1) + lgamma(N0-R1+1));

      if (tmp < loglik) {
         break;
      }

      loglik = tmp;
      best_N0 = round(N0);
      best_mu = mu;

   }

   size_t G = idx.chr->gsize;
   double H0l = 2 * lgamma(G+1) + (L1+R1) * log(.75) +
      (L2+R2) * log(.25) + (2*G-L1-R1) * log(1-pow(.25,n)) -
      (lgamma(G-L1+1) + lgamma(G-R1+1));

   return (uN0_t) { best_mu, best_N0, H0l-loglik };

in_case_of_failure:
   // Return an impossible value.
   return (uN0_t) { 0. / 0., -1,  0. / 0.};

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


//seedp_t
double
quality
(
         aln_t   aln,
   const char *  seq,
         index_t idx
)
{

   double slen = strlen(seq);
   // FIXME: assert is very weak (here only
   // to declare 'uN0' on stack).
   assert(slen < 250);
   assert(slen >= GAMMA);

   uN0_t uN0[50] = {{0}};

   // Assume the worst (many similar duplicates).
   const double worst_case_mu = 0.02;
   const size_t worst_case_N0 = 2147483645;

   int tot = 0;
   double lev = 0;
   int yes_max_evidence_N_is_0 = 1;

   // Estimate N and mu.
   for (int s = 0 ; s <= slen-30 ; s += 10, tot++) {
      /*
      if (aln.read_end < s+10 || aln.read_beg > s+19) {
         // Estimate only at the site of the seed.
         // TODO: check whether this is really the best option.
         uN0[tot] = (uN0_t) { worst_case_mu, worst_case_N0, 0 };
         continue;
      }
      */
      uN0[tot] = estimate_uN0(aln.refseq + s, idx);
      // Failsafe when estimation fails (e.g., because of "N").
      if (uN0[tot].u != uN0[tot].u) {
         uN0[tot] = (uN0_t) { worst_case_mu, worst_case_N0, 0 };
      }
      // The magic value '0.6367' is the maximum possible value
      // of the evidence that N0 = 0 for each segment.
      if (uN0[tot].lev < 0.6367) yes_max_evidence_N_is_0 = 0;
      lev += uN0[tot].lev; // Collect total evidence.
   }

   if (yes_max_evidence_N_is_0 && slen >= 40) {
      // Those are the "super reads".
      const int mm = 1 + tot/2;
      if (aln.score == 0) {
         // Odd number of 30 nt triplets: we place a mm in the
         // event triplets. There are 11 positions on the first
         // and the last triplets, 10 positions on the internal
         // ones. The number of ways to choose the positions is
         // 11^2 * 10^{mm-2}. Those are mm errors (occurrence
         // PROB), the other nucleotides are correct (occurrence
         // 1-PROB). The duplicate has two compensating mutations
         // at those positions (occurrence u/3), the other
         // positions are not mutations (occurrence 1-u). We
         // divide by the probability that the read has score 0,
         // approximately equal to the probability that there is
         // no mutation. We also divide by the probability that
         // p0 is small, approximately 2/3 per segment.
         // Note that we divide 'tot' by 3 to approximate
         // the dependence between consecutive segments.
         // We also count a 1/2 probability that the read is
         // repeated with the specified level of u.
         double u = mm / (double) slen; // Worst value of 'u'.
         return .3 * pow(11*PROB*u/3 / (1-PROB), mm) *
                     pow(1-u,slen-mm) / pow(.666,tot/3.);
         // 0.5 * (11*pu/3)^mm * (1-u)^k-mm * (1-p)^k-mm /
         // (1-p)^k * 2/3^tot/3
      }
      else if (aln.score == 1) {
         // Case 1: the mismatch is in the second (or but-to-last)
         // segment, no further than GAMMA = 17 from the border. An
         // error can occur if the mismatch is an error (frequency p)
         // combined with an uncompensated mutation (2u/3). Other
         // hidden compensated errors are in the even segments.
         // Case 3: the mismatch is somewhere else. The most likely
         // scenario for an error is that the mismatch is a mutation
         // (frequency u) and that it there are two hidden and
         // compensated errors in even semgnents.
         // We need to divide by the probability of occurrence,
         // approximately equal to the probability that the read
         // contains one error at the given position  and that all
         // the segments report N=0 (frequency 2/3 each).
         // Note that we divide 'tot' by 3 to approximate the
         // dependence between consecutive segments.
         // We also count a 1/2 probability that the read is
         // repeated with the specified level of u.
         int errpos = aln.read_beg == 0 ? aln.read_end+1 : aln.read_beg-1;
         int case_1 = (errpos >= 10 && errpos < 17) ||
            (errpos < slen-10 && errpos >= slen-17);
         const int mm = (tot+1)/2;
         if (case_1) {
            double u = mm / (double) slen; // Worst value of 'u'.
            return 0.5 * 2*u/3 * pow(11*PROB*u/3 / (1-PROB), mm-1) *
                  pow(1-u, slen-mm) / pow(.666,tot/3.);
            // 0.5 * 14 * 2*pu/3 * (11*pu/3)^mm-1 * (1-u)^k-mm *
            // (1-p)^k-mm / 14 * p*(1-p)^k-1 * 2/3^tot/3
         }
         else {
            double u = (mm+1) / (double) slen; // Worst value of 'u'.
            return 0.5 * pow(11*PROB*u/3 / (1-PROB), mm) * u*(1-PROB) *
                  pow(1-u,slen-mm-1) / PROB / pow(.666,tot/3.);
            // 0.5 * (k-14) * u*(1-p) * (11*pu/3)^mm * (1-u)^k-mm-1 *
            // (1-p)^k-mm-1 / (k-14) * p*(1-p)^k-1 * 2/3^tot/3
         }
      }
      else if (aln.score == 2) {
         // Assume two uncompensated errors in even segments.
         const int mm = (tot+1)/2;
         double u = mm / (double) slen; // Worst value of 'u'.
         return 2. * mm*(mm-1) * pow(11*PROB*u/3 / (1-PROB), mm) *
               pow((1-PROB)/PROB,2) * pow(1-u,slen-mm) / 
               slen / (slen-1) /pow(.666,tot/3.);
         // 0.5 * (mm choose 2) * (11*2*pu/3)^2 * (10*pu/3)^mm-2 *
         // (1-u)^k-mm * (1-p)^k-mm / k*(k-1)/2*p^2*(1-p)^k-2 * 2/3^tot/3
      }
   }

   // Not a super read: use estimation of N and u.
   qsort(uN0, tot, sizeof(uN0_t), cmpN0);
   int imax = tot-1;
   for ( ; imax > 0 ; imax--) {
      // Take the "worst" N0.
      if (uN0[imax].N0 < worst_case_N0) break;
   }

   int N0 = uN0[imax].N0;
   double u = uN0[imax].u;
   double p0 = 1. / (1. + exp(lev));

   // Probability that the true location is not the best.
   int m = aln.score + 1;
   double pm = exp(lgamma(slen+1)-lgamma(m+1)-lgamma(slen-m+1) +
      m*log(PROB) + (slen-m)*log(1-PROB));
   double pswap = 1 - pow(1 - m*u/3 * pow(1-u,slen-1), N0);
   double ptrue_not_best = pm * pswap;

   ssize_t nerr = (aln.read_beg == 0 || aln.read_end == slen-1) ?
      aln.score-1 :
      aln.score-2;
   ssize_t naln = (aln.read_beg == 0 || aln.read_end == slen-1) ?
      slen - aln.read_end + aln.read_beg - 2 :
      slen - aln.read_end + aln.read_beg - 3;

   if (aln.score == 0) {
      // Not a super read, but since the score is 0, the
      // only way to be wrong is that the true location is
      // not the best. Also condition by the probability
      // that the read has no error.
      double cond = pow(1-PROB, slen);
      double prob = p0 * ptrue_not_best / cond;
      return prob > 1 ? 1 : prob;
   }
   else {
      // Condition by the probability that the read has at
      // the given number of errors.
      double cond = exp(lgamma(slen+1) -lgamma(aln.score+1)
         -lgamma(slen-aln.score+1) +
         aln.score*log(PROB) + (slen-aln.score)*log(1-PROB));
      double poff = auto_mem_seed_offp(slen, u, N0);

      // Count only the evidence from the nucleotides that were
      // actually aligned (those in the seed are already taken
      // into account in the Sesame prior probability). For MEM,
      // count one obligatory error because of the end of the seed.
      // For skip seeds, consider that only one seed was used to
      // discover the locus.

      // Weight of the evidence if mapping is correct...
      double A = nerr * log(PROB) + (naln-nerr) * log(1-PROB);
      // ... if mapping is on a duplicate...
      double PU = PROB + u;
      double B = nerr * log(PU)   + (naln-nerr) * log(1-PU);
      // ... and if mapping is random.
      double C = aln.score * log(.75)  + (naln-aln.score) * log(.25);

      // Here 'term3' uses non-informative prios instead of Sesame
      // priors for the sake of simplicity. In practice, 'term3' is
      // either very close to 0 or very close to 1 and the priors
      // do not matter.
      double term1 = p0 * poff / ( poff + exp(A-B)*(1-poff) );
      double term2 = p0 * ptrue_not_best / cond;
      double term3 = aln.score < slen / 5 ? 0 : 1. / (1. + exp(A-C));

      double max = term1 > term2 ? term1 : term2;
      return max + term3 > 1. ? 1. : max + term3;

   }

}


void
batchmap
(
   const char * indexfname,
   const char * readsfname
)
{

   fprintf(stderr, "loading index... ");
   // Load index files.
   index_t idx = load_index(indexfname);

   fprintf(stderr, "done.\n");
   FILE * inputf = fopen(readsfname, "r");
   if (inputf == NULL) exit_cannot_open(readsfname);

#ifdef DEBUG
   fprintf(stdout, "smmfdp: index=%s, reads=%s, perror=%f, skip-thr=%f\n", indexfname, readsfname, PROB, SKIPQUAL);
#endif
   
   size_t sz = 64;
   ssize_t rlen;
   char * seq = malloc(64);
   exit_error(seq == NULL);

   size_t counter = 0; // Used for randomizing.
   size_t maxlen = 0; // Max 'k' value for seeding probabilities.

   // Read sequence file line by line.
   while ((rlen = getline(&seq, &sz, inputf)) != -1) {
      
      // Remove newline character if present.
      if (seq[rlen-1] == '\n') seq[rlen-1] = '\0';

      // If fasta header, print sequence name.

// XXX BENCHMARK STUFF XXX //
#ifdef FASTOUT
      if (seq[0] == '>') {
         fprintf(stdout, ">%s\n", seq+1);
         continue;
      }
#else
      if (seq[0] == '>') {
         fprintf(stdout, "%s\t", seq+1);
         continue;
      }
#endif

      if (rlen > maxlen) {
         maxlen = rlen;
         // (Re)initialize library.
         sesame_set_static_params(GAMMA, rlen, PROB);
      }

      // Compute L1, L2 and MEMs.
      seed_t l1, l2;
      extend_l1l2(seq, idx, &l1, &l2);

      // Compute seeds.
      wstack_t * seeds = mem_seeds(seq, idx, GAMMA);

      // Return if no seeds were found
      if (seeds->pos == 0) {
	 // Did not find anything.
         free(seeds);
         fprintf(stdout, "%s\tNA\tNA\n", seq);
         continue;
      }

      // Compute N(L1,L2)
      N = new_estimate(l1,l2); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

      // Quick mode: only align longest MEMs
      if (N > QUICK_DUPLICATES) 
	 filter_longest_mem(mems);

      alnstack_t * alnstack = mapread(mems, seq, idx, rlen);

      // Did not find anything.
      if (alnstack->pos == 0) {
         free(alnstack);
         fprintf(stdout, "%s\tNA\tNA\n", seq);
         continue;
      }

      // Pick a top alignment at "random".
      aln_t a = alnstack->aln[counter++ % alnstack->pos];

      int there_is_only_one_best_hit = 1;
      if (alnstack->pos > 1) {
         // See if the hits are actually distinct.
         size_t ref = alnstack->aln[0].refpos;
         for (int i = 1 ; i < alnstack->pos ; i++) {
            if (alnstack->aln[i].refpos != ref) {
               there_is_only_one_best_hit = 0;
               break;
            }
         }
      }

      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      if (there_is_only_one_best_hit)
	 a.qual = N > QUICK_DUPLICATES ? quality_quick(a, seq, idx) : quality(a, seq, idx); 
      else
	 a.qual = 1-1./alnstack->pos;
	 
      // Report mapping results
      char * apos = chr_string(a.refpos, idx.chr);
      fprintf(stdout, "%s\t%s\t%e\n", seq, apos, a.qual);
      free(apos);

      // Free alignments
      for(size_t i = 0; i < alnstack->pos; i++)
         free(alnstack->aln[i].refseq);
      free(alnstack);

   }

   fclose(inputf);
   free_index_chr(idx.chr);
   free(seq);
   sesame_clean();

   // FIXME: Would be nice to do this, but the size is unknown.
   // FIXME: Either forget it, or store the size somewhere.
   //munmap(idx.csa);
   //munmap(idx.bwt);
   //munmap(idx.occ);
   //munmap(idx.lut);

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
      fprintf(stdout, "%s", HELP_MSG);
   } else if (strcmp(argv[1], "--index") == 0) {
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
      
      batchmap(index_path, reads_path);
   }
}
