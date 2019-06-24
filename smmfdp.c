#define _GNU_SOURCE

#include <math.h>

#include "bwt.h"
#include "map.h"
#include "divsufsort.h"
#include "sesame.h"

#define GAMMA 17
#define PROBDEFAULT 0.01
#define SKIPQUALDEFAULT 10

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
         index_t idx,
         int     skip
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
   for (int s = 0 ; s <= slen-30 ; s += 10, tot++) {
      if (aln.read_end < s+10 || aln.read_beg > s+19) {
         // Estimate only at the site of the seed.
         uN0[tot] = (uN0_t) { worst_case_mu, worst_case_N0, 0 };
         continue;
      }
      uN0[tot] = estimate_uN0(aln.refseq + s, idx);
      lev += uN0[tot].lev;
      if (uN0[tot].u != uN0[tot].u) {
         // In case estimation fails (e.g., because of "N").
         uN0[tot] = (uN0_t) { worst_case_mu, worst_case_N0, 0 };
      }
   }

   // Pick the largest value of N0.
   qsort(uN0, tot, sizeof(uN0_t), cmpN0);
   int imax = tot-1;
   for ( ; imax > 0 ; imax--) {
      if (uN0[imax].N0 < worst_case_N0) break;
   }

   int N0 = uN0[imax].N0;
   double u = uN0[imax].u;

   double p0 = 1. / (1. + exp(lev));

   // Probability that the true location is not the best.
   int m = aln.score + 1;
   double pm = exp(lgamma(slen+1)-lgamma(m+1)-lgamma(slen-m+1) +
      m*log(PROB) + (slen-m)*log(1-PROB));
   double pswap = 1 - pow(1 - pow(u/3,m) * pow(1-u,slen-m), N0);
   double ptrue_not_best = pm * pswap;

   ssize_t nerr = skip == 0 ?
      ((aln.read_beg == 0 || aln.read_end == slen-1) ?
         aln.score-1 : aln.score-2) : aln.score;
   ssize_t naln = skip == 0 ?
      ((aln.read_beg == 0 || aln.read_end == slen-1) ?
          slen - aln.read_end + aln.read_beg - 2 :
          slen - aln.read_end + aln.read_beg - 3) : slen-GAMMA;

   if (aln.score == 0 || (aln.score == 1 && naln == 0)) {
      // Perfect score or error on the flank. Also condition
      // on the probability that the read has no error or one
      // error on the flank.
      if (p0 < .5) {
         // Special correction for error on the flank.
         double x = aln.score == 0 ? 1 : 2*PROB;
         // Those are the "super reads". Separate odd vs even.
         // We estimate that 0.3 is the fraction of sites with
         // exactly 1 duplicate.
         if (tot % 2 == 1) {
            return pow(11*PROB,tot-1) / pow(1-PROB,tot-1) *
               pow(u/3,tot-1)*pow(1-u,slen-tot+1) * .3 / x;
         }
         else {
            return 2*11*pow(PROB,tot-2) / pow(1-PROB,tot-2) *
               pow(u/3,tot-2)*pow(1-u,slen-tot+2) * .3 / x;
         }
      }
   }
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
      double poff = skip == 0 ?
         auto_mem_seed_offp(slen, u, N0) :
         auto_skip_seed_offp(slen, skip, u, N0);
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

      return term1 + term2 + term3 > 1. ? 1. : term1 + term2 + term3;

   }


}

aln_t
mapfrag
(
   const char    * seq,
         index_t   idx
)
// Processes read fragments of length at most 50.
{

   int skip = 0; // Try MEM seeding first.
   int best_score = 50;
   aln_t aln[2]  = {{0}};

   aln[0].score = aln[1].score = 9999;

   for (int redo = 0 ; redo < 2 ; redo++) {

      alnstack_t * alnstack = mapread(seq, idx, GAMMA, skip, best_score);
      if (!alnstack) exit(EXIT_FAILURE);

      // Did not find anything.
      if (alnstack->pos == 0) {
         free(alnstack);
         break;
      }

      // Pick first top alignment.
      aln_t a = alnstack->aln[0];
      aln[redo].score    = a.score;
      aln[redo].refpos   = a.refpos;
      aln[redo].refseq   = a.refseq;
      aln[redo].read_beg = a.read_beg;
      aln[redo].read_end = a.read_end;

      best_score = a.score;

      aln[redo].qual = alnstack->pos > 1 ?
         1.0 - 1.0 / alnstack->pos :
         quality(aln[redo], seq, idx, skip);

      // Free alignments
      for(size_t i = 0; i < alnstack->pos; i++)
         free(alnstack->aln[i].refseq);
      free(alnstack);

      // We are done if we found a perfect hit or
      // if the quality is higher than 'SKIPQUAL'.
      // (-4.34 = -10/log(10) amounts to taking log10).
      if (aln[0].score == 0 || -4.34*log(aln[0].qual) >= SKIPQUAL)
         break;

      // Otherwise try skip seeds
      skip = 12;

   }

   // Use MEM alignment if better than skip seed alignment,
   // or if equally good and has discovered several equally
   // good hits in the genome.
   int MEM_is_best = aln[0].score < aln[1].score;
   int no_improvement = aln[0].score == aln[1].score && aln[0].qual < .5;

   return (MEM_is_best || no_improvement) ? aln[0] : aln[1];

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
   fprintf(stdout, "smmfdp: index=%s, reads=%s, perror=%f, skip-thr=%f\n",
         indexfname, readsfname, PROB, SKIPQUAL);
#endif
   
   alnstack_t *best = alnstack_new(10);

   size_t sz = 64;
   ssize_t rlen;
   char * seq = malloc(64);
   exit_error(seq == NULL);

   size_t maxlen = 0; // Max 'k' value for seeding probabilities.

   // Read sequence file line by line.
   while ((rlen = getline(&seq, &sz, inputf)) != -1) {
      
      // Remove newline character if present.
      if (seq[rlen-1] == '\n') seq[rlen-1] = '\0';

      // If fasta header, print sequence name.

      if (seq[0] == '>') {
         fprintf(stdout, "%s\t", seq+1);
         continue;
      }

      rlen = strlen(seq);
      if (rlen > maxlen) {
         maxlen = rlen;
         // (Re)initialize library.
         sesame_set_static_params(GAMMA, rlen, PROB);
      }

      if (rlen <= 50) {
         aln_t aln = mapfrag(seq, idx);

         if (aln.score >= 9999) {
            // Nothing found.
            fprintf(stdout, "%s\tNA\tNA\n", seq);
         }
         else {
            char *apos = chr_string(aln.refpos, idx.chr);
            fprintf(stdout, "%s\t%s\t%e\n", seq, apos, aln.qual);
            free(apos);
         }
      }
      else {
         // Cut the read in fragments.
         int nfrags = 1 + (rlen-1) / 50;
         double span = rlen / 50.0 - 1.0;
         aln_t alns[5]= {0};
         char frag[51] = {0};
         for (int i = 0 ; i < nfrags ; i++) {
            int shift = span * 50 * i;
            strncpy(frag, seq + shift, 50);
            alns[i] = mapfrag(frag, idx);
         }

         // Pick the best.
         int iopt = 0;
         aln_t aln = alns[iopt];
         for (int i = 1 ; i < nfrags ; i++) {
            if (alns[i].qual < aln.qual) {
               aln = alns[i];
               iopt = i;
            }
         }

         if (aln.score >= 9999) {
            // Nothing found.
            fprintf(stdout, "%s\tNA\tNA\n", seq);
            continue;
         }

         // Find agreement among fragments.
         for (int i = 0 ; i < nfrags ; i++) {
            if (i == iopt) continue;
            if (abs(alns[i].refpos - aln.refpos) < rlen) {
               aln.qual *= alns[i].qual;
            }
         }

         // Output.
         char *apos = chr_string(aln.refpos, idx.chr);
         fprintf(stdout, "%s\t%s\t%e\n", seq, apos, aln.qual);
         free(apos);

#if 0
// This code can be used if one wants to know the alignment
// score of the read at the candidate location.

         // Re-align to get new score.
         seed_t seed = {
            .beg = aln.read_beg,
            .end = aln.read_end,
         };
         align_t a = {
            .refpos = aln.refpos - 50*iopt,
            .span = rlen,
            .minscore = 0,
            .seed = &seed,
         };

         int best_score = 9999;
         align(a, seq, idx.dna, idx.occ->txtlen, &best_score, &best);

         char *apos = chr_string(a.refpos, idx.chr);
         fprintf(stdout, "%s\t%s\t%e\n", seq, apos, aln.qual);
         free(apos);
#endif

      }

   }

   free(best);

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
