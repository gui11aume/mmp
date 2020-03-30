#define VERSION "1.0"
#define _GNU_SOURCE

#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>

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
typedef struct txtbuf_t txtbuf_t;
typedef struct read_t read_t;
typedef struct writerarg_t writerarg_t;
typedef struct batch_t batch_t;

struct uN0_t {
   double u;
   size_t N0;
   double p0;
};

struct txtbuf_t {
   size_t size;
   size_t pos;
   char   txt[];
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
   pthread_cond_t  * cond;
   pthread_cond_t  * cond_reader;
};

enum status_t {idle = 0, map = 1, output = 2};

struct batch_t {
   // Data
   index_t           idx;
   size_t            lineid;
   size_t            n_reads;
   txtbuf_t        * readbuf;
   txtbuf_t        * outbuf;
   // Thread signaling
   enum status_t     status;
   int             * act_threads;
   pthread_mutex_t * mutex_sesame;
   pthread_mutex_t * mutex;
   pthread_cond_t  * cond_reader;
   pthread_cond_t  * cond_writer;
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

txtbuf_t *
txtbuf_new
(
 size_t size
)
{
   txtbuf_t * buf = malloc(sizeof(txtbuf_t) + size);
   exit_error(buf == NULL);

   buf->pos = 0;
   buf->size = size;

   return buf;
}

void
buf_pushtxt
(
 char      * txt,
 txtbuf_t ** bufp
)
{
   size_t slen = strlen(txt);
   txtbuf_t * buf = *bufp;

   if (buf->pos + slen >= buf->size) {
      size_t newsize = buf->size + slen + 1;
      *bufp = buf = realloc(buf, sizeof(txtbuf_t) + newsize);
      exit_error(buf == NULL);
      buf->size = newsize;
   }

   memcpy(buf->txt + buf->pos, txt, slen+1);
   buf->pos += slen;
}

void
buf_reset
(
 txtbuf_t * buf
)
{
   buf->pos = 0;
   if (buf->size > 0)
      *buf->txt = 0;
}

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
   // Probability of the event if there is one duplicate.
   double prob_1 = mu * pow(1-mu, 2*m) * pow(1-pow(.25, m), G) +
      (1-pow(1-mu,m)) * pow(1-mu,m) * (pow(1-pow(.25,m+1), G) - pow(1-pow(.25,m), G)) +
      mu * pow(1-mu, 2*m) * (pow(1-pow(.25,m+1), G) - pow(1-pow(.25,m), G)) ;
   // Probability of the event if one end is duplicated.
   double prob_2 = pow(1-mu,m) *
      (pow(1-pow(.25,m+1), G) - pow(1-pow(.25,m), G));
   // Probability of the event if no end is duplicated.
   double prob_3 = (1.-pow(1-pow(.25,m+1), G)) *
      (pow(1-pow(.25,m+1), G) - pow(1-pow(.25,m), G));
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


read_t
next_read
(
        char   ** txt,
   enum fmt_t     format
)
{
   read_t read;
   if (format == fasta) {
      read.phred = "*";
      // Name
      read.name = *txt;
      // Mark Name string end
      read.seq  = strchr(read.name, '\n');
      *read.seq = 0;
      read.seq++;
      // Cut long read names
      strtok(read.name, "\t ");
      // Mark Seq string end and update txt pointer
      *txt = strchr(read.seq, '\n');
      if (*txt) {
	 // Last line may not end with '\n'
	 **txt = 0;
	 *txt += 1;
      } else {
	 *txt = strchr(read.seq, '\0') + 1;
      }
   } else {
      // Name
      read.name = *txt;
      // Mark Name string end
      read.seq  = strchr(read.name, '\n');
      *read.seq = 0;
      read.seq++;
      // Cut long read names
      strtok(read.name, "\t ");
      // Mark Seq string end
      read.phred = strchr(read.seq, '\n');
      *read.phred = 0;
      read.phred = strchr(read.phred+1, '\n') + 1;
      // Mark Phred string end and update txt pointer
      *txt = strchr(read.phred, '\n');
      if (*txt) {
	 // Last line may not end with '\n'
	 **txt = 0;
	 *txt += 1;
      } else {
	 *txt = strchr(read.phred, '\0') + 1;
      }
   }
   return read;
}

void *
batch_map
(
 void * arg
)
{
   batch_t * batch  = (batch_t *)arg;
   index_t idx = batch->idx;
   char  * txt = batch->readbuf->txt;

   // Check read format
   enum fmt_t format = unset;
   if (txt[0] == '>') format = fasta;
   else format = fastq;

   // Output line buffer
   size_t samsz   = 1024;
   char * samline = malloc(samsz);
   exit_error(samline == NULL);
   
   read_t read;
   for (size_t i = 0; i < batch->n_reads; i++) {
      read = next_read(&txt, format);
      size_t rlen = strlen(read.seq);
      
      // Compute L1, L2 and MEMs.
      seed_t L1, L2;
      extend_L1L2(read.seq, rlen, idx, &L1, &L2);

      // Compute seeds.
      wstack_t * seeds = mem_seeds(read.seq, idx, GAMMA);

      // Return if no seeds were found
      if (seeds->pos == 0) {
        // Did not find anything.
        free(seeds);
        // Output in sam format.
	int olen = snprintf(NULL, 0, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",
			    read.name+1, read.seq, read.phred);
	if (olen+1 > samsz) {
	   samsz = olen+1;
	   samline = realloc(samline, samsz);
	   exit_error(samline == NULL);
	}
        sprintf(samline, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",
		read.name+1, read.seq, read.phred);
	buf_pushtxt(samline, &(batch->outbuf));
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

      alnstack_t * alnstack = mapread(seeds, read.seq, idx, rlen, batch->lineid + i);

      // Did not find anything.
      if (alnstack->pos == 0) {
	int olen = snprintf(NULL, 0, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",
			    read.name+1, read.seq, read.phred);
	if (olen+1 > samsz) {
	   samsz = olen+1;
	   samline = realloc(samline, samsz);
	   exit_error(samline == NULL);
	}
        sprintf(samline, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n",
		read.name+1, read.seq, read.phred);
	buf_pushtxt(samline, &(batch->outbuf));
        free(alnstack);
        // Output in sam format.
        // Free seeds
        for (size_t i = 0; i < seeds->pos; i++) {
          seed_t * s = (seed_t *) seeds->ptr[i];
          free(s->sa);
          free(s);
        }
        free(seeds);
        continue;
      }

      // Pick a top alignment at "random".
      aln_t a = alnstack->aln[(batch->lineid + i) % alnstack->pos];

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

      if (there_is_only_one_best_hit) {
      a.qual = uN0.N0 > QUICK_DUPLICATES ?
        quality_low(rlen, longest_mem, uN0) :
	 quality(a, read.seq, idx, uN0, batch->mutex_sesame);
      }
      else {
        a.qual = 1-1./alnstack->pos;
      }
	 
      // Report mapping results
      pos_t pos = get_pos(a.refpos, idx.chr);
      // Output in sam format.
      int bits = pos.strand ? 0 : 16;
      size_t leftpos = pos.strand ? pos.pos : pos.pos - rlen+1;
      int olen = snprintf(NULL, 0, "%s\t%d\t%s\t%ld\t%d\t%ldM\t*\t0\t0\t%s\t%s\tXS:i:%d\n",
			  read.name+1, bits, pos.rname, leftpos, (int) (-10*log10(a.qual)),
			  rlen, read.seq, read.phred, a.score);
      if (olen+1 > samsz) {
	 samsz = olen+1;
	 samline = realloc(samline, samsz);
	 exit_error(samline == NULL);
      }
      sprintf(samline, "%s\t%d\t%s\t%ld\t%d\t%ldM\t*\t0\t0\t%s\t%s\tXS:i:%d\n",
         read.name+1, bits, pos.rname, leftpos, (int) (-10*log10(a.qual)),
         rlen, read.seq, read.phred, a.score);
      buf_pushtxt(samline, &(batch->outbuf));
      
      // Free seeds
      for (size_t i = 0; i < seeds->pos; i++) {
	 seed_t * s = (seed_t *) seeds->ptr[i];
	 free(s->sa);
	 free(s);
      }
      free(seeds);

      // Free alignments
      for(size_t i = 0; i < alnstack->pos; i++)
         free(alnstack->aln[i].refseq);
      free(alnstack);
   }

   pthread_mutex_lock(batch->mutex);
   // Decrease active threads
   *(batch->act_threads) = *(batch->act_threads) - 1;
   pthread_cond_signal(batch->cond_reader);

   // Signal writer monitor
   batch->status = output;
   pthread_cond_signal(batch->cond_writer);
   pthread_mutex_unlock(batch->mutex);

   // Free variables
   free(samline);

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
	 // TODO:
	 // - set in_use to 0
	 // - signal reader
	 
	 // Write output buffer
	 fprintf(stdout, "%s", batch->outbuf->txt);


	 pthread_mutex_lock(warg->mutex);
	 // Free batch and signal reader
	 batch->status = idle;
	 pthread_cond_signal(warg->cond_reader);
	 // Wait until next batch is assigned
	 while (batch->next_batch == (void *)-1)
	    pthread_cond_wait(warg->cond, warg->mutex);
	 pthread_mutex_unlock(warg->mutex);
	 
	 batch = batch->next_batch;
      }

      if (!batch) break;

      pthread_mutex_lock(warg->mutex);
      pthread_cond_wait(warg->cond, warg->mutex);
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
   ssize_t maxlen = 0; // Max 'k' value for seeding probabilities.
   
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
   else if (first == '@') format = fastq;
   else {
      fprintf(stderr, "error: unrecognized input format (must be fasta or fastq)\n");
      exit(1);
   }
   ungetc(first, inputf);
   print_sam_header(idx);

   int lines_read;
   if (format == fasta) lines_read = 2;
   else lines_read = 4;

   // Multithreading variables
   int BATCHSIZE = 10000;
   int act_threads = 0;

   // Read buffer variables
   int b = -1, b_it = 0, success = 0;
   size_t sz = 256, line = 0;
   char   * lineptr = malloc(sz);

   // Create mutex and monitors
   pthread_mutex_t mutex_sesame = PTHREAD_MUTEX_INITIALIZER;
   pthread_mutex_t mutex        = PTHREAD_MUTEX_INITIALIZER;
   pthread_cond_t  cond_reader  = PTHREAD_COND_INITIALIZER, cond_writer = PTHREAD_COND_INITIALIZER;

   // Create batches
   batch_t * last_batch = NULL;
   batch_t ** batch = malloc(2*MAXTHREADS * sizeof(void*));
   exit_error(batch == NULL);
   for (int i = 0; i < 2*MAXTHREADS; i++) {
      batch[i] = malloc(sizeof(batch_t));
      exit_error(batch[i] == NULL);
      batch[i]->idx          = idx;
      batch[i]->readbuf      = txtbuf_new(100*BATCHSIZE);
      batch[i]->outbuf       = txtbuf_new(100*BATCHSIZE);
      batch[i]->status       = idle;
      batch[i]->act_threads  = &act_threads;
      batch[i]->mutex_sesame = &mutex_sesame;
      batch[i]->mutex        = &mutex;
      batch[i]->cond_reader  = &cond_reader;
      batch[i]->cond_writer  = &cond_writer;
      batch[i]->next_batch   = (void *)-1;
   }

   // Create writer thread
   writerarg_t warg = (writerarg_t){batch[0], &mutex, &cond_writer, &cond_reader};
   pthread_t writer_thread;
   pthread_create(&writer_thread, NULL, mt_write, &warg);

   while (1) {
      // Select next buffer
      pthread_mutex_lock(&mutex);
      do {
	 // Find free buffer
	 for (b_it = 0; b_it < 2*MAXTHREADS; b_it++) {
	    b = (b+1)%(2*MAXTHREADS);
	    if (batch[b]->status == idle) {
	       buf_reset(batch[b]->readbuf);
	       buf_reset(batch[b]->outbuf);
	       break;
	    }
	 }
	 
	 // Wait for free buffer (writer may be busy)
	 if (b_it == 2*MAXTHREADS)
	    pthread_cond_wait(&cond_reader, &mutex);
	 
      } while (b_it == 2*MAXTHREADS);
      
      if (last_batch) {
	 // Complete linked list and signal writer
	 last_batch->next_batch = batch[b];
	 pthread_cond_signal(&cond_writer);
      }
      pthread_mutex_unlock(&mutex);

      last_batch = batch[b];
      
      // Read input batch
      size_t n_lines = 0;
      size_t line_id = line;
      while ((success = getline(&lineptr, &sz, inputf)) > 0) {
	 // Set sesame static params
	 if (line%lines_read == 1 && strlen(lineptr)-1 > maxlen) {
	    maxlen = strlen(lineptr)-(lineptr[strlen(lineptr)-1] == '\n');
	    pthread_mutex_lock(&mutex_sesame);
	    sesame_set_static_params(GAMMA, maxlen, PROB);
	    pthread_mutex_unlock(&mutex_sesame);
	 }

	 // Push line to buffer
	 buf_pushtxt(lineptr, &(batch[b]->readbuf));

	 line++; n_lines++;

	 if (n_lines >= BATCHSIZE*lines_read)
	    break;
      }

      if (success < 0 && (errno == EINVAL || errno == ENOMEM)) {
	 fprintf(stderr, "error while reading input file (line %ld)\n", line+1);
	 exit(EXIT_FAILURE);
      }

      // Fill new batch info
      batch[b]->next_batch = success < 1 ? NULL: (void *)-1;
      batch[b]->lineid     = line_id/lines_read;
      batch[b]->n_reads    = n_lines/lines_read;
      batch[b]->status     = map;
      
      // Create new thread
      pthread_mutex_lock(&mutex);
      act_threads++;
      pthread_mutex_unlock(&mutex);
      
      pthread_t thread;
      if (pthread_create(&thread, NULL, batch_map, batch[b])) exit_error(1);
      pthread_detach(thread);

      // Break on last batch
      if (success < 1) break;

      // Wait for condition signal (available threads)
      pthread_mutex_lock(&mutex);
      while (act_threads >= MAXTHREADS)
	 pthread_cond_wait(&cond_reader, &mutex);
      pthread_mutex_unlock(&mutex);
   }

   // Wait for writer.
   pthread_join(writer_thread, NULL);

   // Free stuff
   free(lineptr);
   for (int i = 0; i < 2*MAXTHREADS; i++) {
      free(batch[i]->readbuf);
      free(batch[i]->outbuf);
      free(batch[i]);
   }
   free(batch);
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
