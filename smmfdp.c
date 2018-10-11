#define _GNU_SOURCE

#include <math.h>

#include "bwt.h"
#include "divsufsort.h"
#include "mem_seed_prob.h"

#define GAMMA 17

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


void
build_index
(
   const char * fname
)
{

   // Open fasta file.
   FILE * fasta = fopen(fname, "r");
   if (fasta == NULL) exit_cannot_open(fname);

   // Read and normalize genome
   fprintf(stderr, "reading genome... ");
   char * genome = normalize_genome(fasta);
   fprintf(stderr, "done.\n");

   fprintf(stderr, "creating suffix array... ");
   int64_t * sa = compute_sa(genome);
   fprintf(stderr, "done.\n");

   fprintf(stderr, "creating bwt... ");
   bwt_t * bwt = create_bwt(genome, sa);
   fprintf(stderr, "done.\n");

   fprintf(stderr, "creating Occ table... ");
   occ_t * occ = create_occ(bwt);
   fprintf(stderr, "done.\n");

   fprintf(stderr, "filling lookup table... ");
   lut_t * lut = malloc(sizeof(lut_t));
   fill_lut(lut, occ, (range_t) {.bot=1, .top=strlen(genome)}, 0, 0);
   fprintf(stderr, "done.\n");

   fprintf(stderr, "compressing suffix array... ");
   csa_t * csa = compress_sa(sa);
   fprintf(stderr, "done.\n");

   // Write files
   char buff[256];
   char * data;
   ssize_t ws;
   size_t sz;


   // Write the compressed suffix array file.
   sprintf(buff, "%s.sa", fname);
   int fsar = creat(buff, 0644);
   if (fsar < 0) exit_cannot_open(buff);
   
   ws = 0;
   sz = sizeof(csa_t) + csa->nint64 * sizeof(int64_t);
   data = (char *) csa;
   while (ws < sz) ws += write(fsar, data + ws, sz - ws);
   close(fsar);


   // Write the Burrows-Wheeler transform.
   sprintf(buff, "%s.bwt", fname);
   int fbwt = creat(buff, 0644);
   if (fbwt < 0) exit_cannot_open(buff);

   ws = 0;
   sz = sizeof(bwt_t) + bwt->nslots * sizeof(uint8_t);
   data = (char *) bwt;
   while (ws < sz) ws += write(fbwt, data + ws, sz - ws);
   close(fbwt);


   // Write the Occ table.
   sprintf(buff, "%s.occ", fname);
   int focc = creat(buff, 0644);
   if (focc < 0) exit_cannot_open(buff);

   ws = 0;
   sz = sizeof(occ_t) + occ->nrows * SIGMA * sizeof(blocc_t);
   data = (char *) occ;
   while (ws < sz) ws += write(focc, data + ws, sz - ws);
   close(focc);


   // Write the lookup table
   sprintf(buff, "%s.lut", fname);
   int flut = creat(buff, 0644);
   if (flut < 0) exit_cannot_open(buff);

   ws = 0;
   sz = sizeof(lut_t);
   data = (char *) lut;
   while (ws < sz) ws += write(flut, data + ws, sz - ws);
   close(flut);

   // Clean up.
   free(csa);
   free(bwt);
   free(occ);
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

   return (index_t) { .csa = csa, .bwt = bwt, .occ = occ, .lut = lut };

}


void
analyze_mem
(
   mem_t mem
)
{

   const int n = GAMMA;
   const double MU[6] = {.01, .02, .05, .10, .15, .25};

   const double L1 = mem.left[n] - 1;
   const double R1 = mem.right[n] - 1;
         double L2 = (n+1) * (mem.left[n] - 1);
         double R2 = (n+1) * (mem.right[n] - 1);

   // FIXME: put a sentinel in 'mem.left'.
   for (int i = GAMMA+1 ; mem.left[i] > 0 ; i++) {
      L2 += mem.left[i] - 1;
      R2 += mem.right[i] - 1;
   }

   double loglik = -INFINITY;
   size_t best_N0 = 0.0;
   double best_mu = 0.0;

   for (int iter = 0 ; iter < 6 ; iter++) {

      double mu = MU[iter];

      double L3 = L2 / (1-mu) - L1 / mu;
      double R3 = R2 / (1-mu) - R1 / mu;

      double C = (1 - pow(1-mu,n)) / (n * pow(1-mu,n-1));

      // Compute the number of duplicates.
      double N0 = (L1+R1 + C*(L3+R3)) / 2.0;

      if (N0 < 1.0) N0 = 1.0;

      double tmp = 2*lgamma(N0+1) + (L1+R1) * log(mu)
         + (L2+R2) * log(1-mu)
         + (2*N0 - (L1+L2)) * log(1-pow(1-mu,n))
         - lgamma(N0-L1+1)  - lgamma(N0-L2+1);

      if (tmp < loglik) {
         break;
      }

      loglik = tmp;
      best_N0 = ceil(N0);
      best_mu = mu;

   }

   fprintf(stderr, "N0 = %ld, mu = %.2f\n", best_N0, best_mu);
   double prob = mem_false_pos(50, best_mu, best_N0);
   fprintf(stderr, "False prob: %f\n", prob);

   return;

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

   // Load the genome.
   FILE * fasta = fopen(indexfname, "r");
   if (fasta == NULL) exit_cannot_open(indexfname);
   char * genome = normalize_genome(fasta);

   fprintf(stderr, "done.\n");
   FILE * inputf = fopen(readsfname, "r");
   if (inputf == NULL) exit_cannot_open(readsfname);

   // Initialize MEM seed probablities.
   set_params_mem_prob(GAMMA, 50, .01);

   size_t sz = 64;
   ssize_t rlen;
   char * seq = malloc(64);
   exit_error(seq == NULL);

   // Read sequence file line by line.
   while ((rlen = getline(&seq, &sz, inputf)) != -1) {
      if (seq[rlen-1] == '\n') seq[rlen-1] = '\0';
      alnstack_t * aln = mapread(seq, idx, genome, GAMMA);

      // VERBOSE (DEBUG).
      fprintf(stderr, "Best alignment(s):\n");

      for (int i = 0; i < aln->pos; i++) {
         aln_t a = aln->aln[i];
         fprintf(stderr, "[%d]\n  score: %d\n  pos: %ld\n  MEMs:\n",
               i, a.score, a.refpos);
         for (int j = 0; j < a.nmem; j++) {

            mem_t mem = a.mem[j];
            fprintf(stderr, "   [%d] beg: %ld, end: %ld\n",
                  j, mem.beg, mem.end);

            analyze_mem(mem);

            //fprintf(stderr, "   Left : ");
            //for (int r = GAMMA ; mem.left[r] > 0 ; r++) {
            //   fprintf(stderr, "%zu ", mem.left[r]);
            //}
            //fprintf(stderr, "\n   Right: ");
            //for (int r = GAMMA ; mem.right[r] > 0 ; r++) {
            //   fprintf(stderr, "%zu ", mem.right[r]);
            //}
            //fprintf(stderr, "\n");
         }
      }
   }

   free(seq);
   clean_mem_prob();

   // FIXME: Would be nice to do this, but the size is unknown.
   // FIXME: Either forget it, or store the size somewhere.
   //munmap(idx.csa);
   //munmap(idx.bwt);
   //munmap(idx.occ);
   //munmap(idx.lut);

}


int main(int argc, char ** argv) {

   // Sanity checks.
   if (argc < 2) {
      fprintf(stderr, "First argument must be \"index\" or \"mem\".\n");
      exit(EXIT_FAILURE);
   }

   if (strcmp(argv[1], "index") == 0) {
      if (argc < 3) {
         fprintf(stderr, "Specify file to index.\n");
         exit(EXIT_FAILURE);
      }
      build_index(argv[2]);
   }
   else if (strcmp(argv[1], "mem") == 0) {
      if (argc < 4) {
         fprintf(stderr, "Specify index and read file.\n");
         exit(EXIT_FAILURE);
      }
      batchmap(argv[2], argv[3]);
   }
   else {
      fprintf(stderr, "First argument must be \"index\" or \"mem\".\n");
      exit(EXIT_FAILURE);
   }

}
