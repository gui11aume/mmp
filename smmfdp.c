#define _GNU_SOURCE
#include "divsufsort.h"
#include "bwt.h"


// Error-handling macro.
#define exit_cannot_open(x) \
   do { fprintf(stderr, "cannot open file '%s' %s:%d:%s()\n", (x), \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); } while(0)

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
   fprintf(stderr, "done\n");

   fprintf(stderr, "creating suffix array... ");
   int64_t * sa = compute_sa(genome);
   fprintf(stderr, "done\n");

   fprintf(stderr, "creating BWT... ");
   bwt_t * bwt = create_bwt(genome, sa);
   fprintf(stderr, "done\n");

   fprintf(stderr, "creating Occ table... ");
   occ_t * occ = create_occ(bwt);
   fprintf(stderr, "done\n");

   fprintf(stderr, "filling lookup table... ");
   lut_t * lut = malloc(sizeof(lut_t));
   fill_lut(lut, occ, (range_t) {.bot=1, .top=strlen(genome)}, 0, 0);
   fprintf(stderr, "done\n");

   fprintf(stderr, "compressing suffix array... ");
   csa_t * csa = compress_sa(sa);
   fprintf(stderr, "done\n");

   // Write files
   char buff[256];
   char * data;
   ssize_t ws;
   size_t sz;

   // Write the suffix array file.
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


   // Write the Burrows-Wheeler transform.
   sprintf(buff, "%s.occ", fname);
   int focc = creat(buff, 0644);
   if (focc < 0) exit_cannot_open(buff);

   ws = 0;
   sz = sizeof(occ_t) + occ->nrows * SIGMA * sizeof(blocc_t);
   data = (char *) occ;
   while (ws < sz) ws += write(focc, data + ws, sz - ws);
   close(focc);

   // Clean up.
   free(csa);
   free(bwt);
   free(occ);
   free(lut);

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

}
