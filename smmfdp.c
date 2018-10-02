#define _GNU_SOURCE
#include "divsufsort.h"
#include "bwt.h"

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
   fprintf(stderr, "done\n");

   fprintf(stderr, "creating suffix array... ");
   int64_t * sa = compute_sa(genome);
   fprintf(stderr, "done\n");

   fprintf(stderr, "creating bwt... ");
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
batchmap
(
   const char * indexfname,
   const char * readsfname
)
{

   // Load index files.
   index_t idx = load_index(indexfname);

   // Load the genome.
   FILE * fasta = fopen(indexfname, "r");
   if (fasta == NULL) exit_cannot_open(indexfname);
   char * genome = normalize_genome(fasta);

   FILE * inputf = fopen(readsfname, "r");
   if (inputf == NULL) exit_cannot_open(readsfname);

   size_t sz = 64;
   ssize_t rlen;
   char * seq = malloc(64);
   exit_error(seq == NULL);

   // Read sequence file line by line.
   while ((rlen = getline(&seq, &sz, inputf)) != -1) {
      if (seq[rlen-1] == '\n') seq[rlen-1] = '\0';
      mapread(seq, idx, genome, GAMMA);
   }

   free(seq);

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
