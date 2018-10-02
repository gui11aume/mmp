#define _GNU_SOURCE
#include <assert.h>
#include "divsufsort.h"
#include "bwt.h"


// FIXME: The logic of the 'mem_t' struct is very weak and does not
// FIXME: stand a chance on something else than the test data set.
typedef struct mem_t mem_t;

#define LEN 50

struct mem_t {
   size_t  beg;
   size_t  end;
   range_t range;
   // Decay cascades.
   size_t  fwd[50];
   size_t  bwd[50];
};

// XXX: Here is a global array to avoid 'malloc()', constructors etc.
// XXX: This is a no-go for a real world implementation.
mem_t MEMARRAY[LEN] = {{0}};


const char ENCODE[256] = { ['c'] = 1, ['g'] = 2, ['t'] = 3,
   ['C'] = 1, ['G'] = 2, ['T'] = 3 };
const char REVCMP[256] = { ['g'] = 1, ['c'] = 2, ['a'] = 3,
   ['G'] = 1, ['C'] = 2, ['A'] = 3 };


// SECTION 1. MACROS //

// Error-handling macros.
#define exit_on_memory_error(x) \
   do { if ((x) == NULL) { fprintf(stderr, "memory error %s:%d:%s()\n", \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); }} while(0)


// SECTION 2. GLOBAL CONSTANTS OF INTEREST //

const char ALPHABET[4] = "ACGT";

const uint8_t NONALPHABET[256] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
};

const char REVCOMP[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0,'T',0,'G',0, 0, 0,'C',0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0,'A',0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0,'T',0,'G',0, 0, 0,'C',0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0,'A',0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
};



// SECTION 3. FUNCTION DEFINITIONS //

// SECTION 3.1 INDEXING FUNCTIONS //

int64_t *
compute_sa   // VISIBLE //
(
   const char * txt
)
{

   const size_t txtlen = strlen(txt) + 1;
   int64_t *array = malloc(txtlen * sizeof(uint64_t));
   exit_on_memory_error(array);

   divsufsort((const unsigned char *) txt, array, txtlen);
   return array;

}


csa_t *
compress_sa
(
   int64_t * sa
)
{

   // The first entry of the suffix array is the length of the text.
   size_t txtlen = sa[0] + 1;

   // Compute the number of required bits.
   size_t nbits = 0;
   while (txtlen > ((uint64_t) 1 << nbits)) nbits++;

   // Compute the number of required bytes and 'uint64_t'.
   size_t nnumb = (txtlen + (16-1)) / 16;
   size_t nint64 = (nbits * nnumb + (64-1)) / 64;
   csa_t * csa = calloc(1, sizeof(csa_t) + nint64 * 8);

   csa->nbits = nbits;
   csa->nint64 = nint64;

   // Set a mask for the 'nbits' lower bits.
   csa->bmask = ((uint64_t) 0xFFFFFFFFFFFFFFFF) >> (64-nbits);
   
   uint8_t lastbit = 0;
   size_t  nb = 0;
   // Sample every 16-th value.
   for (size_t pos = 0 ; pos < txtlen ; pos += 16) {
      int64_t current = sa[pos];
      // Store the compact representation.
      csa->bitf[nb] |= current << lastbit;
      // Update bit offset.
      lastbit += nbits;
      // Check if word is full.
      if (lastbit >= 64) {
         lastbit = lastbit - 64;
         // Complete with remainder or set to 0 (if lastbit = 0).
         // This will clear the upper bits of array.
         csa->bitf[++nb] = current >> (nbits - lastbit);
      }
   }

   return csa;

}


bwt_t *
create_bwt
(
   const char    * txt,
   const int64_t * sa
)
{

   // Allocate new 'bwt_t'.
   const size_t txtlen = strlen(txt) + 1;           // Do not forget $.
   const size_t nslots = (txtlen + (4-1)) / 4;
   const size_t extra = nslots * sizeof(uint8_t);

   bwt_t *bwt = calloc(1, sizeof(bwt_t) + extra);
   exit_on_memory_error(bwt);

   bwt->txtlen = txtlen;
   bwt->nslots = nslots;

   // Encode characters with 2 bits.
   for (size_t pos = 0 ; pos < txtlen ; pos++) {
      if (sa[pos] > 0) { 
         uint8_t c = ENCODE[(uint8_t) txt[sa[pos]-1]];
         bwt->slots[pos/4] |= c << 2*(pos % 4);
      }
      else {
         // Record the position of the zero.
         bwt->zero = pos;
      }
   }

   return bwt;

}


void
write_occ_blocks
(
   occ_t    * occ,
   uint32_t * smpl,
   uint32_t * bits,
   size_t     idx    // Index of 'blocc_t' in array.
)
// Write 'SIGMA' smpl/bits blocks to the 'blocc_t' arrays of 'Occ'
// at position 'pos' (the array index and not the position in
// the bwt).
{
   for (int i = 0 ; i < SIGMA ; i++) {
      occ->rows[i * occ->nrows + idx].smpl = smpl[i];
      occ->rows[i * occ->nrows + idx].bits = bits[i];
   }
}


occ_t *
create_occ
(
   const bwt_t * bwt
)
{

   // Allocate new 'Occ_t'.
   const size_t txtlen = bwt->txtlen;
   const size_t nrows = (txtlen + (32-1)) / 32;
   const size_t extra = SIGMA * nrows * sizeof(blocc_t);

   occ_t * occ = malloc(sizeof(occ_t) + extra);
   exit_on_memory_error(occ);

   occ->txtlen = txtlen;
   occ->nrows = nrows;

   uint32_t smpl[SIGMA] = {0};
   uint32_t diff[SIGMA] = {0};
   uint32_t bits[SIGMA] = {0};

   for (size_t pos = 0 ; pos < bwt->txtlen ; pos++) {
      // Extract symbol at position 'i' from bwt.
      uint8_t c = bwt->slots[pos/4] >> 2*(pos % 4) & 0b11;
      if (pos != bwt->zero) {   // (Skip the '$' symbol).
         diff[c]++;
         bits[c] |= (1 << (31 - pos % 32));
      }
      if (pos % 32 == 31) {     // Write every 32 entries.
         write_occ_blocks(occ, smpl, bits, pos/32);
         memcpy(smpl, diff, SIGMA * sizeof(uint32_t));
         bzero(bits, sizeof(bits));
      }
   }

   write_occ_blocks(occ, smpl, bits, (bwt->txtlen-1)/32);

   // Write 'C'.
   occ->C[0] = 1;
   for (int i = 1 ; i < SIGMA+1 ; i++) {
      occ->C[i] = occ->C[i-1] + diff[i-1];
   }

   return occ;

}


size_t
get_rank
(
   const occ_t   * occ,
         uint8_t   c,
         size_t    pos
)
{
   uint32_t smpl = occ->rows[c*occ->nrows + pos/32].smpl;
   uint32_t bits = occ->rows[c*occ->nrows + pos/32].bits;
   // Several options for pop-count have been tested for this
   // implementation. In the end, I chose '__builtin_popcountl' because
   // the performance is good and the code is simple.
   return occ->C[c] + smpl + __builtin_popcountl(bits >> (31 - pos % 32));
}


void
fill_lut
(
         lut_t   * lut,
   const occ_t   * occ,
   const range_t   range,
   const size_t    depth,
   const size_t    kmerid
)
{
   if (depth >= LUTK) {
      lut->kmer[kmerid] = range;
      return;
   }
   for (uint8_t c = 0 ; c < SIGMA ; c++) {
      size_t bot = get_rank(occ, c, range.bot - 1);
      size_t top = get_rank(occ, c, range.top) - 1;
      fill_lut(lut, occ, (range_t) { .bot=bot, .top=top },
            depth+1, c + (kmerid << 2));
   }
}


// SECTION 3.2 QUERY FUNCTIONS //

range_t
backward_search
(
   const char   * query,
   const size_t   len,
   const occ_t  * occ
)
// Used to search a substring using 'occ' and 'C'.
// In case the query is not found, the condition
// 'range.top - range.bot == -1' is true.
{

   range_t range = { .bot = 1, .top = occ->txtlen-1 };
   int offset = 0;

   for ( ; offset < len ; offset++) {
      int c = ENCODE[(uint8_t) query[len-offset-1]];
      range.bot = get_rank(occ, c, range.bot - 1);
      range.top = get_rank(occ, c, range.top) - 1;
      if (range.top < range.bot)
         return range;
   }

   return range;

}


size_t
query_csa
(
   csa_t  * csa,
   bwt_t  * bwt,
   occ_t  * occ,
   size_t   pos
)
{

   if (pos == bwt->zero) return 0;
   if (pos % 16 == 0) {
      // Value is sampled. Extract it.
      size_t idx = pos / 16;
      size_t lo  = csa->nbits * idx;
      size_t hi  = csa->nbits * (idx+1)-1;
      if (lo/64 == hi/64) {
         // Entry fits in a single 'uint64_t'.
         // Use mask to extract the n-bit encoding of the position.
         return (size_t) csa->bitf[lo/64] >> lo % 64 & csa->bmask;
      }
      else {
         // Entry is split between two 'uint64_t'.
         size_t lo_bits = (size_t) csa->bitf[lo/64] >> lo % 64;
         size_t hi_bits = (size_t) csa->bitf[hi/64] << (64-lo) % 64;
         return (lo_bits | hi_bits) & csa->bmask;
      }
   }
   uint8_t c = bwt->slots[pos/4] >> 2*(pos % 4) & 0b11;
   size_t nextpos = get_rank(occ, c, pos) - 1;
   return query_csa(csa, bwt, occ, nextpos) + 1;

}


char *
normalize_genome
(
   FILE * inputf
)
{

   // Read variables.
   size_t sz = 64;
   ssize_t rlen;
   char * buffer = malloc(64);
   exit_on_memory_error(buffer);

   // Genome storage.
   size_t gsize = 0;
   size_t gbufsize = 64;
   char * genome = malloc(64); 
   exit_on_memory_error(genome);

   // Load fasta file line by line and concatenate.
   while ((rlen = getline(&buffer, &sz, inputf)) != -1) {
      if (buffer[0] == '>') continue;
      if (gbufsize < gsize + rlen) {
         while (gbufsize < gsize + rlen) gbufsize *= 2;
         char * rsz = realloc(genome, gbufsize);
         exit_on_memory_error(rsz);
         genome = rsz;
      }
      int one_if_newline = (buffer[rlen-1] == '\n');
      strncpy(genome + gsize, buffer, rlen - one_if_newline);
      gsize += rlen - one_if_newline;
   }

   // Normalize (use only capital alphabet letters).
   for (size_t pos = 0; pos < gsize ; pos++) {
      int iter = 0;
      if (NONALPHABET[(uint8_t) genome[pos]]) {
         // Replace by cycling over (A,C,G,T).
         genome[pos] = ALPHABET[iter++ % 4];
      }
      else {
         // Use only capital letters (important for
         // sorting the suffixes in lexicographic order).
         genome[pos] = toupper(genome[pos]);
      }
   }

   // Realloc buffer.
   char * rsz = realloc(genome, 2*gsize + 1);
   exit_on_memory_error(rsz);
   genome = rsz;

   // Reverse complement.
   size_t div = gsize;
   for (size_t pos = 0 ; pos < div ; pos++)
      genome[div + pos] = REVCOMP[(uint8_t) genome[div-pos-1]];

   gsize = 2*gsize + 1;

   // Add the terminator.
   genome[2*div] = '\0';

   // Clean up.
   free(buffer);

   return genome;

}


int
banded_needleman_wunsch
(
   const char   * seq1,
   const char   * seq2,
   const int      len1,
   const size_t   sz,
   const int      cutoff
)
{

   if (sz > 128) exit(EXIT_FAILURE);

   int x;  // Mis/Match.
   int y;  // Insertion.
   int z;  // Deletion.

   int a1[257] = {0};
   int a2[257] = {0};

   int * LCURR = a1 + 128;
   int * LPREV = a2 + 128;

   for (int i = 0 ; i < len1 ; i++) {
      if (i > sz) {
         x = LPREV[-sz] + (seq1[i] != seq2[i-sz]);
         y = LPREV[-sz+1] + 1;
         LCURR[-sz] = x < y ? x : y;
         x = LPREV[sz] + (seq1[i-sz] != seq2[i]);
         z = LPREV[sz-1] + 1;
         LCURR[-sz] = x < z ? x : z;
      }
      for (int j = -sz+1; j < 0 ; j++) {
         if (i+j < 0) continue;
         x = LPREV[j] + (seq1[i] != seq2[i+j]);
         y = LPREV[j+1] + 1;
         z = LCURR[j-1] + 1;
         LCURR[j] = x < y ? (x < z ? x : z) : (y < z ? y : z);
      }
      for (int j = sz-1 ; j > 0 ; j--) {
         if (i-j < 0) continue;
         x = LPREV[j] + (seq1[i-j] != seq2[i]);
         y = LPREV[j-1] + 1;
         z = LCURR[j+1] + 1;
         LCURR[j] = x < y ? (x < z ? x : z) : (y < z ? y : z);
      }
      x = LPREV[0] + (seq1[i] != seq2[i]);
      y = LCURR[-1] + 1;
      z = LCURR[1] + 1;
      LCURR[0] = x < y ? (x < z ? x : z) : (y < z ? y : z);

      if (LCURR[0] > cutoff) return cutoff + 1 ;

      int * tmp = LPREV;
      LPREV = LCURR;
      LCURR = tmp;

   }

   return LCURR[0];

}


void
align
(
   const index_t   idx,
   const char    * seq,
   const char    * genome,
   const int       nmem
)
{
   for (int i = 0 ; i < nmem ; i++) {
      mem_t mem = MEMARRAY[i];
      // We still kinda need to chain the seeds.
      for (size_t pos = mem.range.bot ; pos <= mem.range.top ; pos++) {
         size_t hitpos = query_csa(idx.csa, idx.bwt, idx.occ, pos);
         fprintf(stderr, "%.*s\n", (int) (mem.end - mem.beg + 1),
               seq + mem.beg);
         fprintf(stderr, "%.*s\n--\n", (int) (mem.end - mem.beg + 1),
               genome + hitpos);
      }
      // When we have the seed(s) that gave the best hit, we need to
      // evaluate the number of copies of that sequence, and the
      // divergence. We said we would do this with direct computation of
      // the log-likelihood like: set mu, compute N and get loglik.
      // For this, we need to keep more dense record of the seeding
      // process for every seed.
   }
}


void
mapread
(
   const char    * seq,
   const index_t   idx,
   const char    * genome,
   const size_t    gamma
)
{

   fprintf(stderr, "\nprocesssing: %s\n", seq);

   int len = strlen(seq);
   assert(len <= LEN);

   int rend = len-1;

   range_t range = {0};
   range_t newrange = {0};

   // Number of MEM seeds.
   int nmem = 0;

   while (1) {

      // Grab a new struct of type 'mem_t'.
      mem_t mem = {0};

      mem.end = rend;

      // Backward <<<
      range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
      int roffset = 0;

      // Look up the beginning (reverse) of the query in lookup table.
      if (rend >= LUTK) {
         size_t merid = 0;
         for ( ; roffset < LUTK ; roffset++) {
            merid = (merid << 2) +
               ENCODE[(uint8_t) seq[rend - roffset]];
         }
         range = idx.lut->kmer[merid];
         mem.fwd[LUTK-1] = range.top - range.bot + 1;
      }

      // Cancel if we went too far already.
      if (range.top - range.bot < 1) {
         range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
         roffset = 0;
      }

      for ( ; roffset <= rend ; roffset++) {
         int c = ENCODE[(uint8_t) seq[rend-roffset]];
         newrange.bot = get_rank(idx.occ, c, range.bot - 1);
         newrange.top = get_rank(idx.occ, c, range.top) - 1;
         mem.fwd[roffset-1] = newrange.top - newrange.bot + 1;
         // Stop if less than 1 hit.
         if (newrange.top < newrange.bot)
            break;
         range = newrange;
      }

      mem.beg = rend - roffset + 1;
      mem.range = range;

      // Forward >>>
      range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
      int foffset = 0;

      if (roffset >= LUTK) {
         size_t merid = 0;
         for ( ; foffset < LUTK ; foffset++) {
            merid = (merid << 2) +
               REVCMP[(uint8_t) seq[rend - roffset + foffset + 1]];
         }
         range = idx.lut->kmer[merid];
         mem.bwd[LUTK-1] = range.top - range.bot + 1;
      }

      // Cancel if we went too far already.
      if (range.top - range.bot < 1) {
         range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
         foffset = 0;
      }

      for ( ; foffset < roffset-1 ; foffset++) {
         int c = REVCMP[(uint8_t) seq[rend - roffset + foffset + 1]];
         range.bot = get_rank(idx.occ, c, range.bot - 1);
         range.top = get_rank(idx.occ, c, range.top) - 1;
         mem.bwd[LUTK-1] = range.top - range.bot + 1;
         // Stop if less than 2 hits.
         if (range.top - range.bot < 1)
            break;
      }

      // Overwrite MEM if it is too short.
      if (roffset >= gamma)
         MEMARRAY[nmem++] = mem;

      if (roffset >= rend)
         break;

      rend += - roffset + foffset - 1;

   }

   align(idx, seq, genome, nmem);

}
