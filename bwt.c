#define _GNU_SOURCE
#include <assert.h>
#include "divsufsort.h"
#include "bwt.h"


// FIXME: The logic of the 'mem_t' struct is very weak and does not
// FIXME: stand a chance on something else than the test data set.
typedef struct memstack_t memstack_t;
typedef struct seed_t seed_t;

#define LEN 50

#define min(x,y) ((x) < (y) ? (x) : (y))
#define min3(x,y,z) (min(min(x,y),z))

struct memstack_t {
   size_t pos;
   size_t max;
   mem_t  mem[];
};

struct seed_t {
   size_t refpos;
   mem_t  mem;
};

memstack_t *  memstack_new (size_t max);
void          mem_push     (mem_t mem, memstack_t ** stackp);
alnstack_t *  alnstack_new (size_t max);
void          aln_push     (aln_t aln, alnstack_t ** stackp);

int           seed_by_refpos (const void * a, const void * b) {
   return ((seed_t *)a)->refpos > ((seed_t *)b)->refpos;
};

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

const char CAPS[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0,'A',0,'C',0, 0, 0,'G',0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0,'T',0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0,'A',0,'C',0, 0, 0,'G',0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0,'T',0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
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


void
recursive_csa_query
(
   csa_t  * csa,
   bwt_t  * bwt,
   occ_t  * occ,
   range_t  range,
   size_t * sa_values,
   size_t   path_offset
)
{
   // Get sampled SA values.
   // Compute offset to the closest %16.
   size_t offset = (16 - range.bot%16)%16;
   size_t pos    = range.bot + offset;
   while (pos <= range.top) {
      // Extract sampled value.
      size_t idx = pos / 16;
      size_t lo  = csa->nbits * idx;
      size_t hi  = csa->nbits * (idx+1)-1;
      if (lo/64 == hi/64) {
         // Entry fits in a single 'uint64_t'.
         // Use mask to extract the n-bit encoding of the position.
         sa_values[offset] = (size_t) csa->bitf[lo/64] >> lo % 64 & csa->bmask;
      }
      else {
         // Entry is split between two 'uint64_t'.
         size_t lo_bits = (size_t) csa->bitf[lo/64] >> lo % 64;
         size_t hi_bits = (size_t) csa->bitf[hi/64] << (64-lo) % 64;
         sa_values[offset] = (lo_bits | hi_bits) & csa->bmask;
      }
      // Add path offset to get original position.
      sa_values[offset] += path_offset;
      // Find new sampled position.
      offset += 16;
      pos    += 16;
   }

   // If all positions are covered, return.
   int    done = 1;
   size_t rlen = range.top - range.bot + 1;
   for (int i = 0; i < rlen; i++) {
      if (!sa_values[i]) {
	 done = 0;
	 break;
      }
   }

   if (done) return;

   // Otherwise extend unsampled paths.
   uint8_t * prev_c  = malloc(rlen * sizeof(uint8_t));
   uint8_t * extend  = calloc(4, sizeof(uint8_t));
   int     * c_count = calloc(4, sizeof(int));
   exit_on_memory_error(prev_c);
   exit_on_memory_error(extend);
   exit_on_memory_error(c_count);

   // For each position in range:
   for (size_t pos = range.bot, i = 0; pos <= range.top; pos++, i++) {
      // 1. Get preceding character from BWT.
      prev_c[i] = bwt->slots[pos/4] >> 2*(pos % 4) & 0b11;
      c_count[prev_c[i]]++;
      // 2. Mark the nucleotide for extension if SA value is missing.
      if (!sa_values[i]) extend[(int)prev_c[i]] = 1;
   }

   // Extend unsampled nucleotides.
   for (int c = 0; c < 4; c++) {
      if(!extend[c]) continue;

      size_t * idx_c       = malloc(c_count[c]*sizeof(size_t));
      size_t * sa_values_c = malloc(c_count[c]*sizeof(size_t));
      exit_on_memory_error(idx_c);
      exit_on_memory_error(sa_values_c);

      for (int i = 0, j = 0; i < rlen; i++) {
	 if (prev_c[i] != c) continue;
	 // Assemble new sa_values vector for the current nucleotide.
	 sa_values_c[j] = sa_values[i];
	 // Store the positions of the current nucleotide in a list.
	 idx_c[j++] = i;
      }

      // Get new range.
      range_t newrange;
      newrange.bot = get_rank(occ, c, range.bot - 1);
      newrange.top = get_rank(occ, c, range.top) - 1;

      // Recursive call.
      recursive_csa_query(csa, bwt, occ, newrange, sa_values_c, path_offset+1);

      // Place sa_values back in their original position.
      for (int i = 0; i < c_count[c]; i++) {
	 sa_values[idx_c[i]] = sa_values_c[i];
      }

      // Free memory.
      free(idx_c);
      free(sa_values_c);
   }

   free(prev_c);
   free(extend);
   free(c_count);
}


size_t *
query_csa_range
(
   csa_t  * csa,
   bwt_t  * bwt,
   occ_t  * occ,
   range_t  range
)
{

   size_t * sa_values = calloc((range.top - range.bot + 1), sizeof(size_t));
   exit_on_memory_error(sa_values);

   recursive_csa_query(csa, bwt, occ, range, sa_values, 0);

   return sa_values;
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
nw
(
   const char   * seq1,
   const char   * seq2,
   const int      len1,
   const int      len2,
   const int      cutoff
)
// Place reference in seq1, you can add extra length to len1 to allocate insertions in seq2.
{
   // Penalties (don't set i_p or d_p smaller than m_p!).
   const int m_p = 1, i_p = 1, d_p = 1;

   int   len = min(len1, len2);
   int * rowp = calloc(len1+2, sizeof(int));
   int * colp = calloc(len2+2, sizeof(int));
   // Comment for SW.
   for (int i = 1; i < len1+2; i++) rowp[i] = rowp[i-1] + d_p;
   for (int j = 1; j < len2+2; j++) colp[j] = colp[j-1] + i_p;

   // Align row and col.
   int * row = rowp + 1;
   int * col = colp + 1;

   int score = cutoff;
   for (int pos = 0; pos < len; pos++) {
      score = cutoff;
      // Update row.
      int i, left, diag;
      for (i = pos, left = col[pos], diag = row[i-1]; i < len1 && row[i-1] <= cutoff; i++) {
	 int next_diag = row[i];
	 left = row[i] = min3(diag + (CAPS[(int)seq1[i]] != CAPS[(int)seq2[pos]])*m_p,
			 row[i] + i_p,
			 left + d_p);
	 diag = next_diag;
	 if (row[i] < score) score = row[i];
      }
      // Uncomment for SW.
      //      row[i] = row[i-1];

      // Update column.
      int j, up;
      for (j = pos+1, up = row[pos], diag = col[pos]; j < len2 && col[j-1] <= cutoff; j++) {
	 int next_diag = col[j];
	 up = col[j] = min3(diag + (CAPS[(int)seq1[pos]] != CAPS[(int)seq2[j]])*m_p,
			    up + i_p,
			    col[j] + d_p);
	 diag = next_diag;
	 if (col[j] < score) score = col[j];
      }
      // Uncomment for SW.
      //      col[j] = col[j-1];
   }

   free(rowp);
   free(colp);

   return score;
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


alnstack_t *
align
(
   const index_t      idx,
   const char       * seq,
   const char       * genome,
   const memstack_t * mems
)
{

   int slen = strlen(seq);
   // Count loci.
   
   size_t nloc = 0;
   for (int i = 0 ; i < mems->pos ; i++) {
      mem_t mem = mems->mem[i];
      nloc += mem.range.top - mem.range.bot + 1;
   }

   // Allocate align positions (one per locus).
   seed_t * seeds = malloc(nloc * sizeof(seed_t));
   exit_on_memory_error(seeds);

   // Query SA.
   for (int i = 0, j = 0; i < mems->pos; i++) {
      mem_t mem = mems->mem[i];
      // We still kinda need to chain the seeds.
      /*
      // THIS USES SINGLE CSA QUERY.
      for (size_t pos = mem.range.bot ; pos <= mem.range.top ; pos++) {
         size_t hitpos = query_csa(idx.csa, idx.bwt, idx.occ, pos);
         seeds[j++] = (seed_t) {hitpos, mem};
      }
      */

      // RANGE CSA QUERY.
      size_t * sa_values = query_csa_range(idx.csa, idx.bwt, idx.occ, mem.range);
      for (size_t k = 0 ; k < mem.range.top - mem.range.bot + 1 ; k++) {
         seeds[j++] = (seed_t) {sa_values[k], mem};
      }
   }

   // Sort by align position in the genome.
   qsort(seeds, nloc, sizeof(seed_t), seed_by_refpos);

   // Best alignment (low score is good).
   int best_score = slen+1;
   size_t last_locus = 0;
   alnstack_t * best = alnstack_new(10);
   
   // Align and chain mems.
   for (int i = 0; i < nloc; i++) {
      seed_t seed = seeds[i];
      mem_t mem = seed.mem;
      if (last_locus <= seed.refpos) { // Align.
         const char * ref = genome + seed.refpos - mem.beg;
         int score = nw(ref,
               seq,
               slen+3, // Allow 3 nucleotides to allocate insertions.
               slen,
               best_score + 1
         );
	 
         // VERBOSE ALIGNMENT (DEBUG)
         fprintf(stderr, "%.*s %.*s %.*s\n",
               (int)mem.beg, seq, (int) (mem.end - mem.beg + 1),
               seq + mem.beg, (int)(slen-mem.end-1), seq + mem.end + 1);
         fprintf(stderr, "%.*s %.*s %.*s\nscore: %d\n--\n",
               (int)mem.beg, ref, (int) (mem.end - mem.beg + 1),
               ref + mem.beg, (int)(slen-mem.end-1), ref + mem.end + 1,
               score);

         // Check align score.
         if (score <= best_score) {
            // Reset chain anchor.
            last_locus = seed.refpos + (slen - mem.beg);
            // Create new alignment.
            aln_t aln;
            aln.score  = score;
            aln.nmem   = 1;
            aln.refpos = seed.refpos;
            aln.refseq = ref;
            aln.mem[0] = mem;

            if (score < best_score) {
               best_score = score;
               // Reset best align stack.
               best->pos = 0;
               aln_push(aln, &best);
            } else {
               // Add alignment to best.
               aln_push(aln, &best);
            }
         }

      } else { // Chain. (Add mem to alignment)
         fprintf(stderr, "* chain seed - alignment skipped\n--\n");
         best->aln[best->pos-1].mem[best->aln[best->pos-1].nmem++] = mem;
      }
   }
      
   // When we have the seed(s) that gave the best hit, we need to
   // evaluate the number of copies of that sequence, and the
   // divergence. We said we would do this with direct computation of
   // the log-likelihood like: set mu, compute N and get loglik.
   // For this, we need to keep more dense record of the seeding
   // process for every seed.

   free(seeds);
   return best;

}


alnstack_t *
mapread
(
   const char    * seq,
   const index_t   idx,
   const char    * genome,
   const size_t    gamma
)
{

   int len = strlen(seq);

   // TODO: the assert is super weak.
   assert(len <= LEN);

   int end = len-1;

   range_t range = {0};
   range_t newrange = {0};

   // Number of MEM seeds.
   memstack_t * mems = memstack_new(50);

   while (1) {

      // Grab a new struct of type 'mem_t'.
      mem_t mem = {0};

      mem.end = end;

      // Backward <<<
      range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
      int mpos = end, mlen = 0;

      // Look up the beginning (reverse) of the query in lookup table.
      if (end >= LUTK - 1) {
         size_t merid = 0;
         for ( ; mlen < LUTK ; mlen++, mpos--) {
            uint8_t c = ENCODE[(uint8_t) seq[end-mlen]];
            merid = c + (merid << 2);
         }
         range = idx.lut->kmer[merid];
         mem.left[LUTK-1] = range.top - range.bot + 1;
      }

      // Cancel if we went too far already.
      if (range.top < range.bot) {
         range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
         mpos = end; mlen = 0;
      }

      for ( ; mpos >= 0 ; mpos--, mlen++) {
         int c = ENCODE[(uint8_t) seq[mpos]];
         newrange.bot = get_rank(idx.occ, c, range.bot - 1);
         newrange.top = get_rank(idx.occ, c, range.top) - 1;
         // Collect (reverse) cascade data.
         mem.left[mlen] = newrange.top - newrange.bot + 1;
         // Stop if no hits.
         if (newrange.top < newrange.bot)
            break;
         range = newrange;
      }

      mem.beg = ++mpos;
      mem.range = range;

      // Forward >>>
      range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
      mlen = 0;

      if (end >= LUTK - 1) {
         size_t merid = 0;
         for ( ; mlen < LUTK ; mpos++, mlen++) {
            uint8_t c = REVCMP[(uint8_t) seq[mpos]];
            merid = c + (merid << 2);
         }
         range = idx.lut->kmer[merid];
         mem.right[LUTK-1] = range.top - range.bot + 1;
      }

      // Cancel if we went too far already.
      if (range.top < range.bot) {
         range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
         mpos = mem.beg; mlen = 0;
      }

      for ( ; mpos <= end ; mpos++, mlen++) {
         int c = REVCMP[(uint8_t) seq[mpos]];
         range.bot = get_rank(idx.occ, c, range.bot - 1);
         range.top = get_rank(idx.occ, c, range.top) - 1;
         // Collect (forward) cascade data.
         mem.right[mlen] = range.top - range.bot + 1;
         // Stop if less than bw_hits.
         if (range.top < range.bot)
            break;
      }

      // Keep MEM if above minimum size.
      if (mlen >= gamma) mem_push(mem, &mems);

      if (mem.beg < 1) break;

      // Find new end position (forward).
      range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
      end = mem.beg - 1;
      while (1) {
         int c = REVCMP[(uint8_t) seq[end]];
         range.bot = get_rank(idx.occ, c, range.bot - 1);
         range.top = get_rank(idx.occ, c, range.top) - 1;
         if (range.top < range.bot) {
            end--;
            break;
         }
         end++;
      }

      if (end + 1 < gamma) break;

   }

   fprintf(stderr, "MEMS found: %ld\n-- Alignments --\n", mems->pos);

   // Return the best alignment(s) in an alignment stack.
   alnstack_t * aln = align(idx, seq, genome, mems);

   free(mems);
   return aln;

}


memstack_t *
memstack_new
(
 size_t max
)
{
   size_t base = sizeof(memstack_t);
   size_t extra = max * sizeof(mem_t);
   memstack_t * stack = malloc(base + extra);
   exit_on_memory_error(stack);

   stack->max = max;
   stack->pos = 0;
   return stack;

}


void
mem_push
(
  mem_t mem,
  memstack_t ** stackp
)
{

   memstack_t * stack = *stackp;

   if (stack->pos >= stack->max) {
      size_t newmax = stack->max*2;
      stack = *stackp = realloc(stack,
            sizeof(memstack_t)+newmax*sizeof(mem_t));
      exit_on_memory_error(stack);
      stack->max = newmax;
   }

   stack->mem[stack->pos++] = mem;
   return;

}


alnstack_t *
alnstack_new
(
 size_t max
)
{
   size_t base = sizeof(alnstack_t);
   size_t extra = max * sizeof(aln_t);
   alnstack_t * stack = malloc(base + extra);
   exit_on_memory_error(stack);

   stack->max = max;
   stack->pos = 0;
   return stack;

}


void
aln_push
(
  aln_t aln,
  alnstack_t ** stackp
)
{
   alnstack_t * stack = *stackp;
   if (stack->pos >= stack->max) {
      size_t newmax = stack->max*2;
      stack = *stackp = realloc(stack,
            sizeof(alnstack_t)+newmax*sizeof(aln_t));
      exit_on_memory_error(stack);
      stack->max = newmax;
   }

   stack->aln[stack->pos++] = aln;
   return;
}
