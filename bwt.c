#define _GNU_SOURCE
#include <assert.h>
#include "divsufsort.h"
#include "bwt.h"


// FIXME: The logic of the 'mem_t' struct is very weak and does not
// FIXME: stand a chance on something else than the test data set.


const char ENCODE[256] = { ['c'] = 1, ['g'] = 2, ['t'] = 3,
			   ['C'] = 1, ['G'] = 2, ['T'] = 3 };
const char REVCMP[256] = { ['g'] = 1, ['c'] = 2, ['a'] = 3,
			   ['G'] = 1, ['C'] = 2, ['A'] = 3 };


// SECTION 1. MACROS //

// Error-handling macros.
#define exit_on_memory_error(x)						\
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

chr_t *
index_load_chr
(
 const char * index_file
)
{

   // Read chromosome index.
   char chr_file[256];
   sprintf(chr_file, "%s.chr", index_file);
   // Files
   int fd = open(chr_file,O_RDONLY);

   // Alloc structure.
   chr_t * chr = malloc(sizeof(chr_t));
   ssize_t b;
   b = read(fd, &(chr->gsize), sizeof(size_t));
   if (b < 1) {
      fprintf(stderr, "error reading chr index file\n");
      exit(EXIT_FAILURE);
   }
   b = read(fd, &(chr->nchr), sizeof(size_t));
   
   size_t * start = malloc(chr->nchr * sizeof(size_t));
   char  ** names = malloc(chr->nchr * sizeof(char *));
   exit_on_memory_error(start);
   exit_on_memory_error(names);

   for (size_t i = 0; i < chr->nchr; i++) {
      b = read(fd, start+i, sizeof(size_t));
      size_t slen;
      b = read(fd, &slen, sizeof(size_t));
      char * name = malloc(slen);
      exit_on_memory_error(name);
      b = read(fd, name, slen);
      names[i] = name;
   }

   chr->start = start;
   chr->name = names;

   close(fd);

   return chr;
}

size_t
bisect_search
(
 size_t   start,
 size_t   end,
 size_t * set,
 size_t   value
)
{
   if (start >= end - 1) {
      if (value < set[start]) return start;
      return (value >= set[end] ? end : start) + 1;
   }
   size_t middle = (start + end) / 2;
   if (set[middle] >= value) return bisect_search(start, middle, set, value);
   else return bisect_search(middle, end, set, value);
}

char *
chr_string
(
 size_t   refpos,
 chr_t  * chr
)
{
   if (chr->nchr == 0)
      return NULL;

   char strand = '+';
   if (refpos > chr->gsize) {
      strand = '-';
      refpos = 2*chr->gsize - refpos - 1;
   }
      
   int chrnum = 0;
   if (chr->nchr > 1) {
      chrnum = bisect_search(0, chr->nchr-1, chr->start, refpos+1)-1;
   }
   size_t chrpos = refpos+1 - (chr->start[chrnum]-1);
   int   slen    = snprintf(NULL, 0, "%s:%ld:%c", chr->name[chrnum], chrpos, strand);
   char * str    = malloc(slen+1);
   sprintf(str, "%s:%ld:%c", chr->name[chrnum], chrpos, strand);
   return str;
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
 FILE   * inputf,
 char   * chrfile
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

   // Chromosome file.
   stack_t * seqnames = NULL;
   stack_t * seqstart = NULL;
   
   if (chrfile) {
      seqnames = stack_new(64);
      seqstart = stack_new(64);
   }
   
   // Load fasta file line by line and concatenate.
   while ((rlen = getline(&buffer, &sz, inputf)) != -1) {
      if (buffer[0] == '>') {
	 if (chrfile) {
	    buffer[rlen-1] = 0;
	    int k = 0;
	    while (buffer[k] != ' ' && buffer[k] != 0) k++;
	    buffer[k] = 0;
	    char * seqname = strdup(buffer+1);
	    push(seqname, &seqnames);
	    size_t * gpos = malloc(sizeof(size_t));
	    *gpos = gsize+1;
	    push(gpos, &seqstart);
	 }
	 continue;
      }
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

   if (chrfile) {
      int fd = open(chrfile, O_WRONLY | O_CREAT, 0664);
      ssize_t b;
      b = write(fd, &gsize, sizeof(size_t));
      if (b < 1) {
	 fprintf(stderr, "error writing chr index file\n");
	 exit(EXIT_FAILURE);
      }
	 
      b = write(fd, &(seqnames->pos), sizeof(size_t));
      for (size_t i = 0; i < seqnames->pos; i++) {
	 b = write(fd, seqstart->ptr[i], sizeof(size_t));
	 char * seqname = (char *) seqnames->ptr[i];
	 size_t slen    = strlen(seqname)+1;
	 b = write(fd, &slen, sizeof(size_t));
	 b = write(fd, seqname, slen);
	 free(seqstart->ptr[i]);
	 free(seqnames->ptr[i]);
      }
      free(seqstart);
      free(seqnames);
      close(fd);
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

stack_t *
stack_new
(
 size_t max
 )
{
   size_t base = sizeof(stack_t);
   size_t extra = max * sizeof(void *);
   stack_t * stack = malloc(base + extra);
   exit_on_memory_error(stack);

   stack->max = max;
   stack->pos = 0;
   return stack;
}

void
push
(
 void     * ptr,
 stack_t ** stackp
 )
{
   stack_t * stack = *stackp;
   if (stack->pos >= stack->max) {
      size_t newmax = stack->max*2;
      stack = *stackp = realloc(stack,
				sizeof(stack_t)+newmax*sizeof(void *));
      exit_on_memory_error(stack);
      stack->max = newmax;
   }

   stack->ptr[stack->pos++] = ptr;
   return;
}
