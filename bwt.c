#define _GNU_SOURCE
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include "bwt.h"


// FIXME: The logic of the 'mem_t' struct is very weak and does not
// FIXME: stand a chance on something else than the test data set.


const char ENCODE[256] = { ['c'] = 1, ['g'] = 2, ['t'] = 3,
			   ['C'] = 1, ['G'] = 2, ['T'] = 3 };
const char REVCMP[256] = { ['g'] = 1, ['c'] = 2, ['a'] = 3,
			   ['G'] = 1, ['C'] = 2, ['A'] = 3 };


// SECTION 1. MACROS //

#define min(x,y) (x) < (y) ? (x) : (y)

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

void bwt2occ
(
 const char * basename
)
{
   char * fn = malloc(strlen(basename)+10);

   // Load bwt
   sprintf(fn, "%s.bwt", basename);
   int fbwt = open(fn, O_RDONLY);
   if (fbwt < 0) exit_cannot_open(fn);

   size_t mmsz = lseek(fbwt, 0, SEEK_END);
   bwt_t *bwt = (bwt_t *) mmap(NULL, mmsz, PROT_READ, MMAP_FLAGS, fbwt, 0);
   exit_error(bwt == NULL);
   close(fbwt);

   // Open OCC stream
   size_t   txtlen = bwt->txtlen + 1; // Include wildcard, not in bwt
   uint64_t word_size = 64;
   uint64_t mark_bits = OCC_INTERVAL_SIZE * word_size;
   const uint64_t nintv = (txtlen + mark_bits - 1) / mark_bits;
   const uint64_t nword = nintv * OCC_INTERVAL_SIZE;
   const uint64_t nmark = nintv + 1;
   const size_t occ_size = nword+nmark;

   sprintf(fn, "%s.occ", basename);
   int occ_fd = creat(fn, 0644);

   uint64_t C[SIGMA+1] = {0};
   
   // Struct header
   ssize_t bw, bt;
   uint64_t occ_mark_intv = OCC_INTERVAL_SIZE;
   bw = write(occ_fd, &txtlen, sizeof(size_t));
   bw = write(occ_fd, &occ_size, sizeof(uint64_t));
   bw = write(occ_fd, &occ_mark_intv, sizeof(uint64_t));
   bw = write(occ_fd, &word_size, sizeof(uint64_t));
   bw = write(occ_fd, &mark_bits, sizeof(uint64_t));
   bw = write(occ_fd, C, (SIGMA+1)*sizeof(uint64_t));

   // Loop variables
   int intv_words = SIGMA*(1+OCC_INTERVAL_SIZE);
   uint64_t * occ_intv = calloc(intv_words, sizeof(uint64_t));
   exit_on_memory_error(occ_intv);
   uint64_t occ_abs[SIGMA] = {0};

   // Read bwt
   size_t i = 0, j = 0;
   for (; i < txtlen; i++, j++) {
      // Fill occ word
      if (i != bwt->zero) {
	 uint8_t c = (bwt->slots[j/4] >> (6-2*(j%4))) & 0x03;
	 occ_abs[c]++;
	 occ_intv[SIGMA + c*OCC_INTERVAL_SIZE + (i%mark_bits)/word_size] |= ((uint64_t)1) << (word_size - 1 - (i % word_size));
      } else {
	 // Wildcard is not present in bwt, leave a blank space in occ table
	 j--;
      }

      // Write OCC word
      if ((i+1) % mark_bits == 0) {
	 // Write full interval
	 bw = 0; bt = intv_words*sizeof(uint64_t);
	 while(bw<bt) bw += write(occ_fd, ((char *)occ_intv)+bw, bt-bw);
	 // Prepare next interval
	 memset(occ_intv, 0, bt);
	 // Copy absolute position markers
	 memcpy(occ_intv, occ_abs, SIGMA*sizeof(uint64_t));
      }

      // Verbose progress
      if (i%10000 == 0)
	 fprintf(stderr, "\r%.2f%%", i*100.0/txtlen);
   }

   // Finish OCC
   if (i % mark_bits > 0) {
      // Write full interval
      bw = 0, bt = intv_words*sizeof(uint64_t);
      while(bw<bt) bw += write(occ_fd, ((char *)occ_intv)+bw, bt-bw);
   }
   // Write last OCC absolute position mark
   bw = 0; bt = SIGMA*sizeof(uint64_t);
   while(bw<bt) bw += write(occ_fd, ((char *)occ_abs)+bw, bt-bw);

   // Compute 'C'.
   C[0] = 1;
   for (int i = 1 ; i < SIGMA+1 ; i++)
      C[i] = C[i-1] + occ_abs[i-1];

   // Seek back in occ_fd to write C table
   lseek(occ_fd, sizeof(size_t) + 4*sizeof(uint64_t), SEEK_SET);
   bw = 0; bt = (SIGMA+1)*sizeof(uint64_t);
   while (bw<bt) bw += write(occ_fd, ((char *)C)+bw, bt-bw);

   munmap(bwt, mmsz);
   close(occ_fd);
   free(occ_intv);
   free(fn);
}

void
bwt2sa
(
 const char * basename
)
{
   char * fn = malloc(strlen(basename)+10);
   // Load sequence data
   sprintf(fn, "%s.dna", basename);
   int fdna = open(fn, O_RDONLY);
   if (fdna < 0) exit_cannot_open(fn);

   size_t mmsz_dna = lseek(fdna, 0, SEEK_END);
   char *dna = (char *) mmap(NULL, mmsz_dna, PROT_READ, MMAP_FLAGS, fdna, 0);
   exit_error(dna == NULL);
   close(fdna);

   // Load occ table
   sprintf(fn, "%s.occ", basename);
   int focc = open(fn, O_RDONLY);
   if (focc < 0) exit_cannot_open(fn);

   size_t mmsz_occ = lseek(focc, 0, SEEK_END);
   occ_t *occ = (occ_t *) mmap(NULL, mmsz_occ, PROT_READ, MMAP_FLAGS, focc, 0);
   exit_error(occ == NULL);
   close(focc);

   // Open CSA stream
   size_t nbits = 0;
   while (occ->txtlen > ((uint64_t) 1 << nbits)) nbits++;
   size_t nnumb = (occ->txtlen + (CSA_SAMP_RATIO-1)) / CSA_SAMP_RATIO;
   size_t nint64 = (nbits * nnumb + (64-1)) / 64;
   uint64_t bmask = ((uint64_t) 0xFFFFFFFFFFFFFFFF) >> (64-nbits);
   
   sprintf(fn,"%s.sa", basename);
   int csa_fd = creat(fn, 0644);

   // CSA struct header
   ssize_t bw, bt;
   size_t csa_ratio = CSA_SAMP_RATIO;
   bw = write(csa_fd, &nbits, sizeof(size_t));
   bw = write(csa_fd, &bmask, sizeof(uint64_t));
   bw = write(csa_fd, &nint64, sizeof(size_t));
   bw = write(csa_fd, &csa_ratio, sizeof(size_t));

   // CSA buffer
   uint64_t * csa_buf = calloc(nint64, sizeof(uint64_t));

   // Perform backward search from the end
   int64_t fm_ptr = 0;

   for (size_t i = occ->txtlen-2; i > 0; i--) {
      int c = dna[i/4] >> (6-2*(i%4)) & 0x03;
      fm_ptr = get_rank(occ, c, fm_ptr - 1);

      // CSA (1:csa_samp_rate compression)
      if (fm_ptr%CSA_SAMP_RATIO == 0) {
	 uint64_t saidx = i & bmask;
	 uint64_t word = ((fm_ptr/CSA_SAMP_RATIO)*nbits)/64;
	 uint64_t bit = ((fm_ptr/CSA_SAMP_RATIO)*nbits)%64;

	 // Write SA value (i) in sampled FM index position
	 csa_buf[word]   |= saidx << bit;
	 if (bit + nbits > 64)  // Split word
	    csa_buf[word+1] |= saidx >> (64-bit);
      }

      // Verbose progress
      if (i%10000 == 0)
	 fprintf(stderr, "\r%.2f%%", 100-i*100.0/occ->txtlen);
   }

   // Write sampled SA to file
   bw = 0; bt = nint64*sizeof(int64_t);
   while(bw < bt) bw += write(csa_fd, ((char *)csa_buf) + bw, bt-bw);

   close(csa_fd);
   munmap(dna, mmsz_dna);
   munmap(occ, mmsz_occ);
   free(csa_buf);
   free(fn);
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
   int fd = open(chr_file, O_RDONLY);

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

void
free_index_chr
(
 chr_t * chr
 )
{
   for (size_t i = 0; i < chr->nchr; i++) {
      free(chr->name[i]);
   }
   free(chr->name);
   free(chr->start);
   free(chr);
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
 chr_t  * chr,
 char   * buffer
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
   snprintf(buffer, 255, "%s:%ld:%c", chr->name[chrnum], chrpos, strand);
   return buffer;
}

pos_t
get_pos
(
 size_t   refpos,
 chr_t  * chr
 )
{
   if (chr->nchr == 0)
      return (pos_t) {0};

   int strand = 1; // Forward.

   if (refpos > chr->gsize) {
      strand = 0; // Reverse.
      refpos = 2*chr->gsize - refpos - 1;
   }
      
   int chrnum = 0;
   if (chr->nchr > 1) {
      chrnum = bisect_search(0, chr->nchr-1, chr->start, refpos+1)-1;
   }
   size_t chrpos = refpos+1 - (chr->start[chrnum]-1);
   return (pos_t) {
      .rname = chr->name[chrnum],
	 .pos = chrpos,
	 .strand = strand
	 };
}

size_t
get_rank
(
 const occ_t  * occ_ptr,
 uint8_t  c,
 size_t   pos_unsigned
 )
{
   int64_t pos = (int64_t)pos_unsigned;
   if (pos == -1) return (size_t)occ_ptr->C[c];
   // Compute interval, word, bit and absolute marker positions.
   int64_t intv_num = pos/occ_ptr->occ_mark_bits;
   int64_t intv_wrd = (pos%occ_ptr->occ_mark_bits)/occ_ptr->occ_word_size;
   int64_t wrdptr = intv_num*SIGMA*(occ_ptr->occ_mark_intv+1) + SIGMA + c*occ_ptr->occ_mark_intv + intv_wrd;
   int64_t mrkptr = (intv_num + (intv_wrd >= occ_ptr->occ_mark_intv/2)) * SIGMA*(occ_ptr->occ_mark_intv + 1) + c;
   int64_t bit = pos%occ_ptr->occ_word_size;

   uint64_t occ = occ_ptr->occ[mrkptr];
   if (wrdptr > mrkptr) {
      int64_t  offset = 0;
      // Sum bit offsets.
      for (uint64_t i = wrdptr-intv_wrd; i < wrdptr; i++)
	 offset += __builtin_popcountl(occ_ptr->occ[i]);
      // Sum partial word.
      offset += __builtin_popcountl(occ_ptr->occ[wrdptr] >> (occ_ptr->occ_word_size - 1 - bit));
      // Returm sum.
      occ += offset;
   } else {
      int64_t  offset = 0;
      // Sum partial word.
      if (bit < occ_ptr->occ_word_size - 1)
	 offset += __builtin_popcountl(occ_ptr->occ[wrdptr] << (bit+1));
      // Sum bit offsets.
      for (uint64_t i = wrdptr + 1; i < wrdptr+(occ_ptr->occ_mark_intv-intv_wrd); i++)
	 offset += __builtin_popcountl(occ_ptr->occ[i]);
      // Return subtraction.
      occ -= offset;
   }
   return (size_t)occ_ptr->C[c] + occ;
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
   if (pos % csa->ratio == 0) {
      // Value is sampled. Extract it.
      size_t idx = pos / csa->ratio;
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
   size_t bwtpos = pos - (pos >= bwt->zero);
   uint8_t c = bwt->slots[bwtpos/4] >> (6-2*(bwtpos % 4)) & 0b11;
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
   // Beginning of reference
   if (range.bot <= bwt->zero && range.top >= bwt->zero) {
      sa_values[bwt->zero - range.bot] = path_offset;
   }
   // Get sampled SA values.
   // Compute offset to the closest %csa->ratio.
   size_t offset = (csa->ratio - range.bot%csa->ratio)%csa->ratio;
   size_t pos    = range.bot + offset;
   while (pos <= range.top) {
      // Extract sampled value.
      size_t idx = pos / csa->ratio;
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
      offset += csa->ratio;
      pos    += csa->ratio;
   }

   // If all positions are covered, return.
   int    done = 1;
   size_t rlen = range.top - range.bot + 1;
   for (int i = 0; i < rlen; i++) {
      if (sa_values[i] > occ->txtlen) {
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
      size_t bwtpos = pos - (pos >= bwt->zero);
      prev_c[i] = (bwt->slots[bwtpos/4] >> (6-2*(bwtpos % 4))) & 0b11;
      c_count[prev_c[i]]++;
      // 2. Mark the nucleotide for extension if SA value is missing.
      if (sa_values[i] > occ->txtlen) extend[(int)prev_c[i]] = 1;
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
   ssize_t nloc = range.top - range.bot + 1;
   size_t * sa_values = malloc(nloc * sizeof(size_t));
   exit_on_memory_error(sa_values);
   // Flag all positions as missing
   for (size_t i = 0; i < nloc; i++)
      sa_values[i] = occ->txtlen+1;

   // Do the recursive search
   recursive_csa_query(csa, bwt, occ, range, sa_values, 0);

   return sa_values;
}

char *
decompress_genome
(
 char * dna,
 size_t pos,
 size_t len
 )
{
   // Allocate output.
   char * seq = malloc(len+1);
   exit_on_memory_error(seq);

   for (size_t i = 0, p = pos; i < len; i++, p++)
      seq[i] = ALPHABET[(dna[p/4] >> (6-((p*2)%8))) & 0b11];

   seq[len] = 0;
   return seq;
}

char *
compress_genome
(
 char * genome,
 size_t gsize
 )
{
   // Allocate nucleotides.
   char * dna = calloc(gsize/4+2, 1);
   exit_on_memory_error(dna);

   for (size_t i = 0; i < gsize; i++)
      dna[i/4] |= (ENCODE[(int)genome[i]] & 0b11) << (6-((i*2)%8));

   // Store size of last byte (if 0, leave last byte empty and flag size=0)
   // this is to make it consistent with .pac file in bwa (to construct bwt).
   dna[gsize/4+1] = (uint8_t)(gsize%4);

   return dna;
}

void
pack_fasta
(
 const char * basename
)
{
   // Open files
   FILE * inputf = fzopen(basename, "r");
   if (inputf == NULL) exit_cannot_open(basename);
   
   char * fn = malloc(strlen(basename)+10);
   sprintf(fn, "%s.chr", basename);
   int fd = creat(fn, 0644);
   if (fd < 0) exit_cannot_open(fn);

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
   wstack_t * seqnames = stack_new(64);
   wstack_t * seqstart = stack_new(64);
   
   // Load fasta file line by line and concatenate.
   while ((rlen = getline(&buffer, &sz, inputf)) != -1) {
      if (buffer[0] == '>') {
	 buffer[rlen-1] = 0;
	 int k = 0;
	 while (buffer[k] != ' ' && buffer[k] != 0) k++;
	 buffer[k] = 0;
	 char * seqname = strdup(buffer+1);
	 push(seqname, &seqnames);
	 size_t * gpos = malloc(sizeof(size_t));
	 *gpos = gsize+1;
	 push(gpos, &seqstart);
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

   // Write chromosome names
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


   // FIXME //
   // Here we assume that if NULL is passed as a second parameter,
   // then the caller wants to read the index and we do nothing special.
   // be capitalized. Otherwise we assume the caller wants to create
   // the index and the genome must This is not good: we should pass a
   // flag to tell whether capitalization should be performed, or change
   // the logic of this part entirely.

   // Allocate packed nucleotides.
   char * dna = calloc((2*gsize)/4+2, 1);
   exit_on_memory_error(dna);

   // Pack forward strand
   ssize_t i = 0;
   for (; i < gsize ; i++) {
      uint8_t c = NONALPHABET[(uint8_t) genome[i]] ? 0 : ENCODE[(int)genome[i]];
      dna[i/4] |= (c & 0b11) << (6-((i*2)%8));
   }

   // Pack reverse strand
   for (ssize_t j = gsize-1; j >= 0; j--, i++) {
      uint8_t c = NONALPHABET[(uint8_t) genome[j]] ? 0 : ENCODE[(int)genome[j]];
      dna[i/4] |= ((3-c) & 0b11) << (6-((i*2)%8)); // Reverse complement (3-c)
   }

   // Store size of last byte (if 0, leave last byte empty and flag size=0)
   // this is to make it consistent with .pac file in bwa (to construct bwt).
   dna[(2*gsize)/4+1] = (uint8_t)((2*gsize)%4);

   // Write the compressed genome
   sprintf(fn, "%s.dna", basename);
   int fdna = creat(fn, 0644);
   if (fdna < 0) exit_cannot_open(fn);

   size_t bw = 0, bt = (2*gsize)/4+2;
   while (bw < bt) bw += write(fdna, dna + bw, bt - bw);
   close(fdna);

   // Clean up.
   free(fn);
   free(dna);
   free(buffer);
   free(genome);
}

wstack_t *
stack_new
(
 size_t max
 )
{
   size_t base = sizeof(wstack_t);
   size_t extra = max * sizeof(void *);
   wstack_t * stack = malloc(base + extra);
   exit_on_memory_error(stack);

   stack->max = max;
   stack->pos = 0;
   return stack;
}

void
push
(
 void      * ptr,
 wstack_t ** stackp
 )
{
   wstack_t * stack = *stackp;
   if (stack->pos >= stack->max) {
      size_t newmax = stack->max*2;
      stack = *stackp = realloc(stack,
				sizeof(wstack_t)+newmax*sizeof(void *));
      exit_on_memory_error(stack);
      stack->max = newmax;
   }

   stack->ptr[stack->pos++] = ptr;
   return;
}

u64stack_t *
u64stack_new
(
 size_t max
)
{
   size_t base = sizeof(u64stack_t);
   size_t extra = max * sizeof(uint64_t);
   u64stack_t * stack = malloc(base + extra);
   exit_on_memory_error(stack);

   stack->max = max;
   stack->pos = 0;
   return stack;
}

FILE *
fzopen
(
 const char * fname,
 const char * mode
)
{
   /*
    * Uses zcat pipe to decompress gzip files
    */

   FILE * inputf = fopen(fname, "r");
   uint8_t b0 = getc(inputf);
   uint8_t b1 = getc(inputf);
   exit_error(inputf == NULL);
   
   // Open file
   if ((b0 == 0x1f) && (b1 == 0x8b)) {
      fclose(inputf);
      // Open pipe for communication
      int pfd[2];
      int err = pipe(pfd);
      exit_error(err < 0);
   
      // Fork process
      int pid = fork();
      exit_error(pid < 0);
      if (pid) {
	 // parent
	 close(pfd[1]);
	 return fdopen(pfd[0], mode);
      } else {
	 // child
	 // dup pipe to stdout
	 close(pfd[0]);
	 close(1);
	 int newfd = dup(pfd[1]);
	 exit_error(newfd != 1);
	 close(pfd[1]);
	 // execv zcat
	 err = execlp("zcat","zcat",fname,NULL);
	 exit_error(err < 0);
	 return 0;
      }
   } else {
      ungetc(b1, inputf);
      ungetc(b0, inputf);
      return inputf;
   }
}
