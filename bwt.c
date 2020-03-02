#define _GNU_SOURCE
#include <assert.h>
#include <unistd.h>
#include <stdio.h>
#include "divsufsort_private.h"
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

void
compute_sa
(
 const char * txt,
 const char * basename
 )
{
   const size_t txtlen = strlen(txt) + 1;
   // Allocate memory for SA
   saidx64_t * array  = calloc(txtlen, sizeof(saidx64_t));
   exit_on_memory_error(array);
   // Compute SA
   divsufsort((const unsigned char *) txt, array, txtlen);

   // Store in disk
   char * fname = malloc(strlen(basename)+100);
   exit_on_memory_error(fname);
   sprintf(fname,"%s.sa.tmp", basename);
   int sa_fd = creat(fname, 0644);

   // Write to file descriptor
   size_t bw = 0, bt = txtlen*sizeof(saidx64_t);
   char * sa_bytes = (char *)array;
   while(bw < bt) bw += write(sa_fd, sa_bytes + bw, bt-bw);

   close(sa_fd);
   free(array);
   free(fname);
   return;
}

int
fill_buffer
(
 u64stack_t * buffer,
 int fd
 )
{
   ssize_t bytes, wr = 0, wt = buffer->max * sizeof(uint64_t);
   while (wr<wt) {
      bytes = read(fd, ((char *)buffer->val)+wr, wt-wr);
      if (bytes < 1) {
	 buffer->max = wr/sizeof(uint32_t);
	 break;
      }
      wr += bytes;
   }
   buffer->pos = 0;
   return buffer->max == 0;
}

int
compare
(
 const char  * str_a,
 const char  * str_b
)
{
   return strcmp(str_a, str_b);
}

void
compute_fmindex
(
 const char * basename,
 const char * txt
 )
{
   const size_t txtlen = strlen(txt) + 1; // Do not forget $.
   ssize_t BUFSIZE = min(txtlen*8, 125000000); // Max 1 GB
   size_t zero = 0;
   ssize_t bw = 0, bt = 0;

   // Allocate buffers
   u64stack_t * sa_buf = u64stack_new(BUFSIZE);

   // Open file descriptors
   char * fname = malloc(strlen(basename)+100);
   exit_on_memory_error(fname);
   sprintf(fname, "%s.sa.tmp", basename);
   int sa_fd = open(fname, O_RDONLY);

   // Open CSA stream
   size_t nbits = 0;
   while (txtlen > ((uint64_t) 1 << nbits)) nbits++;
   size_t nnumb = (txtlen + (CSA_SAMP_RATIO-1)) / CSA_SAMP_RATIO;
   size_t nint64 = (nbits * nnumb + (64-1)) / 64;
   uint64_t bmask = ((uint64_t) 0xFFFFFFFFFFFFFFFF) >> (64-nbits);
   
   sprintf(fname,"%s.sa", basename);
   int csa_fd = creat(fname, 0644);

   // struct header
   size_t csa_ratio = CSA_SAMP_RATIO;
   bw = write(csa_fd, &nbits, sizeof(size_t));
   bw = write(csa_fd, &bmask, sizeof(uint64_t));
   bw = write(csa_fd, &nint64, sizeof(size_t));
   bw = write(csa_fd, &csa_ratio, sizeof(size_t));

   // Open BWT stream
   const size_t nslots = (txtlen + (4-1)) / 4;

   sprintf(fname,"%s.bwt", basename);
   int bwt_fd = creat(fname, 0644);

   // struct header 
   bw = write(bwt_fd, &txtlen, sizeof(size_t));
   bw = write(bwt_fd, &txtlen, sizeof(size_t)); // bwt_zero
   bw = write(bwt_fd, &nslots, sizeof(size_t));

   // Open OCC stream
   uint64_t word_size = 64;
   uint64_t mark_bits = OCC_INTERVAL_SIZE * word_size;
   const uint64_t nintv = (txtlen + mark_bits - 1) / mark_bits;
   const uint64_t nword = nintv * OCC_INTERVAL_SIZE;
   const uint64_t nmark = nintv + 1;
   const size_t occ_size = nword+nmark;

   sprintf(fname, "%s.occ", basename);
   int occ_fd = creat(fname, 0644);
   
   // Struct header
   uint64_t occ_mark_intv = OCC_INTERVAL_SIZE;
   bw = write(occ_fd, &txtlen, sizeof(size_t));
   bw = write(occ_fd, &occ_size, sizeof(uint64_t));
   bw = write(occ_fd, &occ_mark_intv, sizeof(uint64_t));
   bw = write(occ_fd, &word_size, sizeof(uint64_t));
   bw = write(occ_fd, &mark_bits, sizeof(uint64_t));
   for (int i = 0; i <= SIGMA; i++) //C
      bw = write(occ_fd, &zero, sizeof(uint64_t));

   // General loop variables
   // CSA loop variables
   uint8_t csa_lastbit = 0;
   uint64_t csa_word = 0;
   // BWT loop variables
   int64_t bwt_zero;
   size_t bwt_bufsize = 4096;
   uint8_t * bwt_bytes = malloc(bwt_bufsize);
   // OCC loop variables
   int intv_words = SIGMA*(1+OCC_INTERVAL_SIZE);
   uint64_t * occ_intv = calloc(intv_words, sizeof(uint64_t));
   exit_on_memory_error(occ_intv);
   uint64_t occ_abs[SIGMA] = {0};

   size_t i = 0;
   sa_buf->pos = sa_buf->max; // Force buffer filling.
   for (; i < txtlen; i++) {
      if (sa_buf->pos == sa_buf->max)
	 fill_buffer(sa_buf, sa_fd);
      uint64_t saidx = sa_buf->val[sa_buf->pos++];

      // CSA (1:csa_samp_rate compression)
      if (i%CSA_SAMP_RATIO == 0) {
	 // Shift first bits of saidx
	 csa_word |= saidx << csa_lastbit;
	 // Update bit offset
	 csa_lastbit += nbits;
      
	 if (csa_lastbit >= 64) {
	    // Shift bit position
	    csa_lastbit = csa_lastbit - 64;
	    // Write word
	    bw = write(csa_fd, &csa_word, sizeof(uint64_t));
	    // Start new word
	    csa_word = 0;
	    // Write last bits of saidx
	    csa_word |= saidx >> (nbits - csa_lastbit);
	 }
      }
      // BWT & OCC
      uint8_t c = 0;
      if (saidx != 0) {
	 c = ENCODE[(uint8_t) txt[saidx-1]];
	 occ_abs[c]++;
	 occ_intv[SIGMA + c*OCC_INTERVAL_SIZE + (i%mark_bits)/word_size] |= ((uint64_t)1) << (word_size - 1 - (i % word_size));
      } else {
	 bwt_zero = i;
      }
      
      bwt_bytes[(i/4)%bwt_bufsize] |= c << 2*(i % 4);

      // Write BWT word
      if ((i+1) % (4*bwt_bufsize) == 0) {
	 bw = 0; bt = bwt_bufsize;
	 while(bw < bt) bw += write(bwt_fd, bwt_bytes + bw, bt-bw);
	 memset(bwt_bytes, 0, bwt_bufsize);
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

      // Verbose
      if (i % 100000 == 0)
	 fprintf(stderr, "\r%.2f%%", (100.0*i)/txtlen);
   }

   // Finish CSA
   bw = write(csa_fd, &csa_word, sizeof(uint64_t));
   
   // Finish BWT
   if (i % (4*bwt_bufsize) > 0) {
      bw = 0; bt = (i%(4*bwt_bufsize))/4 + ((i%4) > 0);
      while(bw<bt) bw += write(bwt_fd, bwt_bytes + bw, bt-bw);
   }
   
   // Seek back in bwt_fd to write zero position
   lseek(bwt_fd, sizeof(size_t), SEEK_SET);
   bw = write(bwt_fd, &bwt_zero, sizeof(size_t));
   
   // Finish OCC
   if (i % mark_bits > 0) {
      // Write full interval
      ssize_t bw = 0, bt = intv_words*sizeof(uint64_t);
      while(bw<bt) bw += write(occ_fd, ((char *)occ_intv)+bw, bt-bw);
   }
   // Write last OCC absolute position mark
   bw = 0; bt = SIGMA*sizeof(uint64_t);
   while(bw<bt) bw += write(occ_fd, ((char *)occ_abs)+bw, bt-bw);
   
   // Compute 'C'.
   uint64_t C[SIGMA+1];
   C[0] = 1;
   for (int i = 1 ; i < SIGMA+1 ; i++)
      C[i] = C[i-1] + occ_abs[i-1];

   // Seek back in occ_fd to write C table
   lseek(occ_fd, sizeof(size_t) + 4*sizeof(uint64_t), SEEK_SET);
   bw = 0; bt = (SIGMA+1)*sizeof(uint64_t);
   while (bw<bt) bw += write(occ_fd, ((char *)C)+bw, bt-bw);
   
   // Close streams
   close(csa_fd);
   close(bwt_fd);
   close(occ_fd);
   close(sa_fd);

   // Remove tmp file
   sprintf(fname, "%s.sa.tmp", basename);
   remove(fname);

   // Free memory
   free(sa_buf);
   free(bwt_bytes);
   free(fname);
   free(occ_intv);
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
      if (sa_values[i] > bwt->txtlen) {
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
      if (sa_values[i] > bwt->txtlen) extend[(int)prev_c[i]] = 1;
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
      sa_values[i] = bwt->txtlen+1;

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
      seq[i] = ALPHABET[(dna[p/4] >> ((p*2)%8)) & 0b11];

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
   char * dna = calloc(gsize/4+1, 1);
   exit_on_memory_error(dna);

   for (size_t i = 0; i < gsize; i++)
      dna[i/4] |= (ENCODE[(int)genome[i]] & 0b11) << ((i*2)%8);

   return dna;
}

char *
normalize_genome
(
 FILE   * inputf,
 char   * chrfile,
 size_t * gsize_ptr
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


   // FIXME //
   // Here we assume that if NULL is passed as a second parameter,
   // then the caller wants to read the index and we do nothing special.
   // be capitalized. Otherwise we assume the caller wants to create
   // the index and the genome must This is not good: we should pass a
   // flag to tell whether capitalization should be performed, or change
   // the logic of this part entirely.
  
   // Normalize (use only capital alphabet letters).
   for (size_t pos = 0; pos < gsize ; pos++) {
      uint64_t iter = 0;
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

   // Return genome size.
   *gsize_ptr = 2*div;

   // Clean up.
   free(buffer);

   return genome;

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
