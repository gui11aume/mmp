#define _GNU_SOURCE
#include <ctype.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>


#ifndef _BWT_HEADER_
#define _BWT_HEADER_

// Size of the alphabet. Note that the code will break if
// the value is changed. It is indicated here in order
// to highlight where the alphabet size is important.

#define SIGMA 4

// Type aliases.
typedef struct blocc_t  blocc_t;
typedef struct bwt_t    bwt_t;
typedef struct csa_t    csa_t;
typedef struct lut_t    lut_t;
typedef struct occ_t    occ_t;
typedef struct range_t  range_t;

typedef unsigned int uint_t;

// Entries of the Occ table combine an Occ value (the cumulative
// number of occurrences of the character in the BWT up to and
// including the given position), and a bitfield where the value
// is set to 1 if and only if the BWT has the character at this
// position.
//
// Example 'blocc_t' data:
//
//  |  .smpl (uint32_t)  |          .bits (uint32_t)          |
//  |       146883       | 0b00100011100000000110000001000000 |
//
// The next .smpl value in an array of 'blocc_t' must be 146890
// (the popcount of the contiguous .bits value).
//
// To avoid confusion, we will use the term "index" when referring
// to 'blocc_t' structs of an array, and "position" when referring
// to the text or the BWT. Thus, the rank is taken on a position
// and not an index.
//
// There is one 'blocc_t' per 32 letters of the BWT. Since each
// 'blocc_t' occupies 64 bits, an array of 'blocc_t' occupies
// 2 bits per letter of the BWT. Since there is one array per
// symbol of the alphabet, an Occ table occupies '2*SIGMA' bits per
// letter of the BWT (i.e. one byte when 'SIGMA' is 4).
//
// Both .smpl and .bits are retrieved in a single memory reference
// (a 'blocc_t' occupies 8 bytes and cache lines are 64 bytes on
// x86), so both values are available for the price of a single
// cache miss.

struct blocc_t {
   uint_t smpl : 32;    // Sampled Occ values.
   uint_t bits : 32;    // Bitfield Occ values.
};

// Return type for range queries.
struct range_t {
   size_t bot;
   size_t top;
};

// The 'Occ_t' struct contains a size variable 'sz', followed by
// 'SIGMA' arrays of 'blocc_t', where 'SIGMA' is the number of
// letters in the alphabet. Note that 'sz' is not the number of
// 'blocc_t' in the arrays, but the size of the BWT, including the
// termination character.
struct occ_t {
   size_t   txtlen;      // 'strlen(txt) + 1'.
   size_t   C[SIGMA+1];  // The 'C' array.
   size_t   nrows;       // Number of entries.
   blocc_t  rows[0];     // Occ entries.
};

// The compressed suffix array.
struct csa_t {
   size_t    nbits;      // Bits in encoding.
   uint64_t  bmask;      // Bit mask (lower nbits set to 1).
   size_t    nint64;     // Size of the bit field.
   int64_t   bitf[0];    // Bie field.
};

// The Burrow-Wheeler transform.
struct bwt_t {
   size_t   txtlen;      // 'strlen(txt) + 1'.
   size_t   zero;        // Position of '$'.
   size_t   nslots;      // Number of slots.
   uint8_t  slots[0];    // 2-bit characters.
};

// Lookup table 
#define LUTK 12  // Size of k-mers in the LUT.
struct lut_t { range_t kmer[1<<(2*LUTK)]; };


// VISIBLE FUNCTION DECLARATIONS //

csa_t   * compress_sa(int64_t *);
int64_t * compute_sa(const char *);
bwt_t   * create_bwt(const char *, const int64_t *);
occ_t   * create_occ(const bwt_t *);
void      fill_lut(lut_t *, const occ_t *, range_t, size_t, size_t);
char    * normalize_genome(FILE *);

#endif
/*

// ------- Error handling macros ------- //

#define exit_if_null(x) \
   do { if ((x) == NULL) { fprintf(stderr, "memory error %s:%d:%s()\n", \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); }} while(0)

#define exit_cannot_open(x) \
   do { fprintf(stderr, "cannot open file '%s' %s:%d:%s()\n", (x), \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); } while(0)

#define exit_if(x) \
   do { if (x) { fprintf(stderr, "%s %s:%d:%s()\n", #x, \
         __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); }} while(0)


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
*/
