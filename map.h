#define _GNU_SOURCE

#include <assert.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "bwt.h"

#ifndef _map_HEADER_
#define _map_HEADER_

typedef struct alnstack_t alnstack_t;
typedef struct aln_t      aln_t;
typedef struct align_t    align_t;
typedef struct seed_t     seed_t;

struct seed_t {
         size_t    beg;
         size_t    end;
         range_t   range;
         size_t  * sa;
         int       aligned;
};

struct aln_t {
   int      score;
   size_t   refpos;
   char   * refseq;
   size_t   read_beg;
   size_t   read_end;
   double   qual;
};

struct align_t {
   size_t    refpos;
   size_t    span;
   int       minscore;
   seed_t  * seed;
};

struct alnstack_t {
         size_t pos;
         size_t max;
         aln_t  aln[];
};


alnstack_t * mapread (const char *, const index_t, const size_t, const int, const int);
void         extend_l1l2 (const char *, const index_t, seed_t *, seed_t *);
int          filter_longest_mem (wstack_t *);
alnstack_t * alnstack_new (size_t max);
void         align (align_t , const char*, char *, size_t , int *, alnstack_t **);

#endif
