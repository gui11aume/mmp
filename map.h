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

typedef struct aln_t      aln_t;
typedef struct alnstack_t alnstack_t;
typedef struct mem_t      mem_t;

struct mem_t {
   size_t    beg;
   size_t    end;
   range_t   range;
   size_t  * sa;
   int       aligned;
};

struct aln_t {
   int          score;
   int          nmem;
   size_t       refpos;
   const char * refseq;
   mem_t        mem;
};

struct alnstack_t {
   size_t pos;
   size_t max;
   aln_t  aln[];
};


alnstack_t * mapread (const char *, index_t, const char *, size_t);

#endif
