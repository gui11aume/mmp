/*
 * divsufsort_private.h for libdivsufsort
 * Copyright (c) 2003-2008 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
# include "divsufsort.h"
# ifndef SAIDX_T
#  define SAIDX_T
#  define saidx_t saidx64_t
# endif /* SAIDX_T */
# ifndef PRIdSAIDX_T
#  define PRIdSAIDX_T PRIdSAIDX64_T
# endif /* PRIdSAIDX_T */

/*- Constants -*/
#if !defined(UINT8_MAX)
# define UINT8_MAX (255)
#endif /* UINT8_MAX */
#if defined(ALPHABET_SIZE) && (ALPHABET_SIZE < 1)
# undef ALPHABET_SIZE
#endif
#if !defined(ALPHABET_SIZE)
# define ALPHABET_SIZE (UINT8_MAX + 1)
#endif
/* for divsufsort.c */
#define BUCKET_A_SIZE (ALPHABET_SIZE)
#define BUCKET_B_SIZE (ALPHABET_SIZE * ALPHABET_SIZE)
/* for sssort.c */
#if defined(SS_INSERTIONSORT_THRESHOLD)
# if SS_INSERTIONSORT_THRESHOLD < 1
#  undef SS_INSERTIONSORT_THRESHOLD
#  define SS_INSERTIONSORT_THRESHOLD (1)
# endif
#else
# define SS_INSERTIONSORT_THRESHOLD (8)
#endif
#if defined(SS_BLOCKSIZE)
# if SS_BLOCKSIZE < 0
#  undef SS_BLOCKSIZE
#  define SS_BLOCKSIZE (0)
# elif 32768 <= SS_BLOCKSIZE
#  undef SS_BLOCKSIZE
#  define SS_BLOCKSIZE (32767)
# endif
#else
# define SS_BLOCKSIZE (1024)
#endif
/* minstacksize = log(SS_BLOCKSIZE) / log(3) * 2 */
#if SS_BLOCKSIZE == 0
#  define SS_MISORT_STACKSIZE (96)
#elif SS_BLOCKSIZE <= 4096
# define SS_MISORT_STACKSIZE (16)
#else
# define SS_MISORT_STACKSIZE (24)
#endif
#define SS_SMERGE_STACKSIZE (64)

/* for trsort.c */
#define TR_INSERTIONSORT_THRESHOLD (8)
#define TR_STACKSIZE (96)



/*- Macros -*/
#ifndef SWAP
# define SWAP(_a, _b) do { t = (_a); (_a) = (_b); (_b) = t; } while(0)
#endif /* SWAP */
#ifndef MIN
# define MIN(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#endif /* MIN */
#ifndef MAX
# define MAX(_a, _b) (((_a) > (_b)) ? (_a) : (_b))
#endif /* MAX */
#define STACK_PUSH(_a, _b, _c, _d)\
  do {\
    assert(ssize < STACK_SIZE);\
    stack[ssize].a = (_a), stack[ssize].b = (_b),\
    stack[ssize].c = (_c), stack[ssize++].d = (_d);\
  } while(0)
#define STACK_PUSH5(_a, _b, _c, _d, _e)\
  do {\
    assert(ssize < STACK_SIZE);\
    stack[ssize].a = (_a), stack[ssize].b = (_b),\
    stack[ssize].c = (_c), stack[ssize].d = (_d), stack[ssize++].e = (_e);\
  } while(0)
#define STACK_POP(_a, _b, _c, _d)\
  do {\
    assert(0 <= ssize);\
    if(ssize == 0) { return; }\
    (_a) = stack[--ssize].a, (_b) = stack[ssize].b,\
    (_c) = stack[ssize].c, (_d) = stack[ssize].d;\
  } while(0)
#define STACK_POP5(_a, _b, _c, _d, _e)\
  do {\
    assert(0 <= ssize);\
    if(ssize == 0) { return; }\
    (_a) = stack[--ssize].a, (_b) = stack[ssize].b,\
    (_c) = stack[ssize].c, (_d) = stack[ssize].d, (_e) = stack[ssize].e;\
  } while(0)
/* for divsufsort.c */
#define BUCKET_A(_c0) bucket_A[(_c0)]
#if ALPHABET_SIZE == 256
#define BUCKET_B(_c0, _c1) (bucket_B[((_c1) << 8) | (_c0)])
#define BUCKET_BSTAR(_c0, _c1) (bucket_B[((_c0) << 8) | (_c1)])
#else
#define BUCKET_B(_c0, _c1) (bucket_B[(_c1) * ALPHABET_SIZE + (_c0)])
#define BUCKET_BSTAR(_c0, _c1) (bucket_B[(_c0) * ALPHABET_SIZE + (_c1)])
#endif


/*- Private Prototypes -*/
/* sssort.c */
void
sssort(const sauchar_t *Td, const saidx_t *PA,
       saidx_t *first, saidx_t *last,
       saidx_t *buf, saidx_t bufsize,
       saidx_t depth, saidx_t n, saint_t lastsuffix);
/* trsort.c */
void
trsort(saidx_t *ISA, saidx_t *SA, saidx_t n, saidx_t depth);
