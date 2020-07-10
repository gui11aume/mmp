#include "map.h"
#define max(x,y) ((x) > (y) ? (x) : (y))
#define min(x,y) ((x) < (y) ? (x) : (y))

typedef struct seedchain_t seedchain_t;
typedef struct aligncd_t   aligncd_t;

#define MAX_MEM_SEED_LOCI 300
#define MAX_SKIP_CHAIN_SEEDS 10000
#define MAX_MINSCORE_REPEATS 2
#define MAX_ALIGN_MISMATCHES 30
#define MAX_CHAIN_INDEL_RATE 0.1

#ifdef DEBUG
int DEBUG_VERBOSE = 1;
#else
int DEBUG_VERBOSE = 0;
#endif

#define min(x,y) ((x) < (y) ? (x) : (y))
#define min3(x,y,z) (min(min(x,y),z))

struct seedchain_t {
   size_t    pos;
   size_t    max;
   size_t    loci;
   int       minscore;
   int       span;
   seed_t  * seed[];
};

struct aligncd_t {
   size_t     cnt;
   align_t  * align;
};

seedchain_t * seedchain_new (size_t max);
void          seed_push     (seed_t * mem, seedchain_t ** stackp);
void          aln_push     (aln_t aln, alnstack_t ** stackp);

int seed_by_refpos (const void * a, const void * b) {
   return ((align_t *)a)->refpos > ((align_t *)b)->refpos;
};
int seed_by_span (const void * a, const void * b) {
   align_t * sa = (align_t *)a;
   align_t * sb = (align_t *)b;
   return sa->span < sb->span;
};

int seed_by_start (const void * a, const void * b) {
   return (*(seed_t **)a)->beg > (*(seed_t **)b)->beg;
};

int mem_by_loci (const void * a, const void * b) {
   seed_t * sa = *(seed_t **)a;
   seed_t * sb = *(seed_t **)b;
   return (sa->range.top - sa->range.bot) > (sb->range.top - sb->range.bot);
};

int mem_by_span (const void * a, const void * b) {
   seed_t * sa = *(seed_t **)a;
   seed_t * sb = *(seed_t **)b;
   return (sa->end - sa->beg) < (sb->end - sb->beg);
};

int seed_by_first_locus (const void * a, const void * b) {
   const seed_t A = **(seed_t **) a;
   const seed_t B = **(seed_t **) b;
   return (A.sa[0] > B.sa[0]) - (A.sa[0] < B.sa[0]);
}

int SA_by_locus (const void *a, const void *b) {
   const size_t locus_A = *(size_t *)a;
   const size_t locus_B = *(size_t *)b;
   return (locus_A > locus_B) - (locus_A < locus_B);
}


int           minscore_then_span (const void * a, const void * b) {
   int sa = (*(seedchain_t **)a)->minscore;
   int sb = (*(seedchain_t **)b)->minscore;
   if (sa > sb)
      return 1;
   else if (sa < sb)
      return -1;
   else
      return (*(seedchain_t **)a)->span < (*(seedchain_t **)b)->span;
};

int           align_minscore_then_span (const void * a, const void * b) {
   int sa = ((align_t *)a)->minscore;
   int sb = ((align_t *)b)->minscore;
   if (sa > sb)
      return 1;
   else if (sa < sb)
      return -1;
   else
      return ((align_t *)a)->span < ((align_t *)b)->span;
};

// Error-handling macros.
#define exit_on_memory_error(x)						\
   do { if ((x) == NULL) { fprintf(stderr, "memory error %s:%d:%s()\n", \
				   __FILE__, __LINE__, __func__); exit(EXIT_FAILURE); }} while(0)




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


void
recursive_mem_chain
(
   wstack_t  * mems,
   size_t      mem_pos,
   size_t      chain_pos,
   seed_t   ** chain,
   wstack_t ** chain_stack
)
{
   size_t mem_end = ((seed_t *) mems->ptr[mem_pos])->end;
   // For MEMs overlapping MEM[pos].
   for (size_t i = mem_pos; i < mems->pos && ((seed_t *)mems->ptr[i])->beg <= mem_end; i++) {
      // Extend chain.
      chain[chain_pos] = (seed_t *) mems->ptr[i];
      // Get next nonoverlapping MEM.
      size_t j;
      for (j = i+1; j < mems->pos; j++) {
         if (((seed_t *)mems->ptr[j])->beg > ((seed_t *) mems->ptr[i])->end)
            break;
      }
      // We reached the chain end, store chain.
      if (j >= mems->pos) {
         seedchain_t * seedchain = seedchain_new(chain_pos+1);
         // Push mems to chain and compute span.
         int span = 0;
         size_t loci = 0;
         for (size_t k = 0; k <= chain_pos; k++) {
            span += chain[k]->end - chain[k]->beg + 1;
            loci += chain[k]->range.top - chain[k]->range.bot + 1;
            seed_push(chain[k], &seedchain);
         }
         // Push mem chain to chain stack.
         seedchain->span = span;
         seedchain->loci = loci;
         push(seedchain, chain_stack);
      } else {
         recursive_mem_chain(mems, j, chain_pos+1, chain, chain_stack);
      }
   }
}


seed_t *
new_one_locus_seed
(
   seed_t * old,
   size_t   idx
)
{
   seed_t * new = calloc(1, sizeof(seed_t));
   exit_on_memory_error(new);
   range_t range = (range_t) {
      // One locus.
      .bot = old->range.bot + idx,
      .top = old->range.bot + idx,
   };
   new->beg = old->beg;
   new->end = old->end;
   new->range = range;
   new->sa = malloc(sizeof(size_t));
   exit_on_memory_error(new->sa);
   new->sa[0] = old->sa[idx];
   return new;
}


wstack_t *
merge_overlapping_seeds
(
         wstack_t * skip_seeds,
   const index_t    idx
)
{

   // Initialize stack of merged seeds.
   wstack_t * merged = stack_new(64);
   exit_on_memory_error(merged);

   // Get SA values and sort by locus.
   for (int i = 0 ; i < skip_seeds->pos ; i++) {
      seed_t * s = (seed_t *) skip_seeds->ptr[i];
      size_t nloci = s->range.top - s->range.bot + 1;
      if (nloci > 16) {
         // Keep at most 16 loci per seed.
         s->range.top = s->range.bot + 15;
         nloci = 16;
      }
      s->sa = query_csa_range(idx.csa, idx.bwt, idx.occ, s->range);
      // Sort SA values by locus order.
      qsort(s->sa, nloci, sizeof(size_t), SA_by_locus);
   }

   // Temporary arrays for the merge.
   seed_t *  databuffer_1[128] = {0};
   seed_t *  databuffer_2[128] = {0};
   seed_t ** prev = databuffer_1;
   seed_t ** next = databuffer_2;

   // Assume seeds are sorted by start/end position.
   // Start with the leftmost seed and create new one-locus seeds.
   seed_t * leftmost_seed = (seed_t *) skip_seeds->ptr[skip_seeds->pos-1];
   int n_prev_loci = leftmost_seed->range.top - leftmost_seed->range.bot + 1;
   if (n_prev_loci > 128) n_prev_loci = 128;
   for (int n = 0 ; n < n_prev_loci ; n++) {
      // Calls 'malloc' and 'exit_on_memory_error'.
      seed_t * seed = new_one_locus_seed(leftmost_seed, n);
      push(seed, &merged);
      prev[n] = seed;
   }

   for (int i = skip_seeds->pos-2 ; i >= 0 ; i--) {
      seed_t * next_seed = (seed_t *) skip_seeds->ptr[i];
      int n_next_loci = next_seed->range.top - next_seed->range.bot + 1;
      // Left and right pointers.
      int lidx = 0;
      int ridx = 0;
      int inc = 0;
      while (inc < 128 && lidx < n_prev_loci && ridx < n_next_loci) {
         seed_t * prev_seed = prev[lidx];
         ssize_t next_locus = next_seed->sa[ridx];
         const ssize_t dhit = next_seed->beg - prev_seed->beg;
         if (next_locus > prev_seed->sa[0] + dhit) {
            // The one-locus seed already exists.
            if (prev_seed->end >= next_seed->beg) {
               next[inc++] = prev_seed;
            }
            lidx++;
         }
         else if (next_locus < prev_seed->sa[0] + dhit) {
            // Create new one-locus seed.
            seed_t * new = new_one_locus_seed(next_seed, ridx);
            push(new, &merged);
            next[inc++] = new;
            ridx++;
         }
         else {
            // Found a match: extend the one-locus seed.
            prev_seed->end = next_seed->end;
            next[inc++] = prev_seed;
            lidx++;
            ridx++;
         }
      }
      // At most one of the loop I or II is executed.
      // Add remaining loci or seed / loci.
      for ( ; inc < 128 && lidx < n_prev_loci ; lidx++) {
         // Loop I.
         seed_t * seed = prev[lidx];
         if (seed->end >= next_seed->beg) {
            next[inc++] = seed;
         }
      }
      for ( ; inc < 128 && ridx < n_next_loci ; ridx++) {
         // Loop II.
         seed_t * new = new_one_locus_seed(next_seed, ridx);
         exit_on_memory_error(new);
         push(new, &merged);
         next[inc++] = new;
      }
      // Update prev / next.
      seed_t ** tmp = prev;
      prev = next;
      next = tmp;
      n_prev_loci = inc;
   }

   // The 'merged' stack contains the seed / loci.
   return merged;

}

wstack_t *
nonoverlapping_mems
(
   wstack_t * mems
)
{
   // Alloc.
   wstack_t * chain_stack = stack_new(8);

   if (mems->pos > 0) {
      seed_t  ** chain = malloc(mems->pos * sizeof(seed_t *));
      exit_on_memory_error(chain);

      // 1. Sort MEMs by start position.
      qsort(mems->ptr, mems->pos, sizeof(seed_t *), seed_by_start);

      // 2. Recursive call to mem group.
      recursive_mem_chain(mems, 0, 0, chain, &chain_stack);

      free(chain);
   }

   // 3. Return stack of non-overlapping MEM combinations.
   return chain_stack;
}


int
mem_chain_min_score
(
   seedchain_t * chain,
   const int    seqlen
)
{
   // Commented lines remove the MEM masking bug.
   // The code was not removed because the idea can
   // be reused with fixed-length seeds.
   int minscore = 0;

   // Add mismatches at chain ends.

   if (chain->seed[0]->beg > 0)
      //minscore += max(0,chain->mem[0]->beg / gamma - 1) + 1;
      minscore += 1;

   if (chain->seed[chain->pos-1]->end < seqlen-1)
      //minscore += max(0,(seqlen - 2 - chain->mem[chain->pos-1]->end)/gamma - 1) + 1;
      minscore += 1;

   // Add gap mismatches.
   for (int i = 1; i < chain->pos; i++) {
      // A gap implies one mismatch, even if it's a gap of length 0.
      //int gap_size = chain->mem[i]->beg - chain->mem[i-1]->end - 1;
      //minscore += max(0,gap_size / gamma - 1) + 1;
      minscore += 1;
   }

   return minscore;
}

wstack_t *
chain_mems
(
 int        slen,
 wstack_t * mems
)
{
   // Find all non-overlapping MEM combinations.
   wstack_t * chain_stack = nonoverlapping_mems(mems);

   // Compute minimum alignment score given seed distribution.
   for (int i = 0; i < chain_stack->pos; i++) {
      seedchain_t * chain = (seedchain_t *)chain_stack->ptr[i];
      chain->minscore = mem_chain_min_score(chain, slen);
   }

   // Sort mem chains by minscore(inc) then span(dec).
   qsort(chain_stack->ptr, chain_stack->pos, sizeof(seedchain_t *), minscore_then_span);

   // DEBUG VERBOSE
   if (DEBUG_VERBOSE) {
      fprintf(stdout,"MEM chains (%ld):\n", chain_stack->pos);
      for(int i = 0; i < chain_stack->pos; i++) {
	 seedchain_t * c = (seedchain_t *) chain_stack->ptr[i];
	 fprintf(stdout, "Chain [%d] (mems: %ld, span: %d, minscore: %d):\n", i, c->pos, c->span, c->minscore);
	 for (int j = 0; j < 0; j++) {
	    seed_t * m = c->seed[j];
	    fprintf(stdout, "[%d] (%ld, %ld) range: (%ld, %ld)\n", i, m->beg, m->end, m->range.bot, m->range.top);
	 }
      }
      fprintf(stdout,"\n");
   }

   return chain_stack;
}

aligncd_t
chain_skip
(
 size_t     slen,
 int        gamma,
 int        skip,
 wstack_t * seeds,
 index_t    idx
)
{
   // Get all Suffix Arrays.
   size_t nloc = 0;
   for (size_t i = 0; i < seeds->pos; i++) {
      seed_t * seed = (seed_t *)seeds->ptr[i];
      size_t seed_loc = seed->range.top - seed->range.bot + 1;
      if (seed_loc > MAX_SKIP_CHAIN_SEEDS)
	 return (aligncd_t){0, NULL};
      nloc += seed_loc;
   }

   if (nloc == 0)
      return (aligncd_t){0, NULL};

   // Recompute seed positions
   ssize_t     b = slen - gamma;
   size_t   nbeg = b/skip + 1 + (b%skip > 0);
   size_t * sbeg = malloc(nbeg * sizeof(size_t));
   exit_on_memory_error(sbeg);
   for (ssize_t i = nbeg-1; i >= 0; i--) {
      sbeg[i] = b;
      b = max(b - skip, 0);
   }

   // Allocate all suffix array positions
   align_t * loc_list = malloc(nloc*sizeof(align_t));
   exit_on_memory_error(loc_list);

   // Make chained alignment candidates from seed genomic positions
   size_t j = 0;
   for (size_t i = 0; i < seeds->pos; i++) {
      seed_t * seed = (seed_t *)seeds->ptr[i];
      // Get genomic positions
      seed->sa = query_csa_range(idx.csa, idx.bwt, idx.occ, seed->range);
      for (int k = 0; k < seed->range.top - seed->range.bot + 1; k++) {
	 loc_list[j++] = (align_t){seed->sa[k], seed->beg, 1, seed};
      }
   }

   // Sort loci
   qsort(loc_list, j, sizeof(align_t), seed_by_refpos);

   // Chain alignment positions
   int   max_indels = min(skip-1, slen*(MAX_CHAIN_INDEL_RATE));
   size_t    nchain = 0;
   wstack_t * chain = stack_new(seeds->pos);

   for (size_t n = 0; n < nloc; n++) {
      // Skip consumed seeds
      if (loc_list[n].minscore == -1)
	 continue;

      // Append seed to chain
      chain->pos = 0;
      push(loc_list[n].seed, &chain);

      // Find chain
      int max_ref_dist = (slen - loc_list[n].span - gamma + 1 + max_indels);
      ssize_t read_last = 0;
      for (size_t j = n+1; j < nloc; j++) {
	 int gen_dist = ((ssize_t)loc_list[j].refpos - (ssize_t)loc_list[n].refpos);

	 // No more possible chaining
	 if (gen_dist > max_ref_dist)
	    break;

	 // Seed mislocation
	 if (loc_list[n].span >= loc_list[j].span || loc_list[j].span <= read_last)
	    continue;

	 // Compute distance between seeds
	 int read_dist = (ssize_t)loc_list[j].span - (ssize_t)loc_list[n].span;

	 // Chain gap too big
	 if (gen_dist > read_dist + max_indels || gen_dist < read_dist - max_indels)
	    continue;

	 // Append seed to chain
	 push(loc_list[j].seed, &chain);
	 read_last = loc_list[j].span;

	 // Mark seed as consumed
	 loc_list[j].minscore = -1;
      }

      // Chain min score
      size_t  gap = 0;
      size_t  pos = 0;
      size_t span = 0;
      int minscore = 0;
      for (size_t i = 0; i < chain->pos; i++) {
	 seed_t * s = (seed_t *)chain->ptr[i];
	 // Gap found
	 if (s->beg > gap) {
	    // Set pos to first mismatched seed
	    while (sbeg[pos] + gamma <= gap)
	       pos++;
	    // Kill minimum number of seeds to produce this gap
	    while (sbeg[pos] < s->beg) {
	       minscore++;
	       // Skip overlapping seeds
	       size_t seedend = sbeg[pos] + gamma;
	       while (pos < nbeg && sbeg[pos] < s->beg && sbeg[pos] < seedend)
		  pos++;
	    }
	    // Span
	    span += gamma;
	 } else {
	    span += s->end - gap + 1;
	 }
	 gap = s->end+1;
      }

      // Trailing gap
      if (gap < slen) {
	 while (sbeg[pos] + gamma <= gap)
	    pos++;
	 while (pos < nbeg) {
	    minscore++;
	    // Skip overlapping seeds
	    size_t seedend = sbeg[pos] + gamma;
	    while (pos < nbeg && sbeg[pos] < seedend)
	       pos++;
	 }
      }

      // Create alignment
      loc_list[nchain++] = (align_t){loc_list[n].refpos, span, minscore, loc_list[n].seed};
   }

   free(chain);
   free(sbeg);

   loc_list = realloc(loc_list, nchain*sizeof(align_t));
   exit_on_memory_error(loc_list);

   // Sort by minscore then span.
   qsort(loc_list, nchain, sizeof(align_t), align_minscore_then_span);

   return (aligncd_t){nchain, loc_list};

}

void
extend_L1L2
(
 const char   * seq,
 const int      len,
 const index_t  idx,
 seed_t       * L1,
 seed_t       * L2
)
{

   // Preallocate ranges
   range_t newrange = {0};

   int i;
   size_t merid;

   // L1.
   L1->end = len-1;

   // Look up the beginning (reverse) of the query in lookup table.
   merid = 0;
   for (int j = 0 ; j < LUTK ; j++) {
      // Note: every "N" is considered "A".
      uint8_t c = ENCODE[(uint8_t) seq[len-j-1]];
      merid = c + (merid << 2);
   }

   range_t range = idx.lut->kmer[merid];

   if (range.top < range.bot) {
      range = (range_t) {.bot = 1, .top = idx.occ->txtlen-1};
      i = len-1;
   }
   else {
      i = len-1 - LUTK;
   }

   for ( ; i >= 0 ; i--) {
      if (NONALPHABET[(uint8_t)seq[i]]) break;
      int c = ENCODE[(uint8_t) seq[i]];
      newrange.bot = get_rank(idx.occ, c, range.bot - 1);
      newrange.top = get_rank(idx.occ, c, range.top) - 1;
      // Stop if fewer than 2 hits.
      if (newrange.top < newrange.bot + 1)
         break;
      range = newrange;
   }

   L1->beg = i+1;
   L1->range = range;

   // Check L1 result
   if (L1->end - L1->beg + 1 == len) {
      L2->beg   = L1->beg;
      L2->end   = L1->end;
      L2->range = L1->range;
      return;
   }

   // L2.
   L2->beg = 0;

   // Look up the beginning (forward) of the query in lookup table.
   merid = 0;
   for (int j = 0 ; j < LUTK ; j++) {
      // Note: every "N" is considered "A".
      uint8_t c = REVCMP[(uint8_t) seq[j]];
      merid = c + (merid << 2);
   }

   range = idx.lut->kmer[merid];

   if (range.top < range.bot) {
      range = (range_t) {.bot = 1, .top = idx.occ->txtlen-1};
      i = 0;
   }
   else {
      i = LUTK;
   }

   for ( ; i < len ; i++) {
      if (NONALPHABET[(uint8_t)seq[i]]) break;
      int c = REVCMP[(uint8_t) seq[i]];
      newrange.bot = get_rank(idx.occ, c, range.bot - 1);
      newrange.top = get_rank(idx.occ, c, range.top) - 1;
      // Stop if fewer than 2 hist.
      if (newrange.top < newrange.bot + 1)
	 break;
      range = newrange;
   }

   L2->end = i-1;
   L2->range = range;

   return;

}

wstack_t *
mem_seeds
(
 const char    * seq,
 const index_t   idx,
 const size_t    gamma
)
{
   int len = strlen(seq);
   int end = len-1;
   while (end > 0 && NONALPHABET[(uint8_t)seq[end]]) end--;

   if (end == 0)
      return stack_new(1);

   // Initialize mem stack
   wstack_t * mems = stack_new(32);

   range_t range;
   range_t newrange = {0};

   // Iterate over all read positions
   while (1) {
      seed_t mem = {0};
      mem.end = end;

      // Backward <<<
      range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
      int mpos = end, mlen = 0;

      // Query the beginning of the read in lookup table.
      if (end >= LUTK - 1) {
         size_t merid = 0;
         for ( ; mlen < LUTK ; mlen++, mpos--) {
            if (NONALPHABET[(uint8_t) seq[end-mlen]]) {
               range.bot = 1;
               range.top = 0;
               break;
            }
            uint8_t c = ENCODE[(uint8_t) seq[end-mlen]];
            merid = c + (merid << 2);
         }
         range = idx.lut->kmer[merid];
      }

      // Cancel if we went too far already.
      if (range.top < range.bot) {
         range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
         mpos = end; mlen = 0;
      }

      for ( ; mpos >= 0 ; mpos--, mlen++) {
         if (NONALPHABET[(uint8_t)seq[mpos]]) break;
         int c = ENCODE[(uint8_t) seq[mpos]];
         newrange.bot = get_rank(idx.occ, c, range.bot - 1);
         newrange.top = get_rank(idx.occ, c, range.top) - 1;
         // Stop if no hits.
         if (newrange.top < newrange.bot)
            break;
         range = newrange;
      }

      mem.beg = ++mpos;
      mem.range = range;

      // Keep MEM if above minimum length.
      if (mlen >= gamma) {
         seed_t * m = malloc(sizeof(seed_t));
         exit_on_memory_error(m);
         memcpy(m, &mem, sizeof(seed_t));
         push(m, &mems);
      }

      if (mem.beg < 1) break;

      // Find new end position (forward).
      end = mem.beg - 1;
      if (NONALPHABET[(uint8_t) seq[end]]) {
         while (end > 0 && NONALPHABET[(uint8_t) seq[end]]) end--;
      } else {
         range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
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
      }

      if (end + 1 < gamma) break; // No more seeds.

   }

   // DEBUG VERBOSE
   if (DEBUG_VERBOSE) {
      fprintf(stdout,"\nMEMs (%ld):\n", mems->pos);
      for(int i = 0; i < mems->pos; i++) {
         seed_t * m = (seed_t *) mems->ptr[i];
         fprintf(stdout, "[%d] (%ld, %ld) loci: %ld, range: (%ld, %ld)\n",
               i, m->beg, m->end, m->range.top-m->range.bot+1, m->range.bot, m->range.top);
      }
      fprintf(stdout,"\n");
   }

   return mems;

}

wstack_t *
skip_seeds
(
const char    * seq,
const index_t   idx,
const size_t    gamma,
const size_t    skip
)
{
   int len = strlen(seq);
   int end = len-1;
   while (end > 0 && NONALPHABET[(uint8_t)seq[end]]) end--;

   if (end == 0)
      return stack_new(1);

   // Initialize mem stack
   wstack_t * seeds = stack_new(32);

   range_t range;

   // Iterate over all read positions
   while (end >= gamma - 1) {
      seed_t seed = {0};
      seed.end = end;

      // Backward <<<
      range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
      int mpos = end, mlen = 0;

      // Query the beginning of the read in lookup table.
      if (end >= LUTK - 1) {
         size_t merid = 0;
         for ( ; mlen < LUTK ; mlen++, mpos--) {
            if (NONALPHABET[(uint8_t) seq[end-mlen]]) {
               range.bot = 1;
               range.top = 0;
               break;
            }
            uint8_t c = ENCODE[(uint8_t) seq[end-mlen]];
            merid = c + (merid << 2);
         }
         range = idx.lut->kmer[merid];
      }

      // Cancel if we went too far already.
      if (range.top < range.bot) {
         range = (range_t) { .bot = 1, .top = idx.occ->txtlen-1 };
         mpos = end; mlen = 0;
      }

      for ( ; mlen < gamma ; mpos--, mlen++) {
         if (NONALPHABET[(uint8_t)seq[mpos]])
            break;

         int c = ENCODE[(uint8_t) seq[mpos]];
         range.bot = get_rank(idx.occ, c, range.bot - 1);
         range.top = get_rank(idx.occ, c, range.top) - 1;
         // Stop if no hits.
         if (range.top < range.bot)
            break;
      }

      if (mlen == gamma) {
         // Update seed info
         seed.beg = seed.end - gamma + 1;
         seed.range = range;
         // Push seed to stack
         seed_t * m = malloc(sizeof(seed_t));
         exit_on_memory_error(m);
         memcpy(m, &seed, sizeof(seed_t));
         push(m, &seeds);
      }

      // Update end position
      end -= (skip + 1);
      if (end < gamma - 1 && end > gamma - 1 - skip)
         end = gamma - 1;
   }

   return seeds;
}

aligncd_t
mem_alignments
(
 seedchain_t * chain,
 index_t       idx,
 size_t        slen
)
{
   // Allocate alignment candidates.
   align_t * alncd = malloc(chain->loci * sizeof(align_t));
   exit_on_memory_error(alncd);

   // Get all genomic positions.
   size_t nloc = 0;
   for (int j = 0; j < chain->pos; j++) {
      seed_t * seed = chain->seed[j];
      if (seed->aligned) {
         continue;
      }
      seed->aligned = 1;
      // Compute SA positions.
      if (!seed->sa)
         seed->sa = query_csa_range(idx.csa, idx.bwt, idx.occ, seed->range);
      // Make chained alignment candidates from seed genomic positions.
      for (int k = 0; k < seed->range.top - seed->range.bot + 1; k++) {
         alncd[nloc++] = (align_t){seed->sa[k], seed->end - seed->beg + 1, 0, seed};
      }
   }

   if (nloc == 0) {
      if (DEBUG_VERBOSE) {
         fprintf(stdout, "[skip: chain] No alncd after removing alignment duplicates\n");
      }
      free(alncd);
      return (aligncd_t){0, NULL};
   }

   if (DEBUG_VERBOSE)
      fprintf(stdout, "[MEM chain] %ld alncd after removing alignment duplicates\n", nloc);

   // Sort alncd by genomic position.
   qsort(alncd, nloc, sizeof(align_t), seed_by_refpos);

   // Chain alncd to avoid duplicated alignments.
   int cur = 0;
   size_t chain_gap_beg = alncd[cur].refpos + (alncd[cur].seed->end - alncd[cur].seed->beg + 1);
   size_t chain_gap_end = alncd[cur].refpos + (slen - alncd[cur].seed->beg) - 1;

   size_t nchain = 1;
   for (int j = 1; j < nloc; j++) {
      // Chain alncd if they are within 'slen' genomic nucleotides.
      size_t seed_ref_beg = alncd[j].refpos;
      size_t seed_ref_end = seed_ref_beg + (alncd[j].seed->end - alncd[j].seed->beg);
      if (alncd[cur].seed != alncd[j].seed && chain_gap_beg <= seed_ref_beg && seed_ref_end <= chain_gap_end) {
	 // Seed j is within current chain.
	 // Update seed chain span.
	 alncd[cur].span += (alncd[j].seed->end - alncd[j].seed->beg + 1);
	 // Set seed j span to 0.
	 alncd[j].span = 0;
      } else {
	 // Seed j is outside current chain.
	 // Update current chain.
	 cur = j;
	 // Update gap coordinates.
	 chain_gap_beg = alncd[cur].refpos + (alncd[cur].seed->end - alncd[cur].seed->beg + 1);
	 chain_gap_end = alncd[cur].refpos + (slen - alncd[cur].seed->beg) - 1;
	 // Update chain count.
	 nchain++;
      }
   }

   if (DEBUG_VERBOSE)
      fprintf(stdout, "[MEM chain] %ld alncd after chaining\n", nchain);

   // Sort alncd by span and align them.
   qsort(alncd, nloc, sizeof(align_t), seed_by_span);

   return (aligncd_t){nchain, alncd};
}

void
align
(
 align_t       alignment,
 const char  * seq,
 char        * genome,
 size_t        genome_len,
 int         * best_score,
 alnstack_t ** best
)
{
   seed_t * seed = alignment.seed;
   size_t   slen = strlen(seq);

   if (seed->beg > alignment.refpos)
      return;

   if (alignment.refpos + slen - seed->beg >= genome_len)
      return;

   int score = -1;
   if (seed->beg == 0 && seed->end == slen-1) {
      // Do not align perfect seeds (this should not happen because
      // these events are detected earlier).
      score = 0;
      if (DEBUG_VERBOSE)
         fprintf(stdout, "skip alignment: perfect seed -> score: 0\n");

   } else {
      // Get genomic sequence
      char * ref = decompress_genome(genome, alignment.refpos - seed->beg, slen+3);

      score = nw(ref,
            seq,
            slen+3, // Allow 3 nucleotides to allocate insertions.
            slen,
            *best_score + 1
      );

      // VERBOSE ALIGNMENT (DEBUG)
      if (DEBUG_VERBOSE) {
         fprintf(stdout, "locus: %ld (best_score: %d)\n", alignment.refpos, *best_score);
         fprintf(stdout, "%.*s %.*s %.*s\n",
               (int)seed->beg, seq, (int) (seed->end - seed->beg + 1),
               seq + seed->beg, (int)(slen-seed->end-1), seq + seed->end + 1);
         fprintf(stdout, "%.*s %.*s %.*s\nscore: %d\n--\n",
               (int)seed->beg, ref, (int) (seed->end - seed->beg + 1),
               ref + seed->beg, (int)(slen-seed->end-1), ref + seed->end + 1,
               score);
      }

      free(ref);
   }
   // Check align score.
   if (score <= *best_score) {
      // Create new alignment.
      aln_t aln = {0};
      aln.score    = score;
      aln.refpos   = alignment.refpos - seed->beg;
      aln.read_beg = seed->beg;
      aln.read_end = seed->end;
      aln.refseq   = NULL;

      if (score < *best_score) {
         // Reset best align stack.
         *best_score = score;
         (*best)->pos = 0;
      }
      aln_push(aln, best);
   }

   // Bug control on skip seeds
   // TODO: make sure with Eduard that it is OK to remove that bit.
//   if (score < alignment.minscore) {
//      fprintf(stderr, "bug found: score (%d) < expected minscore (%d), sequence: %s\n", score, alignment.minscore, seq);
//      exit(1);
//   }
}

seed_t *
filter_longest_mem
(
 wstack_t * seeds
)
{
  if (seeds->pos == 0) return 0;
  // Find longest MEMs in seeds.
  int n_mems = seeds->pos;
  seed_t * bestseed = seeds->ptr[0];
  int maxlen = bestseed->end - bestseed->beg + 1;
  for (int i = 1; i < n_mems; i++) {
    seed_t * s = seeds->ptr[i];
    int len = s->end - s->beg + 1;
    if (len > maxlen) {
      maxlen = len;
      free(bestseed);
      bestseed = s;
    }
    else {
      free(s);
    }
  }
  seeds->ptr[0] = bestseed;
  seeds->pos = 1;
  return bestseed;
}


alnstack_t *
remap_with_skip_seeds
(
          wstack_t   * skipseeds,
    const alnstack_t * mem_alst,
    const char       * seq,
    const index_t      idx,
    const int          max_mismatches,
    const size_t       circ_num
)
{

   size_t slen = strlen(seq);

   if (DEBUG_VERBOSE) {
      fprintf(stdout, "\n[REMAPPING]\n");
   }

   // Distinguish the MEM-based best score from the new best score.
   const int init_best_score = min(max_mismatches, MAX_ALIGN_MISMATCHES);
   int best_score = init_best_score;

   // Allocate new best hits.
   alnstack_t * best = alnstack_new(10);

   // Merge overlapping seeds (calls 'malloc' and 'exit_on_memory_error').
   wstack_t * merged = merge_overlapping_seeds(skipseeds, idx);

   // Sort one-locus seeds by loci in ascending order.
   qsort(merged->ptr, merged->pos, sizeof(seed_t *), seed_by_first_locus);

   // Chain seeds.

   // The leftmost seed is a position [1] on reads of 50 nt.
   const int n_slots = 1 + (slen-19) / 10;
   const int left_slot = slen-19 - 10 * (n_slots-1);

   // Start with leftmost seed (leftmost in the genome).
   seed_t * old = (seed_t *) merged->ptr[0];
   seed_t * new = old;
   seed_t * longest = old; // Longest seed of the chain.

   // Count minimum number of errors on the left of the seed.
   // The formula assumes that seeds are spaced by 10 nt.
   int minscore = (old->beg - left_slot + 10) / 20;
   // Calls 'malloc' and 'exit_on_memory_error'.
   wstack_t * chains = stack_new(16);
   for (int i = 1 ; i < merged->pos ; i++) {
      new = (seed_t *) merged->ptr[i];
      if (new->sa[0] - old->sa[0] < slen) {
         // New seed is in the same chain.
         // The formula assumes that seeds are spaced by 10*k + 2 nt.
         minscore += (new->beg - old->end + 18) / 20;
         if (new->end - new->beg >= longest->end - longest->beg) {
            // New seed is at least as long as old seed: replace.
            longest = new;
         }
      }
      else {
         // New seed is in another chain. Store previous chain.
         // Count minimum number of errors on the right of the seed.
         // The formula assumes that seeds are spaced by 10 nt.
         minscore += (slen-1 - old->end + 10) / 20;
         // Use field 'aligned' to store the minimum score of the seed.
         longest->aligned = minscore;
         push(longest, &chains);
         minscore = (new->beg - left_slot + 10) / 20;
         longest = new;
      }
      old = new; // Old is the new new.
   }
   // Store last chain.
   minscore += (slen-1 - new->end + 10) / 20;
   longest->aligned = minscore;
   push(longest, &chains);

   // Align seeds.
   int nseen = 0;
   for (size_t i = 0 ; i < chains->pos ; i++) {
      seed_t * seed = (seed_t *) chains->ptr[i];

      if (DEBUG_VERBOSE) {
         fprintf(stdout, "chain[%ld] (main seed) beg: %ld, end: %ld\n",
               i, seed->beg, seed->end);
      }

      // Check if locus can beat the best hit.
      if (seed->aligned >= init_best_score) {
         if (DEBUG_VERBOSE) {
            fprintf(stdout, "minscore %d >= previous best: skipping\n--\n",
                  seed->aligned);
         }
         continue;
      }

      // Check if locus was alreay aligned.
      size_t indexpos = (seed->sa[0] - seed->beg) / slen;
      int skip_alignment = 0;
      // Run through MEM-based alignments.
      for (int j = 0 ; j < 50 && skip_alignment == 0; j++) {
         if (mem_alst->seen[j] == 0)
            break;
         if (indexpos == mem_alst->seen[j])
            skip_alignment = 1;
      }
      if (skip_alignment) {
         if (DEBUG_VERBOSE) {
            fprintf(stdout, "locus %ld already aligned: skipping\n--\n",
                  seed->sa[0]-seed->beg);
         }
         continue;
      }
      else if (nseen < 50) {
         best->seen[nseen++] = indexpos; // Not used in present version.
         if (DEBUG_VERBOSE) {
            fprintf(stdout, "adding locus %ld to alignment list\n--\n",
                  seed->sa[0] - seed->beg);
         }
      }
      else if (DEBUG_VERBOSE) {
         fprintf(stdout, "skipping locus %ld (limit exceeded)\n--\n",
                  seed->sa[0] - seed->beg);
      }

      align_t alignment = (align_t) {
         .refpos = seed->sa[0],
         .span = seed->end - seed->beg + 1,
         .minscore = best_score,
         .seed = seed,
      };

      if (DEBUG_VERBOSE) {
         char buffer[256] = {0};
         fprintf(stdout, "alignment: refpos: %ld (%s), span: %ld\n",
               alignment.refpos,
               chr_string(alignment.refpos, idx.chr, buffer),
               alignment.span
         );
      }

      align(alignment, seq, idx.dna, idx.occ->txtlen, &best_score, &best);

   }

   // Copy genomic sequences for best alignments.
   for (size_t i = 0 ; i < best->pos ; i++) {
      best->aln[i].refseq = decompress_genome(
         idx.dna,
         best->aln[i].refpos,
         slen + best->aln[i].score
      );
   }

   // Free 'merged' and 'chains'.
   for (size_t i = 0 ; i < merged->pos ; i++) {
      seed_t * seed = (seed_t *) merged->ptr[i];
      free(seed->sa);
      free(seed);
   }
   free(merged);
   free(chains); // All seeds of 'chains' were in 'merged'.

   return best;

}



alnstack_t *
mapread
(
   wstack_t      * seeds,
   const char    * seq,
   const index_t   idx,
   const int       max_mismatches,
   const size_t    circ_num
)
{

   size_t slen = strlen(seq);

   if (DEBUG_VERBOSE) {
      fprintf(stdout, "\n[READ] sequence: %s\n",seq);
   }

   int best_score = min(max_mismatches, MAX_ALIGN_MISMATCHES);
   alnstack_t * best = alnstack_new(10);

   // See if we can compute the score without performing any alignment.
   seed_t * rightmost = (seed_t *) seeds->ptr[0];
   if (rightmost->end == slen-1) {
      if (rightmost->beg == 0) {
         if (DEBUG_VERBOSE) {
            fprintf(stdout, "MEM seed aligns perfectly: skipping alignment (score 0)\n");
         }
         rightmost->sa = query_csa_range(idx.csa, idx.bwt, idx.occ, rightmost->range);
         size_t nloci = rightmost->range.top - rightmost->range.bot + 1;
         for (size_t i = 0 ; i < nloci && i < 10; i++) {
            // No need to access the genome and pay a cache miss:
            // the sequence is the same as that of the read.
            char * refseq = strndup(seq, 1024);
            exit_on_memory_error(refseq);
            aln_t aln = (aln_t) {
               .score = 0,
               .refpos = rightmost->sa[i],
               .refseq = refseq,
               .read_beg = 0,
               .read_end = slen-1,
               .qual = 0. / 0.,
            };
            best->aln[i] = aln;
            best->seen[i] = rightmost->sa[i];
         }
         best->pos = nloci < 10 ? nloci : 10;
         return best;
      }
      else {
         // See if there exists a seed from [0] to [rightmost->beg-1].
         // NOTE: we miss the cases where the error is on the flanks of
         // the read, but this case is rare so it would complicate the code
         // with minor speed benefits.
         seed_t * leftmost = (seed_t *) seeds->ptr[seeds->pos-1];
         if (leftmost->beg == 0 && leftmost->end == rightmost->beg-2) {
            // Need to check if the seeds match contiguous regions.
            rightmost->sa = query_csa_range(idx.csa, idx.bwt, idx.occ, rightmost->range);
            leftmost->sa = query_csa_range(idx.csa, idx.bwt, idx.occ, leftmost->range);
            ssize_t dhit = rightmost->beg; // Target distance between loci.
            size_t nloci_right = rightmost->range.top - rightmost->range.bot + 1;
            size_t nloci_left = leftmost->range.top - leftmost->range.bot + 1;
            for (int i = 0 ; i < nloci_right && best->pos < 10 ; i++) {
            for (int j = 0 ; j < nloci_left && best->pos < 10 ; j++) {
               if (rightmost->sa[i] - leftmost->sa[j] == dhit) {
                  if (DEBUG_VERBOSE) {
                     fprintf(stdout, "Pair of MEM seeds aligns with one mismatch: skipping alignment (score 1)\n");
                  }
                  aln_t aln = (aln_t) {
                     .score = 1,
                     .refpos = leftmost->sa[j],
                     .refseq = decompress_genome(idx.dna, leftmost->sa[j], slen+1),
                     .read_beg = 0,
                     .read_end = slen-1,
                     .qual = 0. / 0.,
                  };
                  best->aln[best->pos] = aln;
                  best->seen[best->pos] = leftmost->sa[j];
                  best->pos++;
               }
            }
            }
         }
         if (best->pos > 0) {
            return best;
         }
      }
   }

   // Sort seeds by loci, ascending.
   qsort(seeds->ptr, seeds->pos, sizeof(seed_t *), mem_by_loci);

   // Alignable MEMS.
   size_t n_mem = 0;
   for (size_t i = 0; i < seeds->pos; i++, n_mem++) {
      seed_t * s = (seed_t *)seeds->ptr[i];
      if (s->range.top - s->range.bot >= MAX_MEM_SEED_LOCI)
         break;
   }

   // Compute repeats minscore.
   int repeats_minscore = slen+1;
   for (size_t i = n_mem; i < seeds->pos; i++) {
      seed_t * s = (seed_t *)seeds->ptr[i];
      repeats_minscore = min(repeats_minscore, (s->beg > 0) + (s->end < slen-1));
   }

   if (DEBUG_VERBOSE) {
      fprintf(stdout, "MEMs: %ld, alignable: %ld, repeat-minscore: %d\n",
            seeds->pos, n_mem, repeats_minscore);
   }

   int nseen = 0;

   // Sort alignable seeds by span.
   qsort(seeds->ptr, n_mem, sizeof(seed_t *), mem_by_span);

   // Align seeds from largest to smallest.
   for (size_t i = 0; i < n_mem; i++) {
      seed_t * mem = (seed_t *)seeds->ptr[i];
      int minscore = (mem->beg > 0) + (mem->end < slen-1);

      if (DEBUG_VERBOSE) {
         fprintf(stdout, "MEM[%ld] beg: %ld, end: %ld, span: %ld, minscore: %d, loci: %ld\n",
               i, mem->beg, mem->end, mem->end-mem->beg+1, minscore,
               mem->range.top - mem->range.bot + 1);
      }

      // Avoid alignments with score that is too high.
      if (minscore >= repeats_minscore || minscore > best_score ||
            (minscore == best_score && best->pos >= MAX_MINSCORE_REPEATS)) {
         continue;
      }

      // Get SA values (may have been done in rare cases alrady).
      if (mem->sa == NULL) mem->sa = query_csa_range(idx.csa, idx.bwt, idx.occ, mem->range);

      // Align at each locus.
      ssize_t nloci = mem->range.top - mem->range.bot + 1;
      for (ssize_t v = 0; v < nloci; v++) {
         ssize_t k = (v+circ_num)%nloci;
         // Check if locus was alreay aligned.
         size_t indexpos = (mem->sa[k] - mem->beg) / slen;
         int skip_alignment = 0;
         for (int j = 0 ; j < nseen && skip_alignment == 0; j++) {
            if (indexpos == best->seen[j]) {
               skip_alignment = 1;
            }
         }
         if (skip_alignment) {
            if (DEBUG_VERBOSE) {
               fprintf(stdout, "locus %ld already aligned: skipping\n--\n",
                     mem->sa[k]-mem->beg);
            }
            continue;
         }
         else if (nseen < 50) {
            best->seen[nseen++] = indexpos;
            if (DEBUG_VERBOSE) {
               fprintf(stdout, "adding locus %ld to alignment list\n--\n",
                     mem->sa[k]-mem->beg);
            }
         }
         else if (DEBUG_VERBOSE) {
            fprintf(stdout, "skipping locus %ld (limit alignments exceeded)\n--\n",
                     mem->sa[k]-mem->beg);
         }
         align_t alignment = (align_t){mem->sa[k], mem->end - mem->beg + 1, minscore, mem};
         if (DEBUG_VERBOSE) {
            char buffer[256] = {0};
            fprintf(stdout, "alignment: refpos: %ld (%s), span: %ld, MEM minscore: %d\n",
                  alignment.refpos,
                  chr_string(alignment.refpos, idx.chr, buffer),
                  alignment.span,
                  minscore
            );
         }

         align(alignment, seq, idx.dna, idx.occ->txtlen, &best_score, &best);

         // Stop alignments if necessary conditions met
         if (minscore == best_score && best->pos >= MAX_MINSCORE_REPEATS)
            break;
      }
   }

   // Check if best alignment is 100% unique
   if (best_score >= repeats_minscore) {
      best->pos = 0;
   }

   // Copy genomic sequences for best alignments.
   for (size_t i = 0; i < best->pos; i++)
      best->aln[i].refseq = decompress_genome(idx.dna, best->aln[i].refpos, slen + best->aln[i].score);

   // Free seeds
   //   for (size_t i = 0; i < seeds->pos; i++) {
   //      seed_t * s = (seed_t *) seeds->ptr[i];
   //      free(s->sa);
   //      free(s);
   //   }
   //   free(seeds);

   return best;

}


seedchain_t *
seedchain_new
(
 size_t max
 )
{
   size_t base = sizeof(seedchain_t);
   size_t extra = max * sizeof(seed_t *);
   seedchain_t * stack = malloc(base + extra);
   exit_on_memory_error(stack);

   stack->max = max;
   stack->pos = 0;
   stack->loci = 0;
   stack->minscore = 0;
   stack->span = 0;

   return stack;
}


void
seed_push
(
 seed_t       * mem,
 seedchain_t ** stackp
 )
{

   seedchain_t * stack = *stackp;

   if (stack->pos >= stack->max) {
      size_t newmax = stack->max*2;
      stack = *stackp = realloc(stack,
				sizeof(seedchain_t)+newmax*sizeof(seed_t *));
      exit_on_memory_error(stack);
      stack->max = newmax;
   }

   stack->seed[stack->pos++] = mem;
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
   bzero(stack->seen, 50 * sizeof(ssize_t));
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
