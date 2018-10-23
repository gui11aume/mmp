#include "map.h"

typedef struct memchain_t memchain_t;
typedef struct seed_t     seed_t;

#define LEN 50
#define MEM_MAX_LOCI 1000
#define MAX_MINSCORE_REPEATS 2

#ifdef DEBUG
int DEBUG_VERBOSE = 1;
#else
int DEBUG_VERBOSE = 0;
#endif

#define min(x,y) ((x) < (y) ? (x) : (y))
#define min3(x,y,z) (min(min(x,y),z))

struct memchain_t {
   size_t    pos;
   size_t    max;
   size_t    loci;
   int       minscore;
   int       span;
   mem_t   * mem[];
};

struct seed_t {
   size_t   refpos;
   size_t   span;
   int      minscore;
   mem_t  * mem;
};

memchain_t *  memchain_new (size_t max);
void          mem_push     (mem_t * mem, memchain_t ** stackp);
alnstack_t *  alnstack_new (size_t max);
void          aln_push     (aln_t aln, alnstack_t ** stackp);

int           seed_by_refpos (const void * a, const void * b) {
   return ((seed_t *)a)->refpos > ((seed_t *)b)->refpos;
};
int           seed_by_minscore_then_span (const void * a, const void * b) {
   seed_t * sa = (seed_t *)a;
   seed_t * sb = (seed_t *)b;
   if (sa->minscore > sb->minscore)
      return 1;
   else if (sa->minscore < sb->minscore)
      return -1;
   else
      return sa->span < sb->span;
};

int           mem_by_start (const void * a, const void * b) {
   return (*(mem_t **)a)->beg > (*(mem_t **)b)->beg;
};
int           minscore_then_span (const void * a, const void * b) {
   int sa = (*(memchain_t **)a)->minscore;
   int sb = (*(memchain_t **)b)->minscore;
   if (sa > sb)
      return 1;
   else if (sa < sb)
      return -1;
   else
      return (*(memchain_t **)a)->span < (*(memchain_t **)b)->span;
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
recursive_mem_chain
(
 stack_t  * mems,
 size_t     mem_pos,
 size_t     chain_pos,
 mem_t   ** chain,
 stack_t ** chain_stack
 )
{
   size_t mem_end = ((mem_t *) mems->ptr[mem_pos])->end;
   // For MEMs overlapping MEM[pos].
   for (size_t i = mem_pos; i < mems->pos && ((mem_t *)mems->ptr[i])->beg <= mem_end; i++) {
      // Extend chain.
      chain[chain_pos] = (mem_t *) mems->ptr[i];
      // Get next nonoverlapping MEM.
      size_t j;
      for (j = i+1; j < mems->pos; j++) {
	 if (((mem_t *)mems->ptr[j])->beg > ((mem_t *) mems->ptr[i])->end)
	    break;
      }
      // We reached the chain end, store chain.
      if (j >= mems->pos) {
	 memchain_t * memchain = memchain_new(chain_pos+1);
	 // Push mems to chain and compute span.
	 int span = 0;
	 size_t loci = 0;
	 for (size_t k = 0; k <= chain_pos; k++) {
	    span += chain[k]->end - chain[k]->beg + 1;
	    loci += chain[k]->range.top - chain[k]->range.bot + 1;
	    mem_push(chain[k], &memchain);
	 }
	 // Push mem chain to chain stack.
	 memchain->span = span;
	 memchain->loci = loci;
	 push(memchain, chain_stack);
      } else {
	 recursive_mem_chain(mems, j, chain_pos+1, chain, chain_stack);
      }
   }
}


stack_t *
nonoverlapping_mems
(
 stack_t * mems
 )
{
   // Alloc.
   stack_t * chain_stack = stack_new(8);
   mem_t  ** chain = malloc(mems->pos*sizeof(mem_t *));
   exit_on_memory_error(chain);

   // 1. Sort MEMs by start position.
   qsort(mems->ptr, mems->pos, sizeof(mem_t *), mem_by_start);

   // 2. Recursive call to mem group.
   recursive_mem_chain(mems, 0, 0, chain, &chain_stack);
   free(chain);
   // 3. Return stack of non-overlapping MEM combinations.
   return chain_stack;
}


int
chain_min_score
(
 memchain_t * chain,
 const int    gamma,
 const int    seqlen
 )
{
   int minscore = 0;

   // Add mismatches at chain ends.
   
   if (chain->mem[0]->beg > 0)
      minscore += ((chain->mem[0]->beg - 1) / gamma) + 1;

   if (chain->mem[chain->pos-1]->end < seqlen-1)
      minscore += (seqlen - 2 - chain->mem[chain->pos-1]->end)/gamma + 1;

   // Add gap mismatches.
   for (int i = 1; i < chain->pos; i++) {
      // A gap implies one mismatch, even if it's a gap of length 0.
      minscore++;
      int gap_size = chain->mem[i]->beg - chain->mem[i-1]->end - 1;
      minscore += gap_size > 1;
      minscore += (gap_size-2)/gamma;
   }

   return minscore;
}


alnstack_t *
align
(
 const index_t   idx,
 const char    * seq,
 const char    * genome,
 stack_t       * mems,
 const int       gamma
 )
{
   int slen = strlen(seq);

   // DEBUG VERBOSE
   if (DEBUG_VERBOSE) {
      fprintf(stdout,"\nMEMs (%ld):\n", mems->pos);
      for(int i = 0; i < mems->pos; i++) {
	 mem_t * m = (mem_t *) mems->ptr[i];
	 fprintf(stdout, "[%d] (%ld, %ld) range: (%ld, %ld)\n", i, m->beg, m->end, m->range.bot, m->range.top);
      }
      fprintf(stdout,"\n");
   }

   // Filter repeats > MEM_MAX_LOCI.
   for (size_t i = 0; i < mems->pos; i++) {
      mem_t * m = (mem_t *) mems->ptr[i];
      size_t  mem_loci = m->range.top - m->range.bot + 1;
      if (mem_loci > MEM_MAX_LOCI) {
	 m->range.top = m->range.bot + MEM_MAX_LOCI - 1;
	 if (DEBUG_VERBOSE)
	    fprintf(stdout, "[MEM %ld] too many loci (%ld > MEM_MAX_LOCI(%d))\n", i, mem_loci, MEM_MAX_LOCI);
      }
   }


   // Find all non-overlapping MEM combinations.
   stack_t * chain_stack = nonoverlapping_mems(mems);

   // Compute minimum alignment score given seed distribution.
   for (int i = 0; i < chain_stack->pos; i++) {
      memchain_t * chain = (memchain_t *)chain_stack->ptr[i];
      if (chain->pos)
	 chain->minscore = chain_min_score(chain, gamma, slen);
      else
	 fprintf(stderr, "OJOOOOOOOOOOOOOOOOOOOOOOOOOO AQUI UN CHAIN SENSE MEMS!");
   }

   // Sort mem chains by minscore(inc) then span(dec).
   qsort(chain_stack->ptr, chain_stack->pos, sizeof(memchain_t *), minscore_then_span);

   // DEBUG VERBOSE
   if (DEBUG_VERBOSE) {
      fprintf(stdout,"MEM chains (%ld):\n", chain_stack->pos);
      for(int i = 0; i < chain_stack->pos; i++) {
	 memchain_t * c = (memchain_t *) chain_stack->ptr[i];
	 fprintf(stdout, "Chain [%d] (mems: %ld, span: %d, minscore: %d):\n", i, c->pos, c->span, c->minscore);
	 for (int j = 0; j < 0; j++) {
	    mem_t * m = c->mem[j];
	    fprintf(stdout, "[%d] (%ld, %ld) range: (%ld, %ld)\n", i, m->beg, m->end, m->range.bot, m->range.top);
	 }
      }
      fprintf(stdout,"\n");
   }


   int best_score = slen;
   alnstack_t * best = alnstack_new(10);
   
   // Iterate over sorted MEM chains.
   for (int i = 0; i < chain_stack->pos; i++) {
      memchain_t * chain = (memchain_t *) chain_stack->ptr[i];

      if (DEBUG_VERBOSE) {
	 fprintf(stdout, "[MEM chain %d] mems: %ld, span: %d, loci: %ld, minscore: %d\n", i, chain->pos, chain->span, chain->loci, chain->minscore);
      }

      if (chain->minscore > best_score) {
	 if (DEBUG_VERBOSE)
	    fprintf(stdout, "[break: chain] chain.minscore > best_score\n");
	 break;
      }
      
      if (chain->minscore == best_score && best->pos >= MAX_MINSCORE_REPEATS) {
	 if (DEBUG_VERBOSE)
	    fprintf(stdout, "[break: chain] %d+ alignments with best_score(%d)\n", MAX_MINSCORE_REPEATS, best_score);
	 break;
      }

      // Allocate seeds.
      seed_t * seeds = malloc(chain->loci * sizeof(seed_t));
      exit_on_memory_error(seeds);

      // Get all genomic positions.
      size_t nloc = 0;
      for (int j = 0; j < chain->pos; j++) {
	 mem_t * mem = chain->mem[j];
	 if (mem->aligned) {
	    continue;
	 }
	 mem->aligned = 1;
	 // Compute SA positions.
	 if (!mem->sa)
	    mem->sa = query_csa_range(idx.csa, idx.bwt, idx.occ, mem->range);
	 // Make chained seeds from MEM genomic positions.
	 for (int k = 0; k < mem->range.top - mem->range.bot + 1; k++) {
	    seeds[nloc++] = (seed_t){mem->sa[k], mem->end - mem->beg + 1, 0, mem};
	 }
      }

      if (nloc == 0 && DEBUG_VERBOSE) {
	 fprintf(stdout, "[skip: chain] No seeds after removing duplicates\n");
	 continue;
      }

      if (DEBUG_VERBOSE)
	 fprintf(stdout, "[MEM chain %d] %ld seeds after removing duplicates\n", i, nloc);

      // Sort seeds by genomic position.
      qsort(seeds, nloc, sizeof(seed_t), seed_by_refpos);

      // Chain seeds to avoid duplicated alignments.
      int cur = 0;
      memchain_t * seedchain = memchain_new(8);
      mem_push(seeds[0].mem, &seedchain);

      size_t chain_gap_beg = seeds[cur].refpos + (seeds[cur].mem->end - seeds[cur].mem->beg + 1);
      size_t chain_gap_end = seeds[cur].refpos + (slen - seeds[cur].mem->beg) - 1;

      size_t nchain = 1;
      for (int j = 1; j < nloc; j++) {
	 // Chain seeds if they are within 'slen' genomic nucleotides.
	 size_t seed_ref_beg = seeds[j].refpos;
	 size_t seed_ref_end = seed_ref_beg + (seeds[j].mem->end - seeds[j].mem->beg);
	 if (seeds[cur].mem != seeds[j].mem && chain_gap_beg <= seed_ref_beg && seed_ref_end <= chain_gap_end) {
	    seeds[cur].span += (seeds[j].mem->end - seeds[j].mem->beg + 1);
	    seeds[j].span = 0;
	    seeds[j].minscore = slen;
	    mem_push(seeds[j].mem, &seedchain);
	 } else {
	    // Compute minimum score of seed chain.
	    seeds[cur].minscore = chain_min_score(seedchain, gamma, slen);
	    // Update current seed.
	    cur = j;
	    // Reset seed chain.
	    seedchain->pos = 0;
	    mem_push(seeds[cur].mem, &seedchain);
	    // Update gap coordinates.
	    chain_gap_beg = seeds[cur].refpos + (seeds[cur].mem->end - seeds[cur].mem->beg + 1);
	    chain_gap_end = seeds[cur].refpos + (slen - seeds[cur].mem->beg) - 1;
	    // Update chain count.
	    nchain++;
	 }
      }

      if (DEBUG_VERBOSE)
	 fprintf(stdout, "[MEM chain %d] %ld seeds after chaining\n", i, nchain);
	    
      // Compute minimum score of last seed chain.
      seeds[cur].minscore = chain_min_score(seedchain, gamma, slen);
	    
      free(seedchain);

      // Sort seeds by span and align them.
      qsort(seeds, nloc, sizeof(seed_t), seed_by_minscore_then_span);

      for (int j = 0; j < nloc; j++) {
	 seed_t seed = seeds[j];
	 // Do not align chained seeds.
	 if (!seed.span) break;
	 if (DEBUG_VERBOSE) {
	    char * strpos  = chr_string(seed.refpos, idx.chr);
	    fprintf(stdout, "seed: refpos: %ld (%s), span: %ld, minscore: %d\n", seed.refpos, strpos, seed.span, seed.minscore);
	    free(strpos);
	 }
	 // Stop when minscore > best.
	 if (seed.minscore > best_score) {
	    if (DEBUG_VERBOSE)
	       fprintf(stdout, "[break: align] seed.minscore > best_score\n");
	    break;
	 }
	 if (seed.minscore == best_score && best->pos >= MAX_MINSCORE_REPEATS) {
	    if (DEBUG_VERBOSE)
	       fprintf(stdout, "[break: align] %d+ alignments with best_score(%d)\n", MAX_MINSCORE_REPEATS, best_score);
	    break;
	 }
	 
	 mem_t * mem = seed.mem;
	 const char * ref = genome + seed.refpos - mem->beg;
	 
	 int score;
	 // Do not align perfect seeds.
	 if (mem->beg == 0 && mem->end == slen-1)	{
	    score = 0;
	    if (DEBUG_VERBOSE) 
	       fprintf(stdout, "skip alignment: perfect seed -> score: 0\n");

	 } else {
	    score = nw(ref,
		       seq,
		       slen+3, // Allow 3 nucleotides to allocate insertions.
		       slen,
		       best_score + 1
		       );
	 
	    // VERBOSE ALIGNMENT (DEBUG)
	    if (DEBUG_VERBOSE) {
	       fprintf(stdout, "locus: %ld (seed.minscore: %d, best_score: %d)\n", seed.refpos, seed.minscore, best_score);
	       fprintf(stdout, "%.*s %.*s %.*s\n",
		       (int)mem->beg, seq, (int) (mem->end - mem->beg + 1),
		       seq + mem->beg, (int)(slen-mem->end-1), seq + mem->end + 1);
	       fprintf(stdout, "%.*s %.*s %.*s\nscore: %d\n--\n",
		       (int)mem->beg, ref, (int) (mem->end - mem->beg + 1),
		       ref + mem->beg, (int)(slen-mem->end-1), ref + mem->end + 1,
		       score);
	    }
	 }
         // Check align score.
         if (score <= best_score) {
            // Create new alignment.
            aln_t aln;
            aln.score  = score;
            aln.nmem   = 1;
            aln.refpos = seed.refpos - mem->beg;
            aln.refseq = ref;
            aln.mem    = *mem;

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
      }
      free(seeds);
   }
      
   // When we have the seed(s) that gave the best hit, we need to
   // evaluate the number of copies of that sequence, and the
   // divergence. We said we would do this with direct computation of
   // the log-likelihood like: set mu, compute N and get loglik.
   // For this, we need to keep more dense record of the seeding
   // process for every seed.

   for (int i = 0; i < chain_stack->pos; i++) {
      free(chain_stack->ptr[i]);
   }
   free(chain_stack);
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
   stack_t * mems = stack_new(50);

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
      if (mlen >= gamma) {
	 mem_t * m = malloc(sizeof(mem_t));
	 exit_on_memory_error(m);
	 memcpy(m, &mem, sizeof(mem_t));
	 push(m, &mems);
      }

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

   // Return the best alignment(s) in an alignment stack.
   alnstack_t * aln = align(idx, seq, genome, mems, gamma);

   for (size_t i = 0; i < mems->pos; i++)
      free(mems->ptr[i]);

   free(mems);
   return aln;

}


memchain_t *
memchain_new
(
 size_t max
 )
{
   size_t base = sizeof(memchain_t);
   size_t extra = max * sizeof(mem_t *);
   memchain_t * stack = malloc(base + extra);
   exit_on_memory_error(stack);

   stack->max = max;
   stack->pos = 0;
   stack->loci = 0;
   stack->minscore = 0;
   stack->span = 0;

   return stack;
}


void
mem_push
(
 mem_t       * mem,
 memchain_t ** stackp
 )
{

   memchain_t * stack = *stackp;

   if (stack->pos >= stack->max) {
      size_t newmax = stack->max*2;
      stack = *stackp = realloc(stack,
				sizeof(memchain_t)+newmax*sizeof(mem_t *));
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
