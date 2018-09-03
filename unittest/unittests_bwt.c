#include "unittest.h"
#include "bwt.c"

void
test_compute_sa
(void)
{

   const char txt[] = "GATGCGAGAGATG";

   int64_t *SA = compute_sa(txt);
   test_assert_critical(SA != NULL);

   const int array[] = {13,6,8,10,1,4,12,5,7,9,0,3,11,2};
   for (int i = 0 ; i < 14 ; i++) {
      test_assert(SA[i] == array[i]);
   }

   free(SA);

}


void
test_compress_sa
(void)
{

   // Trivial text of size 257 (including terminator).
   const char txt[] =
      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT";

   int64_t *SA = compute_sa(txt);
   test_assert_critical(SA != NULL);

   csa_t *csa = compress_sa(SA);
   test_assert_critical(csa != NULL);

   test_assert(csa->nbits == 9);
   test_assert(csa->nint64 == 3);
   test_assert(csa->bmask == 0b111111111);

   uint64_t w0 =
     0b1001011111001001111000111111000101111000011111000001111100000000;
//     x|      95|      79|      63|      47|      31|      15|     256
   uint64_t w1 =
     0b1101100111101011111101010111101001111101000111100111111100110111;
//     xx|     207|     191|     175|     159|     143|     127|    111
   uint64_t w2 =
     0b0000000000000000000000000000000000000000111111110111011110110111;
//                                            |     255|     239|   223

   test_assert(csa->bitf[0] == w0);
   test_assert(csa->bitf[1] == w1);
   test_assert(csa->bitf[2] == w2);

   free(csa);
   free(SA);

}


void
test_create_bwt
(void)
{

   const char txt[] = "GATGCGAGAGATG";

   int64_t *SA = compute_sa(txt);
   test_assert_critical(SA != NULL);

   bwt_t *BWT = create_bwt(txt, SA);
   test_assert_critical(BWT != NULL);

   test_assert(BWT->txtlen == 14);
   test_assert(BWT->nslots == 4);
   test_assert(BWT->zero == 10);

   // The BWT is GGGG GGTC AA$T AA.
   test_assert(BWT->slots[0] == 0b10101010);
   test_assert(BWT->slots[1] == 0b01111010);
   test_assert(BWT->slots[2] == 0b11000000);
   test_assert(BWT->slots[3] == 0b00000000);

   free(BWT);
   free(SA);

}


void
test_write_occ_blocks
(void)
{

   size_t sz = sizeof(occ_t) + SIGMA * sizeof(blocc_t);
   occ_t *occ = calloc(1, sz);
   test_assert_critical(occ != NULL);

   occ->nrows = 1; // Important.

   uint32_t smpl[4] = {2,3,4,5};
   uint32_t bits[4] = {6,7,8,9};

   write_occ_blocks(occ, smpl, bits, 0);

   test_assert(occ->rows[0].smpl == 2);
   test_assert(occ->rows[1].smpl == 3);
   test_assert(occ->rows[2].smpl == 4);
   test_assert(occ->rows[3].smpl == 5);

   test_assert(occ->rows[0].bits == 6);
   test_assert(occ->rows[1].bits == 7);
   test_assert(occ->rows[2].bits == 8);
   test_assert(occ->rows[3].bits == 9);

   free(occ);

}


void
test_create_occ
(void)
{

   const char txt[] = "GATGCGAGAGATG";

   int64_t *SA = compute_sa(txt);
   test_assert_critical(SA != NULL);

   bwt_t *BWT = create_bwt(txt, SA);
   test_assert_critical(BWT != NULL);

   free(SA);

   occ_t *occ = create_occ(BWT);
   test_assert_critical(occ != NULL);

   test_assert(occ->txtlen == 14);
   test_assert(occ->nrows == 1);

   test_assert(occ->rows[0].smpl == 0);
   test_assert(occ->rows[1].smpl == 0);
   test_assert(occ->rows[2].smpl == 0);
   test_assert(occ->rows[3].smpl == 0);

   // The BWT is GGGGGGTCAA$TAA.
   test_assert(occ->rows[0].bits == 0b00000000110011000000000000000000);
   test_assert(occ->rows[1].bits == 0b00000001000000000000000000000000);
   test_assert(occ->rows[2].bits == 0b11111100000000000000000000000000);
   test_assert(occ->rows[3].bits == 0b00000010000100000000000000000000);

   test_assert(occ->C[0] == 1);
   test_assert(occ->C[1] == 5);
   test_assert(occ->C[2] == 6);
   test_assert(occ->C[3] == 12);
   test_assert(occ->C[4] == 14);

   free(occ);
   free(BWT);

}


void
test_get_rank
(void)
{

   const char txt[] = "GATGCGAGAGATG";

   int64_t *SA = compute_sa(txt);
   test_assert_critical(SA != NULL);

   bwt_t *BWT = create_bwt(txt, SA);
   test_assert_critical(BWT != NULL);

   free(SA);

   occ_t *occ = create_occ(BWT);
   test_assert_critical(occ != NULL);

   const int rankA[] = {1,1,1,1,1,1,1,1,2,3,3,3,4,5};
   const int rankC[] = {5,5,5,5,5,5,5,6,6,6,6,6,6,6};
   const int rankG[] = {7,8,9,10,11,12,12,12,12,12,12,12,12,12};
   const int rankT[] = {12,12,12,12,12,12,13,13,13,13,13,14,14,14};

   const int *rank[] = {rankA, rankC, rankG, rankT};

   for (int j = 0 ; j < 4 ; j++) {
   for (int i = 0 ; i < 14 ; i++) {
      test_assert(get_rank(occ, j, i) == rank[j][i]);
   }
   }

   free(occ);
   free(BWT);

}


void
test_fill_lut
(void)
{

   const char txt[] = "GATGCGAGAGATG";

   int64_t *SA = compute_sa(txt);
   test_assert_critical(SA != NULL);

   bwt_t *BWT = create_bwt(txt, SA);
   test_assert_critical(BWT != NULL);

   free(SA);

   occ_t *occ = create_occ(BWT);
   test_assert_critical(occ != NULL);

   lut_t *lut = calloc(1, sizeof(lut_t));
   test_assert_critical(lut != NULL);

   // k-mer "ATGCGAGAGAT" has ID 3285612
   fill_lut(lut, occ, (range_t) {.bot=4, .top=4} , 11, 3285612);
   // k-mer "GATGCGAGAGAT" has ID 13142450
   test_assert(lut->kmer[13142450].bot == 10);
   test_assert(lut->kmer[13142450].top == 10);

   free(lut);
   free(occ);
   free(BWT);

}


void
test_backward_search
(void)
{

   const char txt[] = "GATGCGAGAGATG";

   int64_t *SA = compute_sa(txt);
   test_assert_critical(SA != NULL);

   bwt_t *BWT = create_bwt(txt, SA);
   test_assert_critical(BWT != NULL);

   free(SA);

   occ_t *occ = create_occ(BWT);
   test_assert_critical(occ != NULL);

   range_t range = backward_search("GATGCGAGAGAT", 12, occ);

   test_assert(range.bot == 10);
   test_assert(range.top == 10);

   free(occ);
   free(BWT);

}

void
test_query_csa
(void)
{

   const char txt[] = "GATGCGAGAGATG";

   int64_t *SA = compute_sa(txt);
   test_assert_critical(SA != NULL);

   bwt_t *BWT = create_bwt(txt, SA);
   test_assert_critical(BWT != NULL);

   csa_t *csa = compress_sa(SA);
   test_assert_critical(csa != NULL);
   free(SA);

   occ_t *occ = create_occ(BWT);
   test_assert_critical(occ != NULL);

   const int array[] = {13,6,8,10,1,4,12,5,7,9,0,3,11,2};

   for (int i = 0 ; i < 14 ; i++) {
      test_assert(query_csa(csa, BWT, occ, i) == array[i]);
   }

   free(csa);
   free(occ);
   free(BWT);

}

// Test cases for export.
const test_case_t test_cases_bwt[] = {
   {"compute_sa",         test_compute_sa},
   {"compress_sa",        test_compress_sa},
   {"create_bwt",         test_create_bwt},
   {"write_occ_blocks",   test_write_occ_blocks},
   {"create_occ",         test_create_occ},
   {"get_rank",           test_get_rank},
   {"fill_lut",           test_fill_lut},
   {"backward_search",    test_backward_search},
   {"query_csa",          test_query_csa},
   {NULL, NULL},
};
