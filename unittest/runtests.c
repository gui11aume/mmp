#include "unittest.h"

char open_banner[] = "";
char close_banner[] = "---\n";

int
main(
   int     argc,
   char ** argv
)
{

   // Import test cases.
   extern test_case_t test_cases_bwt[];

   // Register test cases.
   const test_case_t *list_of_test_cases[] = {
      test_cases_bwt,
      NULL,
   };

   // Display banner.
   fprintf(stderr, "%s", open_banner);

   // Run the tests.
   int nfailed = run_unittest(argc, argv, list_of_test_cases);

   // Display banner.
   fprintf(stderr, "%s", close_banner);

   return nfailed;

}
