#include "unittest.h"

//     Constants     //
#define UTEST_BUFFER_SIZE 1024
#define MAX_N_ERROR_MSG 20

// Private functions // 
void   redirect_stderr (void);
void * run_test (void *); 
char * sent_to_stderr (void);
void   unit_test_clean (void);
void   unit_test_init (void);
void   unredirect_stderr (void);
void   test_case_clean (void);
void   test_case_init (void);

//     Global variables     //
static double ALLOC_ERR_PROB;
static int    ALLOC_FAIL_COUNTER;
static FILE * DEBUG_DUMP_FILE;
static int    N_ERROR_MSG;
static int    ORIG_STDERR_DESCRIPTOR;
static char * STDERR_BUFFER;
static int    STDERR_OFF;
static int    TEST_CASE_FAILED;

static int    SHOWALL = 0;


//     Function definitions     //

void
update_display_failed
(void)
{
   fprintf(stderr, "[FAIL]\n");
}

void
update_display_success
(void)
{
   fprintf(stderr, "  [OK]\n");
}

void
terminate_thread
(
   int sig
)
{

   // If test was not failed so far, update display.
   if (!TEST_CASE_FAILED) {
      update_display_failed();
   }

   fprintf(stderr, "caught SIGTERM (interrupting)\n");

   // Label the test case as failed. //
   TEST_CASE_FAILED = 1;

   // Clean it all. //
   unredirect_stderr();
   test_case_clean();

   // Return to main thread.
   pthread_exit(NULL);

}


void *
run_test
(
   void * data
)
{

   // Unpack argument. //
   const test_case_t test_case = *(test_case_t *) data;

   // -- Initialize -- //

   // Register signal handlers. This will allow execution
   // to return to main thread in case of crash.
   signal(SIGSEGV, terminate_thread);
   
   // Get test name and fixture (test function). //
   fixture_t fixture = test_case.fixture;

   // Initialize variable, empty buffers... //
   test_case_init();
   unredirect_stderr();

   //  -- Run the test -- //
   (*fixture)();
   
   // -- Clean -- //
   unredirect_stderr();
   test_case_clean();

   return NULL;

}


void
unit_test_init
(void)
{

   DEBUG_DUMP_FILE = fopen(".inspect.gdb", "w");
   if (DEBUG_DUMP_FILE == NULL) {
      fprintf(stderr, "unittest error: %s:%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

}

void
unit_test_clean
(void)
{

   fprintf(DEBUG_DUMP_FILE, "b run_unittest\n");
   fprintf(DEBUG_DUMP_FILE, "run\n");
   fclose(DEBUG_DUMP_FILE);

}


void
test_case_init
(void) 
{

   reset_alloc();

   N_ERROR_MSG = 0;

   TEST_CASE_FAILED = 0;
   STDERR_OFF = 0;

   ORIG_STDERR_DESCRIPTOR = dup(STDERR_FILENO);
   if (ORIG_STDERR_DESCRIPTOR == -1) {
      fprintf(stderr, "unittest error (%s:%d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

   STDERR_BUFFER = calloc(UTEST_BUFFER_SIZE, sizeof(char));
   if (STDERR_BUFFER == NULL) {
      fprintf(stderr, "unittest error (%s:%d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

}


void
test_case_clean
(void)
{

   free(STDERR_BUFFER);
   STDERR_BUFFER = NULL;

}
   

void
redirect_stderr
(void)
{
   if (STDERR_OFF) return;
   // Flush stderr, redirect to /dev/null and set buffer.
   fflush(stderr);
   int temp = open("/dev/null", O_WRONLY);
   if (temp == -1) {
      fprintf(stderr, "unittest error: %s:%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   if (dup2(temp, STDERR_FILENO) == -1) {
      fprintf(stderr, "unittest error: %s:%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   memset(STDERR_BUFFER, '\0', UTEST_BUFFER_SIZE * sizeof(char));
   if (setvbuf(stderr, STDERR_BUFFER, _IOFBF, UTEST_BUFFER_SIZE) != 0) {
      fprintf(stderr, "unittest error: %s:%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   if (close(temp) == -1) {
      fprintf(stderr, "unittest error: %s:%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   fprintf(stderr, "$");
   fflush(stderr);
   STDERR_OFF = 1;
}


void
unredirect_stderr
(void)
{
   if (!STDERR_OFF) return;
   fflush(stderr);
   if (dup2(ORIG_STDERR_DESCRIPTOR, STDERR_FILENO) == -1) {
      // Could not restore stderr. No need to give an error
      // message because it will not appear.
      exit(EXIT_FAILURE);
   }
   if (setvbuf(stderr, NULL, _IONBF, 0) != 0) {
      fprintf(stderr, "unittest error: %s:%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   STDERR_OFF = 0;
}


void
fail_non_critical
(
   const char * assertion,
   const char * file,
         int    lineno,
   const char * function
)
// Handle failure of non critical assert statements.
// Set the test status to failed and update the display
// and print error message.
{

   // If first failed assertion of the test
   // case, update display for the user.
   if (!TEST_CASE_FAILED) {
      update_display_failed();
   }

   TEST_CASE_FAILED = 1;

   // If stderr is redirected, we will need to
   // take it back to display the error message.
   int toggle_stderr = STDERR_OFF;
   if (toggle_stderr) unredirect_stderr();

   // Don't show more than 'MAX_N_ERROR_MSG', unless
   // user passed the --showall or -a option.
   if (N_ERROR_MSG++ < MAX_N_ERROR_MSG || SHOWALL) {
      fprintf(stderr, "%s:%d: `%s'\n", file, lineno, assertion);
   }
   else if (N_ERROR_MSG == MAX_N_ERROR_MSG + 1) {
      fprintf(stderr, "more than %d failed assertions...\n",
            MAX_N_ERROR_MSG);
   }

   // Flush stderr and put it back as it was (ie redirect
   // it if it was redirected, or leave it as is otherwise).
   fflush(stderr);
   if (toggle_stderr) redirect_stderr();

}


void
fail_critical
(
   const char * assertion,
   const char * file,
         int    lineno,
   const char * function
)
// Handle failure of critical (fatal) assert statements.
// Set the test status to failed and update the display
// and print error message.
{
   // If test was not failed so far, update display.
   if (!TEST_CASE_FAILED) {
      update_display_failed();
   }

   TEST_CASE_FAILED = 1;

   // That's the end of the test case. We can
   // take back stderr without worries.
   unredirect_stderr();
   fprintf(stderr, "assertion failed in %s, %s:%d: `%s' (CRITICAL)\n",
         function, file, lineno, assertion);
   fflush(stderr);

}


void
fail_debug_dump
(
   const char * file,
         int    lineno,
   const char * function
)
// When a test fails, write a gdb file with a breakpoint at
// the line where the failured happened.
{
   fprintf(DEBUG_DUMP_FILE, "b %s:%d\n", file, lineno);
   fflush(DEBUG_DUMP_FILE);
}


char *
caught_in_stderr
(void)
{
   return STDERR_BUFFER;
}



#ifdef __APPLE__
// Intercept calls to 'malloc()' on MacOS //
#include <malloc/malloc.h>
void *(*SYSTEM_MALLOC) (malloc_zone_t *, size_t);
void *(*SYSTEM_CALLOC) (malloc_zone_t *, size_t, size_t);
void *(*SYSTEM_REALLOC) (malloc_zone_t *, void *, size_t);

void
mac_specific_initialization
(void)
{
   malloc_zone_t *zone = malloc_default_zone();
   SYSTEM_MALLOC = zone->malloc;
   SYSTEM_CALLOC = zone->calloc;
   SYSTEM_REALLOC = zone->realloc;
}

void *
fail_countdown_malloc(malloc_zone_t *zone, size_t size)
{
   if (ALLOC_FAIL_COUNTER >= 0) ALLOC_FAIL_COUNTER--;
   return ALLOC_FAIL_COUNTER < 0 ? NULL : SYSTEM_MALLOC(zone, size);
}

void *
fail_countdown_calloc(malloc_zone_t *zone, size_t nitems, size_t size)
{
   if (ALLOC_FAIL_COUNTER >= 0) ALLOC_FAIL_COUNTER--;
   return ALLOC_FAIL_COUNTER < 0 ? NULL : \
      SYSTEM_CALLOC(zone, nitems, size);
}

void *
fail_countdown_realloc(malloc_zone_t *zone, void *ptr, size_t size)
{
   if (ALLOC_FAIL_COUNTER >= 0) ALLOC_FAIL_COUNTER--;
   return ALLOC_FAIL_COUNTER < 0 ? NULL : \
      SYSTEM_REALLOC(zone, ptr, size);
}

void *
fail_prone_malloc(malloc_zone_t *zone, size_t size)
{
   return drand48() < ALLOC_ERR_PROB ? NULL : SYSTEM_MALLOC(zone, size);
}

void *
fail_prone_calloc(malloc_zone_t *zone, size_t nitems, size_t size)
{
   return drand48() < ALLOC_ERR_PROB ? NULL : \
      SYSTEM_CALLOC(zone, nitems, size);
}

void *
fail_prone_realloc(malloc_zone_t *zone, void *ptr, size_t size)
{
   return drand48() < ALLOC_ERR_PROB ? NULL : \
      SYSTEM_REALLOC(zone, ptr, size);
}


void
set_alloc_failure_rate_to(double p)
{

   ALLOC_ERR_PROB = p;
   malloc_zone_t *zone = malloc_default_zone();
   zone->malloc = fail_prone_malloc;
   zone->calloc = fail_prone_calloc;
   zone->realloc = fail_prone_realloc;
}

void
set_alloc_failure_countdown_to(int count)
{
   ALLOC_FAIL_COUNTER = count;
   malloc_zone_t *zone = malloc_default_zone();
   zone->malloc = fail_countdown_malloc;
   zone->calloc = fail_countdown_calloc;
   zone->realloc = fail_countdown_realloc;
}

void
reset_alloc
(void)
{
   ALLOC_ERR_PROB = 0.0;
   malloc_zone_t *zone = malloc_default_zone();
   zone->malloc = SYSTEM_MALLOC;
   zone->calloc = SYSTEM_CALLOC;
   zone->realloc = SYSTEM_REALLOC;
}
#else
// Intercept 'malloc()' calls in a Linux environment //
extern void *__libc_malloc(size_t);
extern void *__libc_calloc(size_t, size_t);
extern void *__libc_realloc(void *, size_t);

void *(*MALLOC_CALL)(size_t) = __libc_malloc;
void *(*REALLOC_CALL)(void *, size_t) = __libc_realloc;
void *(*CALLOC_CALL)(size_t, size_t) = __libc_calloc;

// Overwrite the symbols //
void * malloc (size_t size) {
   return (*MALLOC_CALL)(size);
}
void * realloc (void *ptr, size_t size) {
   return (*REALLOC_CALL)(ptr, size);
}
void * calloc (size_t nitems, size_t size) {
   return (*CALLOC_CALL)(nitems, size);
}

void *
fail_countdown_malloc(size_t size)
{
   if (ALLOC_FAIL_COUNTER >= 0) ALLOC_FAIL_COUNTER--;
   return ALLOC_FAIL_COUNTER < 0 ? NULL : __libc_malloc(size);
}

void *
fail_countdown_calloc(size_t nitems, size_t size)
{
   if (ALLOC_FAIL_COUNTER >= 0) ALLOC_FAIL_COUNTER--;
   return ALLOC_FAIL_COUNTER < 0 ? NULL : __libc_calloc(nitems, size);
}

void *
fail_countdown_realloc(void *ptr, size_t size)
{
   if (ALLOC_FAIL_COUNTER >= 0) ALLOC_FAIL_COUNTER--;
   return ALLOC_FAIL_COUNTER < 0 ? NULL : __libc_realloc(ptr, size);
}

void *
fail_prone_malloc(size_t size)
{
   return drand48() < ALLOC_ERR_PROB ? NULL : __libc_malloc(size);
}

void *
fail_prone_calloc(size_t nitems, size_t size)
{
   return drand48() < ALLOC_ERR_PROB ? NULL : __libc_calloc(nitems, size);
}

void *
fail_prone_realloc(void *ptr, size_t size)
{
   return drand48() < ALLOC_ERR_PROB ? NULL : __libc_realloc(ptr, size);
}


void
set_alloc_failure_rate_to(double p)
{
   ALLOC_ERR_PROB = p;
   MALLOC_CALL = fail_prone_malloc;
   CALLOC_CALL = fail_prone_calloc;
   REALLOC_CALL = fail_prone_realloc;
}

void
set_alloc_failure_countdown_to(int count)
{
   ALLOC_FAIL_COUNTER = count;
   MALLOC_CALL = fail_countdown_malloc;
   CALLOC_CALL = fail_countdown_calloc;
   REALLOC_CALL = fail_countdown_realloc;
}

void
reset_alloc
(void)
{
   ALLOC_ERR_PROB = 0.0;
   MALLOC_CALL = __libc_malloc;
   CALLOC_CALL = __libc_calloc;
   REALLOC_CALL = __libc_realloc;
}
#endif





int
run_unittest
(
         int            argc,
         char        ** argv,
   const test_case_t  * test_case_list[]
)
{

#ifdef __APPLE__
   mac_specific_initialization();
#endif

	// Parse command line options.
   while(1) {
      int option_index = 0;
      static struct option long_options[] = {
         {"showall", no_argument, &SHOWALL,  1},
         {0, 0, 0, 0}
      };

      int c = getopt_long(argc, argv, "a",
            long_options, &option_index);

      // Done parsing named options. //
      if (c == -1) break;

      switch (c) {
      case 0:
         break;
		case 'a':
			SHOWALL = 1;
			break;
		default:
			fprintf(stderr, "cannot parse command line arguments\n");
			return EXIT_FAILURE;
		}
	}

   // Initialize test. //
   unit_test_init();

   pthread_t tid;

   int nbad = 0;

   // Run test cases in sequential order.
   for (int j = 0 ; test_case_list[j] != NULL ; j++) {

      const test_case_t *tests = test_case_list[j];
      for (int i = 0 ; tests[i].fixture != NULL ; i++) {

         // Test display. //
         fprintf(stderr, "%s%*s", tests[i].test_name,
                  35 - (int) strlen(tests[i].test_name), "");

         // Run test case in thread. //
         pthread_create(&tid, NULL, run_test, (void *) &tests[i]);
         pthread_join(tid, NULL);

         if (TEST_CASE_FAILED) {
				// Display was updated if test failed. 
            nbad++;
         }
         else {
            update_display_success();
         }

      }

   }

   // Clean and return. //
   unit_test_clean();
   return nbad;

}
