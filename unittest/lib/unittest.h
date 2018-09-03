//  # unittests_your_module.c
//  void test_your_function (void) { int x = 1; test_assert(x); }
//  const test_case_t test_cases_your_module[] = {
//     {"your_module/your_function",  test_your_function},
//     {NULL, NULL},
//  };

//  # runtests.c
//  extern test_case_t test_cases_your_module[];
//  const test_case_t *list_of_test_cases[] = {
//     test_cases_your_module,
//     NULL,
//  };
//  int nfailed = run_unittest(argc, argv, list_of_test_cases);

#include <getopt.h>
#include <fcntl.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#ifndef _UNITTEST_HEADER
#define _UNITTEST_HEADER

//     Type definitions      // 
struct test_case_t;

typedef struct test_case_t test_case_t;
typedef void (*fixture_t)(void);

struct test_case_t {
   char      * test_name;
   fixture_t   fixture;
};


//     Function declarations     //

// Main function to run tests //
int run_unittest (int, char **, const test_case_t **);

// Interception of 'malloc()' and 'realloc()' //
void * malloc (size_t size);
void * realloc (void *ptr, size_t size);
void * calloc (size_t nitems, size_t size);
void   set_alloc_failure_rate_to (double);
void   set_alloc_failure_countdown_to (int);
void   reset_alloc (void);

// Interception of stderr //
char * caught_in_stderr (void);
void   redirect_stderr (void);
void   unredirect_stderr (void);


// Functions that must be defined here but that should not be used //
void fail_critical (const char *, const char *, int, const char *);
void fail_non_critical (const char *, const char *, int, const char *);
void fail_debug_dump (const char *, int, const char *);


//    Assert macros    //
#define test_assert(expr)  do { \
  if (expr)  { (void)0; } \
  else { \
     fail_debug_dump(__FILE__, __LINE__, __func__); \
     fail_non_critical(#expr, __FILE__, __LINE__, __func__); \
}} while (0)

#define test_assert_critical(expr) do { \
  if (expr)  { (void)0; } \
  else { \
     fail_debug_dump(__FILE__, __LINE__, __func__); \
     fail_critical(#expr, __FILE__, __LINE__, __func__); \
     return; \
}} while (0)

#define test_assert_stderr(expr) do { \
  if (strncmp(expr, caught_in_stderr(), strlen(expr)) == 0) { (void)0; } \
  else { \
     fail_debug_dump(__FILE__, __LINE__, __func__); \
     fail_non_critical(#expr, __FILE__, __LINE__, __func__); \
}} while (0)

#endif
