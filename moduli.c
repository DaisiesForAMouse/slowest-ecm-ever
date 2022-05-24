#include "factor.h"

int main(int argc, char **argv) {
  done = 0;
  int n_threads;
  int verbose;
  switch (argc) {
  case 1:
    printf("enter number to factor on cmd line call\n");
    return -1;
    break;
  case 2:
    n_threads = 2;
    verbose = 0;
    break;
  case 3:
    n_threads = strtol(*(argv + 2), NULL, 10);
    verbose = 0;
    break;
  case 4:
    n_threads = strtol(*(argv + 2), NULL, 10);
    verbose = **(argv + 3) == 'v';
    break;
  default:
    printf("too many arguments");
    break;
  }

  return get_factor(*(argv + 1), n_threads, verbose);
}