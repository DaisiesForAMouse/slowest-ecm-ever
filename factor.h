#ifndef FACTOR_H
#define FACTOR_H

#include "ec.h"
#include "gmp.h"
#include "primes.h"
#include "pthread.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "unistd.h"

#define POLLARD_ITER_MAX 10000000
#define ECM_CURVE_MAX 25
#define ECM_INIT_ITER 25000

struct ecm_info {
  int done;
  int verbose;
  gmp_randstate_t rstate;
  mpz_t n;
};

clock_t tic, toc;

int trial_division(const mpz_t n) {
  for (int i = 0; i < NUM_PRIMES; ++i) {
    if (mpz_divisible_ui_p(n, PRIMES[i])) {
      mpz_t q;
      mpz_init(q);
      mpz_divexact_ui(q, n, PRIMES[i]);
      gmp_printf("SOLUTION: %Zd = %u * %Zd\n", n, PRIMES[i], q);
      mpz_clear(q);
      return 0;
    }
  }

  return -1;
}

int rho(const mpz_t n, gmp_randstate_t rstate, int verbose) {
  mpz_t tort, hare, z, temp, p, q;
  mpz_inits(tort, hare, NULL);
  mpz_urandomm(tort, rstate, n);
  mpz_urandomm(hare, rstate, n);
  unsigned long power = 1, lam = 1;
  mpz_init_set_si(z, 1);
  mpz_init(temp);
  mpz_init(q);

  mpz_init_set_si(p, 1);

  unsigned long counter = 0;

  unsigned long iter_count = 0;

  while (!mpz_cmp_si(p, 1)) {
    if (power == lam) {
      mpz_set(tort, hare);
      power <<= 1;
      lam = 0;
    }

    mpz_mul(hare, hare, hare);
    mpz_add_ui(hare, hare, 1);
    mpz_mod(hare, hare, n);
    lam += 1;

    mpz_sub(temp, tort, hare);
    mpz_mul(z, temp, z);

    if (++counter == 100) {
      mpz_gcd(p, n, z);
      mpz_set_ui(z, 1);
      counter = 0;
    }

    if (!(counter % 5))
      mpz_mod(z, z, n);

    if (iter_count++ == POLLARD_ITER_MAX) {
      if (verbose)
        printf("hit maximum pollard iterations, moving to ecm\n");
      return -1;
    }
  }

  mpz_divexact(q, n, p);
  gmp_printf("SOLUTION: %Zd = %Zd * %Zd\n", n, p, q);

  mpz_clears(tort, hare, z, temp, p, q, NULL);

  return 0;
}

int ecm(const mpz_t n, long int iter_max, int curve_count,
        gmp_randstate_t rstate, int verbose) {
  struct ec curve;
  mpz_init_set(curve.n, n);
  mpz_inits(curve.a, curve.b, NULL);

  struct point p = {0, 0};
  mpz_inits(p.x, p.y, NULL);
  mpz_urandomm(p.x, rstate, n);
  mpz_urandomm(p.y, rstate, n);
  mpz_urandomm(curve.a, rstate, n);
  p.init = 1;

  mpz_t tmp;
  mpz_init(tmp);
  mpz_powm_ui(curve.b, p.y, 2, n);
  mpz_powm_ui(tmp, p.x, 3, n);
  mpz_sub(curve.b, curve.b, tmp); // b = y^2 - x^3 mod n
  mpz_mul(tmp, curve.a, p.x);
  mpz_sub(curve.b, curve.b, tmp);
  mpz_mod(curve.b, curve.b, n);

  int bound = 0;
  unsigned long long prod = 1;
  while (++bound < iter_max) {
    if (prod < 1 << 23)
      prod *= bound;
    else {
      int ret_code = point_window(&p, &p, &curve, prod);
      if (ret_code == 1)
        return 0;
      prod = 1;
    }
  }

  if (curve_count < ECM_CURVE_MAX) {
    return ecm(n, iter_max, curve_count + 1, rstate, verbose);
  } else {
    int new_max = iter_max << 1;
    if (verbose)
      printf("thread id %ld | increasing b to %d\n", (long)pthread_self(),
             new_max);
    return ecm(n, new_max, 0, rstate, verbose);
  }

  return -1;
}

void *ecm_default(void *th) {
  struct ecm_info *info = (struct ecm_info *)th;
  ecm(info->n, ECM_INIT_ITER, 0, info->rstate, info->verbose);
  if (info->verbose)
    printf("thread id %ld | done!\n", (long)pthread_self());
  info->done = 1;
  return NULL;
}

int get_factor(char* str, int n_threads, int verbose) {
  mpz_t n;
  mpz_init_set_str(n, str, 10);
  gmp_randstate_t rstate_main;
  gmp_randinit_default(rstate_main);
  gmp_randseed_ui(rstate_main, time(NULL));

  tic = time(NULL);
  if (mpz_probab_prime_p(n, 50) > 0) {
    gmp_printf("SOLUTION: %Zd = %Zd * 1\n", n, n);
  } else if (mpz_perfect_square_p(n)) {
    mpz_t p;
    mpz_init(p);
    mpz_sqrt(p, n);
    gmp_printf("SOLUTION: %Zd = %Zd * %Zd\n", n, p, p);
  } else {
    int trial_fail = trial_division(n);
    if (trial_fail) {
      int rho_fail = rho(n, rstate_main, verbose);
      if (rho_fail) {
        pthread_t threads[n_threads];
        struct ecm_info infos[n_threads];
        int done = 0;
        for (int i = 0; i < n_threads; ++i) {
          infos[i].done = 0;
          infos[i].verbose = verbose;
          mpz_init_set(infos[i].n, n);
          gmp_randinit_default(infos[i].rstate);
          gmp_randseed_ui(infos[i].rstate, time(NULL));
          for (int j = 0; j <= i; ++j)
            done += infos[j].done;
          if (done)
            break;
          pthread_create(threads + i, NULL, ecm_default, (void *)(infos + i));
          if (verbose)
            printf("thread %i created\n", i);
          usleep(1250 * 1000);
        }

        while (!done) {
          usleep(10 * 1000);
          for (int i = 0; i < n_threads; ++i)
            done += infos[i].done;
        }

        for (int i = 0; i < n_threads; ++i) {
          if (!infos[i].done)
            pthread_cancel(threads[i]);
        }
      }
    }
  }

  toc = time(NULL);

  if (verbose)
    printf("total time spent: %lus\n", toc - tic);

  return 0;
}

#endif