#ifndef EC_H
#define EC_H

#include "gmp.h"
#include "math.h"
#include "pthread.h"
#include "stdio.h"

struct ec {
  mpz_t a;
  mpz_t b;
  mpz_t n;
};

struct point {
  int infty;
  int init;
  mpz_t x;
  mpz_t y;
};

void print_point(const struct point *p) {
  gmp_printf("x: %Zd, y: %Zd\n", p->x, p->y);
}

void free_point(struct point *p) {
  if (p->init) {
    mpz_clears(p->x, p->y, NULL);
    p->init = 0;
  }
}

void copy_point(struct point *ret, const struct point *p) {
  if (ret == p)
    return;

  ret->infty = p->infty;
  mpz_set(ret->x, p->x);
  mpz_set(ret->y, p->y);
  ret->init = 1;
}

int point_dbl(struct point *ret, const struct point *p,
              const struct ec *curve) {
  copy_point(ret, p);

  if (p->infty)
    return 0;

  mpz_t lambda;

  mpz_init(lambda);
  mpz_add(lambda, p->y, p->y); // lambda = 2y

  if (!mpz_invert(lambda, lambda, curve->n)) { // lambda = (2y)^-1
    mpz_gcd(lambda, lambda, curve->n);
    if (!mpz_cmp(lambda, curve->n))
      return -1;
    else {
      mpz_t q;
      mpz_init(q);
      mpz_divexact(q, curve->n, lambda);
      gmp_printf("SOLUTION: %Zd = %Zd * %Zd\n", curve->n, lambda, q);
      mpz_clears(q, lambda, NULL);
      return 1;
    }
    return 0;
  }

  mpz_t tmp_1, tmp_2;
  mpz_inits(tmp_1, tmp_2, NULL);

  mpz_mul(tmp_1, p->x, p->x);
  mpz_mul_ui(tmp_1, tmp_1, 3);
  mpz_mod(tmp_1, tmp_1, curve->n);
  mpz_add(tmp_1, tmp_1, curve->a); // tmp_1 = 3x^2 + a
  mpz_mul(lambda, tmp_1, lambda);  // lambda = (2y)^-1 * (3x^2 + a)
  mpz_add(tmp_1, p->x, p->x);      // tmp_1 = 2x
  mpz_set(tmp_2, p->x);
  mpz_mul(ret->x, lambda, lambda);
  mpz_sub(ret->x, ret->x, tmp_1); // ret->x = lambda^2 - 2x
  mpz_mod(ret->x, ret->x, curve->n);
  mpz_sub(tmp_1, tmp_2, ret->x); // tmp_1 = x - x'
  mpz_set(tmp_2, p->y);
  mpz_mul(ret->y, lambda, tmp_1); // ret->y = lambda(x - x')
  mpz_sub(ret->y, ret->y, tmp_2);
  mpz_mod(ret->y, ret->y, curve->n);

  mpz_clears(lambda, tmp_1, tmp_2, NULL);

  return 0;
}

int point_dbl_re(struct point *ret, const struct point *p,
                 const struct ec *curve, int pow) {
  copy_point(ret, p);
  for (int i = 0; i < pow; ++i) {
    int ret_code = point_dbl(ret, ret, curve);
    if (ret_code)
      return ret_code;
  }

  return 0;
}

int point_add(struct point *ret, const struct point *p, const struct point *q,
              const struct ec *curve) {
  if (p->infty) {
    copy_point(ret, q);
    return 0;
  } else if (q->infty) {
    copy_point(ret, p);
    return 0;
  }

  mpz_t tmp;
  mpz_init(tmp);
  mpz_neg(tmp, q->y);
  if (!mpz_cmp(p->y, tmp)) {
    ret->infty = 1;
    mpz_clear(tmp);
    return 0;
  }
  mpz_clear(tmp);
  ret->infty = 0;

  copy_point(ret, p);

  mpz_t lambda, tmp_1;

  mpz_inits(lambda, tmp_1, NULL);
  mpz_sub(lambda, q->y, p->y);
  mpz_sub(tmp_1, q->x, p->x);
  if (!mpz_invert(tmp_1, tmp_1, curve->n)) {
    mpz_gcd(tmp_1, tmp_1, curve->n);
    if (!mpz_cmp(tmp_1, curve->n))
      return -1;
    else {
      mpz_t r;
      mpz_init(r);
      mpz_divexact(r, curve->n, tmp_1);
      gmp_printf("SOLUTION: %Zd = %Zd * %Zd\n", curve->n, tmp_1, r);
      mpz_clears(r, lambda, tmp_1, NULL);
      return 1;
    }
  }

  mpz_t tmp_2;
  mpz_init(tmp_2);

  mpz_mul(lambda, lambda, tmp_1); // lambda = (x_q - x_p)^-1 * (y_q - y_p)
  mpz_mod(lambda, lambda, curve->n);
  mpz_mul(tmp_2, lambda, lambda);
  mpz_sub(tmp_2, tmp_2, p->x);
  mpz_sub(tmp_2, tmp_2, q->x); // tmp_2 = lambda^2 - x_p - x_q
  mpz_mod(tmp_2, tmp_2, curve->n);
  mpz_sub(tmp_1, p->x, tmp_2);
  mpz_mul(tmp_1, lambda, tmp_1);
  mpz_sub(ret->y, tmp_1, p->y); // ret->y = lambda(x_p - x_r) - y_p
  mpz_mod(ret->y, ret->y, curve->n);
  mpz_set(ret->x, tmp_2);

  mpz_clears(lambda, tmp_1, tmp_2, NULL);

  return 0;
}

int point_mul(struct point *ret, const struct point *p, const struct ec *curve,
              unsigned long long k) {
  struct point q = {1, 0};
  mpz_inits(q.x, q.y, NULL);
  q.init = 1;

  struct point tmp = {0, 0};
  copy_point(&tmp, p);

  while (1) {
    if (k & 1) {
      int ret_code = point_add(&q, &q, &tmp, curve);
      if (ret_code)
        return ret_code;
    }

    if (!k)
      break; // save one point doubling

    int ret_code = point_dbl(&tmp, &tmp, curve);
    if (ret_code)
      return ret_code;
    k >>= 1;
  }

  free_point(&tmp);

  copy_point(ret, &q);

  return 0;
}

int point_window(struct point *ret, const struct point *p,
                 const struct ec *curve, unsigned long long k) {
  if (k < (1 << 7))
    return point_mul(ret, p, curve, k);

  struct point window[8];

  struct point pre_comp = {0, 0};
  mpz_inits(pre_comp.x, pre_comp.y, NULL);
  pre_comp.init = 1;
  int ret_code = point_dbl_re(&pre_comp, p, curve, 3);
  if (ret_code)
    return ret_code;

  for (int i = 0; i < 8; ++i) {
    mpz_inits(window[i].x, window[i].y, NULL);
    window[i].init = 1;
  }

  for (int i = 0; i < 8; ++i) {
    copy_point(window + i, &pre_comp);
    ret_code = point_add(&pre_comp, &pre_comp, p, curve);
    if (ret_code)
      return ret_code;
  }

  free_point(&pre_comp);
  struct point q = {1, 0};
  mpz_inits(q.x, q.y, NULL);
  q.init = 1;

  unsigned long long k_copy = k;
  int m = 0;
  while (k_copy >>= 1)
    m += 1;

  for (int i = m; i >= 0; --i) {
    int bit = (k >> i) & 1;
    if (!bit) {
      ret_code = point_dbl(&q, &q, curve);
      if (ret_code)
        return ret_code;
    } else {
      if (i < 3) {
        int t = k & 0b111;
        struct point tmp = {0, 0};
        mpz_inits(tmp.x, tmp.y, NULL);
        tmp.init = 1;
        ret_code = point_mul(&tmp, p, curve, t);
        if (ret_code) {
          return ret_code;
        }
        ret_code = point_add(&q, &q, &tmp, curve);
        if (ret_code) {
          return ret_code;
        }
      } else {
        int bits = (k >> (i - 3)) & 0b1111;
        ret_code = point_dbl_re(&q, &q, curve, 4);
        if (ret_code)
          return ret_code;
        ret_code = point_add(&q, &q, window + (bits - 8), curve);
        if (ret_code)
          return ret_code;
      }

      i -= 3;
    }
  }

  for (int i = 0; i < 8; ++i)
    free_point(window + i);

  ret->init = 0;
  copy_point(ret, &q);

  free_point(&q);

  return 0;
}

#endif