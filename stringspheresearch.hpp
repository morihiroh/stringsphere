#pragma once
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "stringdistance.hpp"

typedef unsigned long long int CountType;
typedef long double FloatType;

template <typename DistanceTag>
void estimate(
    const int k,  // the size of alphabet
    const std::vector<int>& center, // center string
    const int radius, // radius
    DistanceTag, // distance function
    const CountType iterations, // least # iterations 
    const CountType maxiter, // max # iterations 
    const int ell, // specify ell (-1: all)
    const bool is_only_u // estimate only u
    ) {

  if (log2((double)k)*(center.size() + radius) > sizeof(CountType)*8) {
    std::cerr << "The type will cause overflow error.\n";
    return;
  }

  CountType random_A, random_B;
  if (k == 2 or k == 4) {
    random_A = 5;
    random_B = 3;
  } else if (k == 3) {
    random_A = 4;
    random_B = 2;
  } else {
    std::cerr << "The case of k >= 5 will be implemented in the future.\n";
    return;
  }

  int s = center.size();

  if (radius == 0 or not is_only_u) {
    std::cout << k << "\t" << s << "\t0\t1";
    if (not is_only_u) {
      std::cout << "\t1";
    }
    std::cout << "\n";
  }

  if (radius >= 1) {
    std::vector<int> sn(center.size() + radius);
    std::vector<CountType> ur(radius + 1);

    CountType numiter = 0;
    ur[0] = 1;
    for (int r = is_only_u ? radius : 1; r <= radius; r++) {
      int c = s - r;
      if (c < 0) c = 0;
      int lmax = s + r;
      if (ell >= 0) {
        c = ell;
        lmax = ell;
      }
      CountType urlsum = 0;
      CountType kl = 1;
      for (int i = 0; i < c; i++) {
        kl *= k;
      }
      bool iterflag = true;
      for (int l = c; iterflag and l <= lmax; l++, kl *= k) {
        sn.resize(l);
        CountType least_iterations = std::min(iterations, kl);
        CountType x = 0;
        CountType rnd = std::time(nullptr) % kl;
        FloatType fkl = kl;
        FloatType coeff = 400.0 * fkl * (1-1/fkl);
        for (CountType n = 1; n <= kl; n++) {
          rnd = (rnd * random_A + random_B) % kl;
          std::lldiv_t tmp;
          tmp.quot = rnd;
          for (int i = 0; i < l; i++) {
            tmp = std::lldiv(tmp.quot, k);
            sn[i] = tmp.rem;
          }
          auto d = measure(center, sn, DistanceTag());
          if (within(d, r, typename DistanceTag::CompareType())) {
            x++;
          }
          numiter++;
          if (maxiter > 0 and numiter > maxiter) {
            iterflag = false;
            FloatType p = (FloatType)x/(FloatType)n;
            urlsum += (CountType)(0.5 + p * fkl);
            break;
          } else if (n >= least_iterations and x > 0) {
            FloatType fn = n;
            FloatType p = (FloatType)x/fn;
            if (fn >= coeff * p * (1-p) * (fkl - fn)) {
              //std::cout << "l=" << l << "\tx=" << x << "\tu=" << p*fkl << "\tp=" << p << "\t" << (FloatType)n/(FloatType)kl << "\tn=" << n << "\tkl=" << kl << "\n";
              urlsum += (CountType)(0.5 + p * fkl);
              break;
            }
          }
        }
      }
      std::cout << k << "\t" << s << "\t" << r << "\t" << urlsum;
      if (ell >= 0) {
        std::cout << "\t(ell=" << ell << ")";
      }
      if (not is_only_u) {
        ur[r] = urlsum;
        std::cout << "\t" << urlsum - ur[r-1];
      }
      std::cout << "\n";
    }
  }
}

template <typename DistanceTag>
void searchall(
    const int k,  // the size of alphabet
    const std::vector<int>& center, // center string
    const int radius, // radius
    DistanceTag // distance function
    ) {
  std::vector<int> sn(center.size() + radius);
  CountType xn = 0;
  CountType yn = 0;
  int lstart = center.size() - radius;
  int lend = center.size() + radius;
  if (lstart < 0) {
    lstart = 0;
  }
  for (int l = lstart; l <= lend; l++) {
    sn.resize(l);
    for (int i = 0; i < l; i++) {
      sn[i] = 0;
    }
    bool cf = true;
    while (cf) {
      int d = measure(center, sn, DistanceTag());
      if (within(d, radius, typename DistanceTag::CompareType())) {
        xn++;
        if (d == radius) {
          yn++;
        }
      }
      for (int i = 0; cf and i < l; i++) {
        sn[i]++;
        if (sn[i] == k) {
          cf = true;
          sn[i] = 0;
        } else {
          cf = false;
        }
      }
      cf = not cf;
    }
  }

  std::cout << k << "\t" << center.size() << "\t" << radius << "\t" << xn << "\t" << yn << "\n";
}

