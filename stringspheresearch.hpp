#pragma once
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "stringdistance.hpp"

typedef unsigned long long int CountType;
typedef long double FloatType;

inline void radix_convert(std::vector<int>& s, const int k, const int x) {
  std::lldiv_t tmp;
  tmp.quot = x;
  int len = s.size();
  for (int i = 0; i < len; i++) {
    tmp = std::lldiv(tmp.quot, k);
    s[i] = tmp.rem;
  }
}

template <typename DistanceTag>
void estimate(
    const int method,
    const double gamma,
    const double epsilon,
    const int k,  // the size of alphabet
    const std::vector<int>& center, // center string
    const int radius, // radius
    DistanceTag, // distance function
    const CountType iterations, // least # iterations 
    const CountType maxiter // max # iterations 
                            //const int ell, // specify ell (-1: all)
                            //const bool is_only_u // estimate only u
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
  FloatType gammasq = gamma * gamma;

  std::cout << k << "\t" << s << "\t0\t1";
  std::cout << "\t1";
  std::cout << "\n";

  if (radius >= 1) {
    if (s == 0) { // Case 2
      CountType kl = k;
      CountType u = 1;
      for (int r = 1; r <= radius; r++, kl *= k) {
        u += kl;
        std::cout << k << "\t" << s << "\t" << r << "\t" << u << "\t" << kl << "\n";
      }
      return;
    }

    std::vector<int> sn(center.size() + radius);
    std::vector<CountType> ur(radius + 1);

    CountType numiter = 0;
    ur[0] = 1;  // for r = 0
    for (int r = 1; r <= radius; r++) {
      const int lmax = s + r;
      CountType urlsum = 0;
      CountType kl = 1;
      bool case3flag = false;
      int c = s - r;
      if (c <= 0) { // Case 3
        case3flag = true;
        c = r + 1;
        for (int i = 0; i < c; i++) {
          urlsum += kl;
          kl *= k;
        }
      } else {  // Case 4
        for (int i = 0; i < c; i++) {
          kl *= k;
        }
      }
      bool iterflag = true;
      for (int l = c; iterflag and l <= lmax; l++, kl *= k) { // Component
        sn.resize(l);
        CountType least_iterations = std::min(iterations, kl);
        CountType x = 0;
        CountType rnd = std::time(nullptr) % kl;  // initialization for random numbers
        FloatType fkl = kl;
        for (CountType n = 1; n <= kl; n++) {
          rnd = (rnd * random_A + random_B) % kl; // linear congruential generator
          radix_convert(sn, k, rnd);              // convert a number 'rnd' to a string

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
            FloatType pfkl = p * fkl;
            FloatType gsq = 0.0;
            if (method == 2) {  // absolute error
              if (case3flag) {
                if (epsilon >= s) {
                  double t = epsilon / s - 0.5;
                  gsq = t * t / s;
                } else {
                  gsq = 0.25 / s;
                }
              } else {
                double fr = 2.0 * r + 1.0;
                if (epsilon >= fr) {
                  double t = epsilon / fr - 0.5;
                  gsq = t * t / fr;
                } else {
                  gsq = 0.25 / fr;
                }
              }
            } else {            // relative error
              if (case3flag) {
                if (epsilon >= s / pfkl) {
                  FloatType t = epsilon * pfkl / s - 0.5;
                  gsq = t * t / s;
                } else {
                  gsq = 0.25 / s;
                }
              } else {
                double fr = 2.0 * r + 1.0;
                if (epsilon >= fr / pfkl) {
                  FloatType t = epsilon * pfkl / fr - 0.5;
                  gsq = t * t / fr;
                } else {
                  gsq = 0.25 / fr;
                }
              }
            }
            if (fn >= fkl / (1.0 + gsq * (fkl - 1.0) / (gammasq * pfkl * (1.0 - p) * fkl))) {
              urlsum += (CountType)(0.5 + pfkl);
              break;
            }
          }
        }
      }
      std::cout << k << "\t" << s << "\t" << r << "\t" << urlsum;
      ur[r] = urlsum;
      std::cout << "\t" << urlsum - ur[r-1];
      if (not iterflag) {
        std::cout << "\t(Warning: Forcibly finished because of too many trials (>=" << maxiter << "))";
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

