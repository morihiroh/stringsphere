#pragma once
#include <vector>
#include <cmath>
#include <ctime>
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
    const int lomit, // omit smaller l
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

    ur[0] = 1;
    for (int r = is_only_u ? radius : 1; r <= radius; r++) {
      int c = s - r;
      if (c < 0) c = 0;
      if (lomit >= 0) {
        if (c < s + r - lomit) {
          c = s + r - lomit;
        }
      }
      CountType urlsum = 0;
      CountType kl = 1;
      for (int i = 0; i < c; i++) {
        kl *= k;
      }
      for (int l = c; l <= s + r; l++, kl *= k) {
        sn.resize(l);
        CountType least_iterations = std::min(iterations, kl);
        CountType x = 0;
        CountType rnd = std::time(nullptr) % kl;
        for (CountType n = 1; n <= kl; n++) {
          rnd = (rnd * random_A + random_B) % kl;
          CountType tmp = rnd;
          for (int i = 0; i < l; i++) {
            sn[i] = tmp % k;
            tmp /= k;
          }
          auto d = measure(center, sn, DistanceTag());
          if (within(d, r, typename DistanceTag::CompareType())) {
            x++;
          }
          FloatType fkl = kl;
          FloatType p = (FloatType)x/(FloatType)n;
          FloatType t = fkl - 1.0/(400.0*p*(1-p)*fkl/(fkl - 1) + 1.0/fkl);
          if (x > 0 and n >= least_iterations and n >= t) {
            //std::cout << "l=" << l << "\tx=" << x << "\tu=" << p*fkl << "\tp=" << p << "\t" << (FloatType)n/(FloatType)kl << "\tn=" << n << "\tkl=" << kl << "\n";
            urlsum += (CountType)(0.5 + p * fkl);
            break;
          }
        }
      }
      std::cout << k << "\t" << s << "\t" << r << "\t" << urlsum;
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

