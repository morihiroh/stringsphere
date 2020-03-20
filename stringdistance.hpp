#pragma once
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <boost/math/distributions/normal.hpp>
#include <boost/multi_array.hpp>

#define DISABLE_JAROWINKLER

struct CompareDistance {};
struct CompareSimilarity {};
struct ExtendedHammingDistance {
  typedef CompareDistance CompareType;
};
struct LongestCommonSubsequence {
  typedef CompareSimilarity CompareType;
};
struct LevenshteinDistance {
  typedef CompareDistance CompareType;
};
struct DamerauLevenshteinDistance {
  typedef CompareDistance CompareType;
};

int measure(const std::vector<int>& x, const std::vector<int>& y, ExtendedHammingDistance) {
  int xsize = x.size();
  int ysize = y.size();
  int minsize = ysize;
  int cost = xsize - ysize;
  if (cost < 0) {
    minsize = xsize;
    cost = -cost;
  }
  for (int i = 0; i < minsize; i++) {
    if (x[i] != y[i]) {
      cost++;
    }
  }
  return cost;
}

int measure(const std::vector<int>& x, const std::vector<int>& y, LongestCommonSubsequence) {
  int xsize = x.size();
  int ysize = y.size();
  boost::multi_array<int, 2> d(boost::extents[xsize + 1][ysize + 1]);
  for (int i = 0; i <= xsize; i++) {
    d[i][0] = 0;
  }
  for (int j = 1; j <= ysize; j++) {
    d[0][j] = 0;
  }
  for (int i = 1; i <= xsize; i++) {
    for (int j = 1; j <= ysize; j++) {
      if (x[i - 1] == y[j - 1]) {
        d[i][j] = d[i - 1][j - 1] + 1;
      } else {
        d[i][j] = std::max(d[i - 1][j], d[i][j - 1]);
      }
    }
  }
  return d[xsize][ysize];
}

int measure(const std::vector<int>& x, const std::vector<int>& y, LevenshteinDistance) {
  const int InsertCost = 1;
  const int DeleteCost = 1;
  const int SubstCost = 1;
  int xsize = x.size();
  int ysize = y.size();
  boost::multi_array<int, 2> d(boost::extents[xsize + 1][ysize + 1]);

  for (int i = 0; i <= xsize; i++) {
    d[i][0] = i;
  }
  for (int j = 1; j <= ysize; j++) {
    d[0][j] = j;
  }
  for (int i = 1; i <= xsize; i++) {
    for (int j = 1; j <= ysize; j++) {
      int cost = (x[i - 1] == y[j - 1]) ? 0 : SubstCost;
      d[i][j] = std::min(std::min(
            d[i - 1][j] + InsertCost,
            d[i][j - 1] + DeleteCost),
          d[i - 1][j - 1] + cost);
    }
  }
  return d[xsize][ysize];
}

int measure(const std::vector<int>& x, const std::vector<int>& y, DamerauLevenshteinDistance) {
  const int InsertCost = 1;
  const int DeleteCost = 1;
  const int SubstCost = 1;
  int xsize = x.size();
  int ysize = y.size();
  boost::multi_array<int, 2> d(boost::extents[xsize + 1][ysize + 1]);

  for (int i = 0; i <= xsize; i++) {
    d[i][0] = i;
  }
  for (int j = 1; j <= ysize; j++) {
    d[0][j] = j;
  }
  for (int i = 1; i <= xsize; i++) {
    for (int j = 1; j <= ysize; j++) {
      int cost = (x[i - 1] == y[j - 1]) ? 0 : SubstCost;
      d[i][j] = std::min(std::min(
            d[i - 1][j] + InsertCost,
            d[i][j - 1] + DeleteCost),
          d[i - 1][j - 1] + cost);
      if (i > 1 and j > 1 and x[i - 1] == y[j - 2] and x[i - 2] == y[j - 1]) {
        d[i][j] = std::min(d[i][j], d[i - 2][j - 2] + cost);
      }
    }
  }
  return d[xsize][ysize];
}

bool within(const int d, const int radius, CompareDistance) {
  return d <= radius;
}

bool within(const int d, const int radius, CompareSimilarity) {
  return d >= radius;
}

#ifndef DISABLE_JAROWINKLER
double operator()(const std::vector<int>& x, const std::vector<int>& y, JaroWinklerDistance) {
  int xsize = x.size();
  int ysize = y.size();
  if (xsize == 0) {
    if (ysize == 0)
      return 1.0;
    else
      return 0.0;
  }
  int match_dist = std::max(xsize, ysize)/2 - 1;
  std::vector<bool> x_matches(xsize, false), y_matches(ysize, false);
  double m = 0.0;
  double t = 0.0;
  for (int i = 0; i < xsize; i++) {
    const int start = std::max(0, i - match_dist);
    const int end = std::min(i + match_dist + 1, ysize);
    for (int j = start; j < end; j++) {
      if (y_matches[j]) continue;
      if (x[i] != y[j]) continue;
      x_matches[i] = true;
      y_matches[j] = true;
      m++;
      break;
    }
  }
  if (m == 0) 
    return 0.0;

  int j = 0;
  for (int i = 0; i < xsize; i++) {
    if (not x_matches[i]) continue;
    while (not y_matches[j]) {
      j++;
    }
    if (x[i] != y[j]) {
      t++;
    }
    j++;
  }
  t /= 2.0;
  return (m/xsize + m/ysize + (m - t)/m)/3.0;
}
#endif

typedef unsigned long long int CountType;
typedef double FloatType;

template <typename DistanceTag>
void estimate_strong(
    const int k,  // the size of alphabet
    const std::vector<int>& center, // center string
    const int radius, // radius
    DistanceTag, // distance function
    const double alpha, // quantile of gaussian distribution
    const CountType iterations, // B: convergence test
    const bool is_only_u,
    const CountType nmax = std::numeric_limits<CountType>::max() // max # randomly generated strings
    ) {
  std::random_device rd;
  std::mt19937 mt(rd());

  std::uniform_int_distribution<int> agen(0, k-1);

  boost::math::normal gauss(0.0, 1.0);

  CountType s = center.size();
  std::vector<int> sn(center.size() + radius);
  std::vector<CountType> ur(radius+1);

  ur[0] = 1;
  for (int r = is_only_u ? radius : 1; r <= radius; r++) {
    int c = s - r;
    if (c < 0) c = 0;
    CountType urlsum = 0;
    CountType kl = 1;
    for (int i = 0; i < c; i++) {
      kl *= k;
    }
    for (CountType l = c; l <= s + r; l++, kl *= k) {
      sn.resize(l);
      CountType x = 0;
      CountType b = 0;
      CountType prevurln = 0;
      for (CountType n = 1; n < nmax; n++) {
        for (CountType i = 0; i < l; i++) {
          sn[i] = agen(mt);
        }
        auto d = measure(center, sn, DistanceTag());
        if (within(d, r, typename DistanceTag::CompareType())) {
          x++;
        }
        CountType urln = (CountType) (0.5 + kl * ((FloatType)x / (FloatType)n));
        if (urln == prevurln) {
          b++;
          if (b > iterations) {
            urlsum += urln;
            break;
          }
        } else {
          b = 0;
          prevurln = urln;
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

template <typename DistanceTag>
void estimate_confidence(
    const int k,  // the size of alphabet
    const std::vector<int>& center, // center string
    const int radius, // radius
    DistanceTag, // distance function
    const double alpha, // quantile of gaussian distribution
    const CountType lowerbound, // N: the lower bound of # randomly generated strings
    const bool is_only_u,
    const CountType nmax = std::numeric_limits<CountType>::max() // max # randomly generated strings
    ) {
  std::random_device rd;
  std::mt19937 mt(rd());

  std::uniform_int_distribution<int> agen(0, k-1);

  boost::math::normal gauss(0.0, 1.0);
  FloatType z = quantile(gauss, 1.0 - alpha * 0.5);
  FloatType z2 = z * z; // z^2
  FloatType z24 = z2 * 4; // 4z^2
  FloatType z4 = z2 * z2; // z^4

  CountType s = center.size();
  std::vector<int> sn(center.size() + radius);
  std::vector<CountType> ur(radius+1);
  ur[0] = 1;
  for (int r = is_only_u ? radius : 1; r <= radius; r++) {
    int c = s - r;
    if (c < 0) c = 0;
    CountType urlsum = 0;
    CountType kl = 1;
    for (int i = 0; i < c; i++) {
      kl *= k;
    }
    for (CountType l = c; l <= s + r; l++, kl *= k) {
      sn.resize(l);
      CountType x = 0;
      for (CountType n = 1; n < nmax; n++) {
        for (CountType i = 0; i < l; i++) {
          sn[i] = agen(mt);
        }
        auto d = measure(center, sn, DistanceTag());
        if (within(d, r, typename DistanceTag::CompareType())) {
          x++;
        }
        FloatType xn = (FloatType)x / (FloatType)n;
        CountType urln = (CountType) (0.5 + kl * xn);
        CountType e = sqrt(z24 * x * (1.0 - xn) + z4) * kl - z2;
        if (n > lowerbound and n > e) {
          urlsum += urln;
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

