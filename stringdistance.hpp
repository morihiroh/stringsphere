#pragma once
#include <string>
#include <vector>
#include <random>
#include <boost/math/distributions/normal.hpp>
#include <boost/multi_array.hpp>

#define DISABLE_JAROWINKLER

struct StringDistance {
  virtual int operator()(const std::vector<int>& x, const std::vector<int>& y) = 0;
};

struct ExtendedHammingDistance : public StringDistance {
    int operator()(const std::vector<int>& x, const std::vector<int>& y) {
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
};

struct LongestCommonSubsequenceDistance : public StringDistance {
    int operator()(const std::vector<int>& x, const std::vector<int>& y) {
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
};

struct LevenshteinDistance : public StringDistance {
    int operator()(const std::vector<int>& x, const std::vector<int>& y) {
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
};

struct DamerauLevenshteinDistance : public StringDistance {
    int operator()(const std::vector<int>& x, const std::vector<int>& y) {
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
};

#ifndef DISABLE_JAROWINKLER
struct JaroWinklerDistance {
  double operator()(const std::vector<int>& x, const std::vector<int>& y) {
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
};
#endif

#if 0
static void print(const std::vector<int>& s) {
  for (auto c : s) {
    std::cout << c;
  }
  std::cout << "\n";
}
#endif

template <typename DistancePred>
void count(int& tildeu, int& tildev, // estimated values
    int& hatu, int& hatv, // strong eastimated values
    const int k,  // the size of alphabet
    const std::vector<int>& center, // center string
    const int radius, // radius
    DistancePred distfunc, // distance function
    const double alpha, // quantile of gaussian distribution
    const int lowerbound, // the lower bound of # randomly generated strings
    const int iterations, // convergence test
    const int nmax = INT_MAX  // max # randomly generated strings
    ) {
  std::random_device rd;
  std::mt19937 mt(rd());
  int lstart = center.size() - radius;
  int lend = center.size() + radius;

  double kmc = std::pow(k, lend) - ((lstart >= 1) ? std::pow(k, lstart - 1) : 0.0); // k^(|s|+r) - c

  if (lstart < 0) {
    lstart = 0;
  }
  std::uniform_int_distribution<int> lgen(lstart, lend);
  std::uniform_int_distribution<int> agen(0, k-1);

  boost::math::normal gauss(0.0, 1.0);
  double z = quantile(gauss, 1.0 - alpha * 0.5);
  double z2 = z * z; // z^2
  double z24 = z2 * 4; // 4z^2
  double z4 = z2 * z2; // z^4
 
  //print(center);

  double xn = 0;
  double yn = 0;
  bool uflag = true;
  bool vflag = true;
  int bu = iterations;
  int bv = iterations;
  int preun = INT_MAX;
  int prevn = INT_MAX;
  std::vector<int> sn(center.size() + radius);
  for (int n = 1; (bu > 0 or bv > 0 or uflag or vflag) and n < nmax; n++) {
    int l = lgen(mt);
    sn.resize(l);
    for (int i = 0; i < l; i++) {
      sn[i] = agen(mt);
    }
    //print(sn);
    auto d = distfunc(center, sn);
    //std::cout << "d: " << d << "\n";
    if (d <= radius) {
      xn++;
      if (d == radius) {
        yn++;
      }
    }
    double tmp = kmc / n;
    int un = (int)(0.5 + tmp * xn);
    int vn = (int)(0.5 + tmp * yn);
    if (uflag) {
      if ((n > lowerbound) and (n > kmc * sqrt(z24 * xn * (1.0 - xn/n) + z4) - z2)) {
        tildeu = un;
        uflag = false;
      }
    }
    if (vflag) {
      if ((n > lowerbound) and (n > kmc * sqrt(z24 * yn * (1.0 - yn/n) + z4) - z2)) {
        tildev = vn;
        vflag = false;
      }
    }
    if (bu > 0) {
      if (un == preun) {
        bu--;
        if (bu == 0) {
          hatu = un;
        }
      } else {
        bu = iterations;
      }
    }
    if (bv > 0) {
      if (vn == prevn) {
        bv--;
        if (bv == 0) {
          hatv = vn;
        }
      } else {
        bv = iterations;
      }
    }
    preun = un;
    prevn = vn;
  }
  if (bu > 0 or bv > 0 or uflag or vflag) {
    std::cerr << "# iterations reached max\n";
  }
}

template <typename DistancePred>
void searchall(int& u, int& v, // estimated values
    const int k,  // the size of alphabet
    const std::vector<int>& center, // center string
    const int radius, // radius
    DistancePred distfunc, // distance function
    //StringDistance& distfunc,
    const double alpha, // quantile of gaussian distribution
    const int nmax = INT_MAX  // max # randomly generated strings
    ) {
  std::vector<int> sn(center.size() + radius);
  int xn = 0;
  int yn = 0;
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
      //print(sn);
      int d = distfunc(center, sn);
      //std::cout << "d: " << d << "\n";
      if (d <= radius) {
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

  u = xn;
  v = yn;
}

