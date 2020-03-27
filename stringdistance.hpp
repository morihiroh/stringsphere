#pragma once
#include <vector>
#include <algorithm>
#include <boost/multi_array.hpp>

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
