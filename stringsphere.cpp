#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>
#include "stringdistance.hpp"
#include <typeinfo>

std::vector<int> parsestring(const std::string& s)
{
  std::vector<std::string> t;
  boost::algorithm::split(t, s, boost::is_any_of(","));
  std::vector<int> res(t.size());
  for (size_t i = 0; i < t.size(); i++) {
    res[i] = boost::lexical_cast<int>(t[i]);
  }
  return res;
}

int main(int argc, char **argv)
{
  using namespace std;
  namespace po = boost::program_options;
  po::options_description opt("Allowed options");

  opt.add_options()
    ("distance,d", po::value<int>()->default_value(2), ": 0: Extended Hamming 1: Longest common subsequence 2: Levenshtein 3: Damerau-Levenshtein"
#ifndef DISABLE_JAROWINKLER
" 4: Jaro-Winkler"
#endif
)
    ("k,k", po::value<int>()->default_value(2), ": the size of alphabet A")
    ("string,s", po::value<int>()->default_value(3), ": the length of center string")
    ("radius,r", po::value<int>()->default_value(2), ": radius")
    ("method,m", po::value<int>()->default_value(0), ": 0: strong consistent estimate 1: confidence interval estimate")
    ("iterations,i", po::value<int>()->default_value(100), ": # of iterations (B) that estimated values are equal (convergence test).")
    ("alpha,a", po::value<double>()->default_value(0.01), ": 100xalpha/2 percent point of gaussian")
    ("lowerbound,l", po::value<int>()->default_value(10000), ": the lower bound (N) of # randomly generated strings")
    ("exhaustive,e", ": search all instead of the stochastic method")
    ("onlyu,u", ": find only u")
    ("quiet,q", ": display only results")
    ("help,h", ": show this help message");

  po::variables_map argmap;
  po::store(po::parse_command_line(argc, argv, opt), argmap);
  po::notify(argmap);

  if (argmap.count("help")) {
    cerr << "Usage: " << argv[0] << " [option]" << endl << opt << endl;

    return EXIT_FAILURE;
  }

  int distancetype = argmap["distance"].as<int>();
  int method = argmap["method"].as<int>();
  double alpha = argmap["alpha"].as<double>();
  int k = argmap["k"].as<int>();
  std::vector<int> center(argmap["string"].as<int>(), 0);
  int radius = argmap["radius"].as<int>();
  CountType lowerbound = argmap["lowerbound"].as<int>();
  CountType iterations = argmap["iterations"].as<int>();
  bool is_only_u = argmap.count("onlyu");

  auto swfunc = [](int distancetype, auto func) {
    switch (distancetype) {
      case 0:
        func(ExtendedHammingDistance());
        break;
      case 1:
        func(LongestCommonSubsequence());
        break;
      case 2:
        func(LevenshteinDistance());
        break;
      case 3:
        func(DamerauLevenshteinDistance());
        break;
#ifndef DISABLE_JAROWINKLER
      case 4:
        func(JaroWinklerDistance());
        break;
#endif
    }
  };

  auto countfunc = [&](auto dist) {
    if (argmap.count("exhaustive")) {
      searchall(k, center, radius, dist);
    } else if (method == 0) {
      estimate_strong(k, center, radius, dist, alpha, iterations, is_only_u);
    } else if (method == 1) {
      estimate_confidence(k, center, radius, dist, alpha, lowerbound, is_only_u);
    }
  };

  if (not argmap.count("quiet")) {
    std::cout << "k\ts\tr\tu\tv\n";
  }
  if (radius == 0) {
    std::cout << k << "\t" << center.size() << "\t" << 0 << "\t" << 1 << "\t" << 1 << "\n";
  } else {
    swfunc(distancetype, countfunc);
  }

  return EXIT_SUCCESS;
}

