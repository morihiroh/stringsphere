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
    ("alpha,a", po::value<double>()->default_value(0.01), ": 100xalpha/2 percent point of gaussian")
    ("k,k", po::value<int>()->default_value(2), ": the size of alphabet A")
    ("center,c", po::value<std::string>()->default_value("0,0,0"), ": center string in integers from 0 to k-1 with comma delimiter")
    ("radius,r", po::value<int>()->default_value(2), ": radius")
    ("lowerbound,l", po::value<int>()->default_value(10000), ": the lower bound of # randomly generated strings")
    ("iterations,i", po::value<int>()->default_value(100), ": # of iterations that estimated values are equal (convergence test).")
    ("searchall,s", ": search all instead of the stochastic method")
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
  double alpha = argmap["alpha"].as<double>();
  int k = argmap["k"].as<int>();
  auto center = parsestring(argmap["center"].as<std::string>());
  int radius = argmap["radius"].as<int>();
  int lowerbound = argmap["lowerbound"].as<int>();
  int iterations = argmap["iterations"].as<int>();

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

  if (argmap.count("searchall")) {
    int u, v;
    auto searchallfunc = [&](auto dist) {
        searchall(u, v, k, center, radius, dist, alpha);
    };

    swfunc(distancetype, searchallfunc);

    if (argmap.count("quiet")) {
      std::cout << u << '\t' << v << '\t';
    } else {
      std::cout << "u = " << u << "\tv = " << v << "\n";
    }
  } else {
    int tildeu, tildev, hatu, hatv;
    auto countfunc = [&](auto dist) {
      count(tildeu, tildev, hatu, hatv, k, center, radius, dist, alpha, lowerbound, iterations);
    };

    swfunc(distancetype, countfunc);

    if (argmap.count("quiet")) {
      std::cout << tildeu << "\t" << tildev << "\t" << hatu << "\t" << hatv << "\t";
    } else {
      std::cout << "~u = " << tildeu << "\t~v = " << tildev << "\t ^u = " << hatu << "\t ^v = " << hatv << "\n";
    }
  }

  return EXIT_SUCCESS;
}

