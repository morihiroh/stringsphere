#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>
#include "stringdistance.hpp"
#include "stringspheresearch.hpp"

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
    ("distance,d", po::value<int>()->default_value(2), ": 0: Extended Hamming 1: Longest common subsequence 2: Levenshtein 3: Damerau-Levenshtein")
    ("k,k", po::value<int>()->default_value(2), ": the size of alphabet A")
    ("string,s", po::value<int>()->default_value(3), ": the length of center string")
    ("radius,r", po::value<int>()->default_value(2), ": radius")
    ("method,m", po::value<int>()->default_value(0), ": 0: random selection method 1: exhaustive search method")
    ("iterations,i", po::value<int>()->default_value(10), ": least # of iterations")
    ("lomit,l", po::value<int>()->default_value(-1), ": omit smaller l (-1: no omit)")
    ("onlyu,u", ": estimate only u of the radius")
    ("quiet,q", ": display only results")
    ("help,h", ": show this help message");

  po::variables_map argmap;
  po::store(po::parse_command_line(argc, argv, opt), argmap);
  po::notify(argmap);

  if (argmap.count("help")) { 
    cerr << "Usage: " << argv[0] << " [option]" << endl << opt << endl;

    return EXIT_FAILURE;
  }

  int method = argmap["method"].as<int>();
  int distancetype = argmap["distance"].as<int>();
  int k = argmap["k"].as<int>();
  std::vector<int> center(argmap["string"].as<int>(), 0);
  int radius = argmap["radius"].as<int>();
  CountType iterations = argmap["iterations"].as<int>();
  int lomit = argmap["lomit"].as<int>();
  bool is_only_u = (method == 0) and argmap.count("onlyu");
  bool is_quiet = argmap.count("quiet");

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
    }
  };

  auto countfunc = [&](auto dist) {
    switch (method) {
      case 0:
        estimate(k, center, radius, dist, iterations, lomit, is_only_u);
        break;
      case 1:
        searchall(k, center, radius, dist);
        break;
    } 
  };

  if (not is_quiet) {
    switch (method) {
      case 0:
        std::cout << "Random selection method";
        break;
      case 1:
        std::cout << "Exhaustive search method";
        break;
    }
    std::cout << " with ";
    switch (distancetype) {
      case 0:
        std::cout << "extended Hamming distance\n";
        break;
      case 1:
        std::cout << "longest common subsequence\n";
        break;
      case 2:
        std::cout << "Levenshtein distance\n";
        break;
      case 3:
        std::cout << "Damerau-Levenshtein distance\n";
        break;
    }
    std::cout << "k\ts\tr\tu";
    if (not is_only_u) {
      std::cout << "\tv";
    }
    std::cout << "\n";
  }

  swfunc(distancetype, countfunc);

  return EXIT_SUCCESS;
}

