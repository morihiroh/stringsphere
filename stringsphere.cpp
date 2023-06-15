#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>
#include <chrono>
#include "stringdistance.hpp"
#include "stringspheresearch.hpp"

int main(int argc, char **argv)
{
  using namespace std;
  namespace po = boost::program_options;
  po::options_description opt("Allowed options");

  opt.add_options()
    ("distance,d", po::value<int>()->default_value(2), ": 0: Extended Hamming 1: Longest common subsequence 2: Levenshtein 3: Damerau-Levenshtein")
    ("k,k", po::value<int>()->default_value(2), ": the size |A| of alphabet A")
    ("string,s", po::value<int>()->default_value(0), ": the center string in decimal number. Ex) 3 means 011 in |A|=2 and the length = 3.")
    ("size,z", po::value<int>()->default_value(3), ": the length (size) of center string")
    ("radius,r", po::value<int>()->default_value(2), ": radius")
    ("method,m", po::value<int>()->default_value(2), std::string(": 1: exhaustive search method 2: random selection with absolute error 3: random selection with relative error").c_str())
    ("gamma,g", po::value<double>()->default_value(1.0), ": specify gamma >= 1")
    ("epsilon,p", po::value<double>()->default_value(0.1), ": specify epsilon > 0")
    ("iterations,i", po::value<CountType>()->default_value(10), ": least # of iterations")
    ("maxiter,j", po::value<CountType>()->default_value(0), ": max # of iterations (0: unlimited)")
    ("ell,l", po::value<int>()->default_value(-1), ": specify ell (max{s-r,0}<=ell<=s+r) (-1: all)")
    ("onlyu,u", ": estimate only u of the radius")
    ("elapsed,e", ": print the elapsed time in milliseconds")
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
  double gamma = argmap["gamma"].as<double>();
  double epsilon = argmap["epsilon"].as<double>();
  int distancetype = argmap["distance"].as<int>();
  int k = argmap["k"].as<int>();
  std::vector<int> center(argmap["size"].as<int>(), 0);
  radix_convert(center, k, argmap["string"].as<int>());
  int radius = argmap["radius"].as<int>();
  CountType iterations = argmap["iterations"].as<CountType>();
  CountType maxiter = argmap["maxiter"].as<CountType>();
  int ell = argmap["ell"].as<int>();
  bool is_only_u = (method == 0) and (argmap.count("onlyu") or (ell >= 0));
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
      case 1:
        searchall(k, center, radius, dist);
        break;
      case 2:
      case 3:
        estimate(method, gamma, epsilon, k, center, radius, dist, iterations, maxiter, ell, is_only_u);
        break;
    } 
  };

  if (not is_quiet) {
    switch (method) {
      case 1:
        std::cout << "Exhaustive search method";
        break;
      case 2:
        std::cout << "Random selection method with absolute error";
        break;
      case 3:
        std::cout << "Random selection method with relative error";
        break;
    }
    std::cout << " in ";
    switch (distancetype) {
      case 0:
        std::cout << "extended Hamming distance";
        break;
      case 1:
        std::cout << "longest common subsequence";
        break;
      case 2:
        std::cout << "Levenshtein distance";
        break;
      case 3:
        std::cout << "Damerau-Levenshtein distance";
        break;
    }
    std::cout << " for center s = \"";
    for (int i = center.size() - 1; i >= 0; i--) {
      std::cout << center[i];
    }
    std::cout << "\"\n|A|\t|s|\tradius\tu";
    if (not is_only_u) {
      std::cout << "\tv";
    }
    std::cout << "\n";
  }

  auto start = std::chrono::system_clock::now();

  swfunc(distancetype, countfunc);
  
  if (argmap.count("elapsed")) {
    auto end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    std::cout << elapsed << "\n";
  }

  return EXIT_SUCCESS;
}

