# stringsphere

1. Description

stringsphere is a softwawe implemented in C++ for computing estimates of the volumes of balls and spheres of strings in the metric space of all strings composed of a given alphabet provided with a specified type of edit distance.

Please see the following paper for how and why this softwawe was developed.

```
Hitoshi Koyano and Morihiro Hayashida, Volume formula and growth rates of the balls of strings under the edit distances, submitted.
```

In this paper, the volume formula for the balls of strings under the extended Hamming distance was derived. 
Furthermore, a randomized algorithm was formulated for computing estimates of the volumes of balls of strings efficiently, and by performing numerical experiments using the algorithm, a conjecture was presented on the volume formula for the balls of strings under the Lavenshtein distance.

stringsphere computes estimates of the volumes of balls and spheres of strings using the abovementioned randomized algorithm.

2. Compile

```
g++ -O3 -lboost_program_options stringsphere.cpp
```

It can be compiled by g++ version 12.3 ([GNU C++ compiler](https://gcc.gnu.org)) with [boost program options library](https://www.boost.org).

3. Usage

We can use the following options to this software.

 - Options

```
  -d [ --distance ] arg (=2)    : 0: Extended Hamming 1: Longest common
                                subsequence 2: Levenshtein 3:
                                Damerau-Levenshtein
  -k [ --k ] arg (=2)           : the size |A| of alphabet A
  -s [ --string ] arg (=0)      : the center string in decimal number. Ex) 3
                                means 011 in |A|=2 and the length = 3.
  -z [ --size ] arg (=3)        : the length (size) of center string
  -r [ --radius ] arg (=2)      : radius
  -m [ --method ] arg (=0)      : 0: random selection method 1: exhaustive
                                search method
  -i [ --iterations ] arg (=10) : least # of iterations
  -j [ --maxiter ] arg (=0)     : max # of iterations (0: unlimited)
  -l [ --ell ] arg (=-1)        : specify ell (max{s-r,0}<=ell<=s+r) (-1: all)
  -u [ --onlyu ]                : estimate only u of the radius
  -e [ --elapsed ]              : print the elapsed time in milliseconds
  -q [ --quiet ]                : display only results
  -h [ --help ]                 : show this help message
```

 - Example

```
stringsphere -m 0 -k 2 -z 5 -s 3 -r 3

Random selection method with Levenshtein distance for center s = 00011
|A|     |s|     radius  u       v
2       5       0       1       1
2       5       1       15      14
2       5       2       84      69
2       5       3       273     189
```

This example shows the results of $U_{2,d_L}(00011,3)$ and $V_{2,d_L}(00011,3)$ by the randomized algorithm on $A=\\{0,1\\}$ for $\phi=00011$. 

