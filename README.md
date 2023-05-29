# stringsphere

1. Description

stringsphere is a softwawe implemented in C++ for computing estimates of the volumes of balls and spheres of strings in the metric space of all strings composed of a given alphabet provided with a specified type of edit distance.

Please see the following paper for how and why this softwawe was developed.

> Hitoshi Koyano and Morihiro Hayashida, Volume formula and growth rates of the balls of strings under the edit distances, submitted.

In this paper, the volume formula for the balls of strings under the extended Hamming distance was derived. 
Furthermore, a randomized algorithm was formulated for computing estimates of the volumes of balls of strings efficiently, and by performing numerical experiments using the algorithm, a conjecture was presented on the volume formula for the balls of strings under the Lavenshtein distance.

stringsphere computes estimates of the volumes of balls and spheres of strings using the abovementioned randomized algorithm.

2. Compile

```sh
g++ -O3 -lboost_program_options stringsphere.cpp
```

3. Usage

stringsphere -h (display help message)

## Example

$u_i(x^2,\sum k)$
