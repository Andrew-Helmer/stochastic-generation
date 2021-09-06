/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * This file implements the "generate_samples" command-line utility.
 * Running "make generate_samples" will build the generate_samples command-line
 * util, which can be used to stochastically generate different sample
 * sequences.
 *
 * The required arguments are:
 * --algorithm= or --a= i.e. which sampler to use
 *    Valid algorithms are:
 *    { ssobol, pmj02, shalton, sfaure(03|05|07|011), uniform }
 * --n= is "number of samples"
 * --nd= is "number of dimensions"
 *
 * e.g.
 *    ./generate_samples --a=shalton --n=500 --nd=5
 *
 * Optional arguments are:
 * --owen will use full owen-scrambling rather than correlated swapping for the
 *    sfaure and shalton sequences.
 * --shuffle to perform progressive shuffling, which
 *    will decorrelate multiple runs of the sequence, while maintaining good
 *    progressive properties.
 * --bn2d will use best-candidate sampling to improve minimum distances for the
 *    first two dimensions only. Note that for the higher base sequences, it
 *    will NOT select the strata or delta-strata based on minimum distances.
 *    That choice is still random.
 */
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "sampling/bn_utils.h"
#include "sampling/ssobol.h"
#include "sampling/sfaure.h"
#include "sampling/shalton.h"
#include "sampling/utils.h"

using std::string;

namespace {
// Print out n-dimensional floating point arrays as tuples, one point per line.
// For example, a sequence of N 3D points would be...
// (x_0, y_0, z_0)
// (x_1, y_1, z_1)
// (x_2, y_2, z_2)
// ...
// (x_(n-1), y_(n-1), z_(n-1))
void printNdPoints(int num_samples, int nd, double* samples) {
  for (int i = 0; i < num_samples; i++) {
    printf("(");
    for (int d = 0; d < nd; d++) {
      if (d > 0) printf(", ");
      printf("%.17g", samples[i*nd + d]);
    }
    printf(")\n");
  }
}

// If the string "arg" starts with the string "arg_name", take the rest of
// arg and put it in *val.
inline void maybeGetStringArg(
    const char* arg_name, const char* arg, string* val) {
  const int len = strlen(arg_name);
  if (strncmp(arg_name, arg, len) == 0) {
    *val = string(arg+len);
  }
}

// If the string "arg" starts with the string "arg_name", take the rest of
// arg, parse it as an integer, and put it in *val.
inline void maybeGetIntArg(const char* arg_name, const char* arg, int* val) {
  const int len = strlen(arg_name);
  if (strncmp(arg_name, arg, len) == 0) {
    *val = atoi(arg+len);
  }
}

// If the string "arg" starts with the string "arg_name", set *val to true.
inline void maybeGetBoolArg(const char* arg_name, const char* arg, bool* val) {
  const int len = strlen(arg_name);
  if (strncmp(arg_name, arg, len) == 0) {
    *val = true;
  }
}
}  // namespace

int main(int argc, const char** argv) {
  typedef void (*sample_fn)(int, int, bool, int, bool, double*);
  static const std::unordered_map<string, sample_fn> kSeqMap = {
    {"ssobol",    &sampling::GetStochasticSobolSamples},
    {"pmj02",       &sampling::GetPMJ02Samples},
    {"sfaure03",  &sampling::GetStochasticFaure03Samples},
    {"sfaure05",  &sampling::GetStochasticFaure05Samples},
    {"sfaure07",  &sampling::GetStochasticFaure07Samples},
    {"sfaure011", &sampling::GetStochasticFaure011Samples},
    {"shalton",   &sampling::GetStochasticHaltonSamples},
  };

  argc--;
  argv++;

  string sequence;
  int num_samples = -1;
  int num_dims = -1;
  bool owen = false;
  bool shuffle = false;
  bool bn2d = false;
  for (int i = 0; i < argc; i++) {
    maybeGetStringArg("--seq=", argv[i], &sequence);
    maybeGetIntArg("--n=", argv[i], &num_samples);
    maybeGetIntArg("--nd=", argv[i], &num_dims);
    maybeGetBoolArg("--owen", argv[i], &owen);
    maybeGetBoolArg("--shuffle", argv[i], &shuffle);
    maybeGetBoolArg("--bn2d", argv[i], &bn2d);
  }

  ASSERT(num_samples > 0,
         "Need to set --n to a value greater than zero.");
  ASSERT(num_dims >= 1,
         "Need to set --nd to a value >= 1.");
  ASSERT(!bn2d || num_dims >= 2,
         "Need at least two dimensions to use --bn2d flag.");

  const int num_candidates = bn2d ? sampling::kBestCandidateSamples : 1;

  auto find_seq = kSeqMap.find(sequence);
  if (find_seq != kSeqMap.end()) {
    sample_fn get_samples = find_seq->second;
    auto samples = std::make_unique<double[]>(num_samples*num_dims);
    auto start = std::chrono::high_resolution_clock::now();
    (*get_samples)(num_samples, num_dims, shuffle, num_candidates, owen,
                   samples.get());
    auto end = std::chrono::high_resolution_clock::now();
    const double milliseconds =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start)
            .count() / 1000.0;
    std::cerr.precision(4);
    std::cerr << num_samples << " points generated in "
              << milliseconds << "ms\n";
    printNdPoints(num_samples, num_dims, samples.get());

    return 0;
  }

  printf("%s isn't a valid sequence.\n", sequence.c_str());
  return 1;
}
