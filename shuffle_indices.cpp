/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * Running "make shuffle_indices" will build a shuffle_indices command-line
 * util, which can output arrays of shuffled indices that can shuffle a base-b
 * (t,s)-sequence into another base-b (t,s)-sequence, maintaining its
 * progressive properties, but decorrelating it from other runs. See Section 5.2
 * of the paper.
 *
 * The required arguments are:
 * --n= is the "number of samples"
 * --b= or --base= is the base to perform the shuffling in.
 *
 * e.g.
 *    ./shuffle_indices --n=9 --b=3
 *
 * Optional arguments are:
 * --prog or --progressive will perform progressive shuffling rather than full
 *    shuffling. This will fully decorrelate stochastic sequences, but
 *    power-of-b prefixes of the original sequence will be maintained.
 *
 * Note that if --progressive is not used, --n must be a power-of-b, n=b^m for
 * some m. If --progressive is used, it can be a constant multiple of a
 * power-of-b, n = k*b^m, where 1 <= k <= b.
 */
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "sampling/shuffling.h"
#include "sampling/utils.h"

using std::string;

namespace {
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
  argc--;
  argv++;

  int num_samples = -1;
  int base = -1;
  bool progressive = false;
  for (int i = 0; i < argc; i++) {
    maybeGetIntArg("--n=", argv[i], &num_samples);
    maybeGetIntArg("--b=", argv[i], &base);
    maybeGetIntArg("--base=", argv[i], &base);
    maybeGetBoolArg("--prog", argv[i], &progressive);
    maybeGetBoolArg("--progressive", argv[i], &progressive);
  }

  ASSERT(num_samples > 0,
         "--n must be set to a value greater than zero.");
  ASSERT(base > 0,
         "--b must be set to a value greater than zero.");
  std::vector<int> shuffled_indices;
  if (progressive) {
    switch (base) {
      case 2:
        shuffled_indices =
            sampling::GetProgressiveShuffledIndices<2>(num_samples);
        break;
      case 3:
        shuffled_indices =
            sampling::GetProgressiveShuffledIndices<3>(num_samples);
        break;
      case 5:
        shuffled_indices =
            sampling::GetProgressiveShuffledIndices<5>(num_samples);
        break;
      case 7:
        shuffled_indices =
            sampling::GetProgressiveShuffledIndices<7>(num_samples);
        break;
      case 11:
        shuffled_indices =
            sampling::GetProgressiveShuffledIndices<11>(num_samples);
        break;
      default:
        ASSERT(false, "base " << base << " not supported. "
                      "Must be 2, 3, 5, 7, or 11.");
        return 1;
    }
  } else {
    switch (base) {
      case 2:
        shuffled_indices = sampling::GetShuffledIndices<2>(num_samples);
        break;
      case 3:
        shuffled_indices = sampling::GetShuffledIndices<3>(num_samples);
        break;
      case 5:
        shuffled_indices = sampling::GetShuffledIndices<5>(num_samples);
        break;
      case 7:
        shuffled_indices = sampling::GetShuffledIndices<7>(num_samples);
        break;
      case 11:
        shuffled_indices = sampling::GetShuffledIndices<11>(num_samples);
        break;
      default:
        ASSERT(false, "base " << base << " not supported. "
                      "Must be 2, 3, 5, 7, or 11.");
        return 1;
    }
  }

  for (int i = 0; i < num_samples; i++) {
    std::cout << shuffled_indices[i] << " ";
  }
  std::cout << std::endl;
  return 1;
}
