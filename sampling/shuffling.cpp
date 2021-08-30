/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * This file just implements a specialized version of progressive shuffling for
 * base 2, which will be much faster and cleaner.
 */
#include "shuffling.h"

#include <vector>

#include "rng.h"
#include "utils.h"

namespace sampling {

template <>
std::vector<int> GetShuffledIndices<2>(const int length) {
  std::vector<int> randomized_indices(length);
  std::vector<int> bit_reversed_indices(length);
  RNG rng;
  bit_reversed_indices[0] = 0;
  randomized_indices[0] = rng.GetUniformInt(0, length);
  int interval_width = length / 2;
  for (int prev_len = 1; prev_len < length; prev_len *= 2) {
    for (int i = 0; i < prev_len && (prev_len+i) < length; i++) {
      bit_reversed_indices[i+prev_len] =
          (bit_reversed_indices[i] ^ interval_width);
      randomized_indices[i+prev_len] =
          (randomized_indices[i] ^ interval_width)
           ^ (rng.GetUniformInt(0, interval_width));
    }
    interval_width /= 2;  // Or interval_width >>= 1, if you prefer.
  }

  // We reindex the array by itself, using the randomized indices.
  for (int i = 0; i < length; i++) {
    bit_reversed_indices[i] = randomized_indices[bit_reversed_indices[i]];
  }
  return bit_reversed_indices;
}

}  // namespace sampling
