/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * This header defines and  implements various functions to do hierarchical
 * shuffling of arrays, basically Owen-scrambling on the array indices.
 */
#ifndef SAMPLING_SHUFFLING_H
#define SAMPLING_SHUFFLING_H

#include <vector>

#include "utils.h"
#include "rng.h"

namespace sampling {
/*
 * GetShuffledIndices<base> is the core of the shuffling functionality.
 * Starting from an array of indices, e.g. {0, 1, 2, 3, ...}, it performs an
 * Owen-scrambling in the given base on those indices.
 *
 * See https://andrew-helmer.github.io/tree-shuffling/ for an explanation of
 * how it works in base-2. Fundamentally it calculates both an Owen-scrambled
 * van der Corput sequence, using the Stochastic Generation approach, and then
 * "reverses" the scrambling using an unscrambled van der Corput sequence.
 */
template <unsigned b>
int GetIntervalOffset(const uint32_t base_seed, const int interval,
                      const int pass) {
  uint32_t seed = CombineHashes(base_seed, Hash(interval));
  int strata_offset = (Permute(pass, b - 1, seed) + 1);
  if (pass > 0) strata_offset -= (Permute(pass - 1, b - 1, seed) + 1);
  return (b + strata_offset) % b;
}

template <unsigned base>
std::vector<int> GetShuffledIndices(const int length) {
  std::vector<int> randomized_indices(length);
  std::vector<int> digit_reversed_indices(length);
  RNG rng;

  digit_reversed_indices[0] = 0;
  randomized_indices[0] = rng.GetUniformInt(0, length);
  int interval_width = length / base;
  for (int prev_len = 1; prev_len < length; prev_len *= base) {
    ASSERT(prev_len*base <= length,
           "The length of the array to be shuffled needs to be a power "
           "of the base.");

    uint32_t offset_rnd_seed = rng.GetUniformInt();
    for (int pass = 0; pass < base-1; pass++) {
      for (int i = 0; i < prev_len; i++) {
        const int idx = prev_len*(pass+1) + i;
        const int prev_idx = idx - prev_len;
        int prev_interval = digit_reversed_indices[prev_idx] / interval_width;
        int next_interval = AddLSDigit<base>(prev_interval, 1);
        digit_reversed_indices[idx] = interval_width*next_interval;

        prev_interval = randomized_indices[prev_idx] / interval_width;
        int higher_interval = prev_interval / base;
        int interval_offset =
            GetIntervalOffset<base>(offset_rnd_seed, higher_interval, pass);
        next_interval = AddLSDigit<base>(prev_interval, interval_offset);
        randomized_indices[idx] =
            rng.GetUniformInt(0, interval_width) + interval_width*next_interval;
      }
    }
    interval_width /= base;
  }

  // We reindex the array by itself, using the randomized indices.
  for (int i = 0; i < length; i++) {
    digit_reversed_indices[i] = randomized_indices[digit_reversed_indices[i]];
  }
  return digit_reversed_indices;
}

// Declared here for a higher performance base-2 specialization in
/// shuffling.cpp.
template <>
std::vector<int> GetShuffledIndices<2>(const int length);

// Implements Progressive Shuffling from the paper. It only shuffles points
// that were generated in the same "pass" in their base.
template <int base>
std::vector<int> GetProgressiveShuffledIndices(const int num_samples) {
  std::vector<int> indices(num_samples);
  for (int i = 0; i < base; i++) indices[i] = i;
  for (int power = base; power < num_samples; power *= base) {
    for (int pass = 1; pass < base && (pass*power < num_samples); pass++) {
      ASSERT((pass+1)*power <= num_samples,
             "Progressive shuffling is only implemented for sequences "
             "that are a constant multiple of the base. Increase the sample "
             "count to " << ((pass+1)*power) << " samples.");
      std::vector<int> shuffled_indices = GetShuffledIndices<base>(power);
      for (int i = 0; i < power; i++) {
        indices[i + pass*power] = shuffled_indices[i] + pass*power;
      }
    }
  }
  return indices;
}

/*
 * This function will perform Progressive Shuffling on a set of N-dimensional
 * samples in a given base.
 */
template <int base>
void ProgressiveShuffleSamples(const int num_samples, const int nd,
                               double* samples) {
  // Copy samples into a temp vector.
  std::vector<double> temp_samples(num_samples*nd);
  for (int i = 0; i < num_samples*nd; i++) temp_samples[i] = samples[i];

  // Generates shuffled indices.
  std::vector<int> shuffled_indices =
      GetProgressiveShuffledIndices<base>(num_samples);

  // Writes out shuffled samples.
  for (int i = 0; i < num_samples; i++) {
    const int idx = shuffled_indices[i];
    for (int d = 0; d < nd; d++) samples[idx*nd + d] = temp_samples[i*nd + d];
  }
}

}  // namespace sampling

#endif  // SAMPLING_SHUFFLING_H
