/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * This file implements Stochastic Sobol' and pmj02 sampling, along with
 * optional best-candidate sampling (for the first two dimensions).
 */
#include "ssobol.h"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

#include "bn_utils.h"
#include "rng.h"
#include "shuffling.h"
#include "utils.h"
#include "xor_values.h"

namespace sampling {
namespace {
inline void GetSobolPoint(const std::vector<int>& strata,
                          const double i_strata,
                          const int nd,
                          RNG* rng,
                          double* sample) {
  for (int d = 0; d < nd; d++) {
    sample[d] = (rng->GetUniformFloat() + strata[d]) * i_strata;
  }
}
}  // namespace

void GetStochasticSobolSamples(const int num_samples,
                               const int nd,
                               const bool shuffle,
                               const int n_candidates,
                               const bool owen,
                               double* samples) {
  ASSERT(nd <= MAX_SOBOL_DIM, "Stochastic Sobol' only works up to "
                                  << MAX_SOBOL_DIM << "dimensions.");
  ERRIF(owen, "--owen flag is meaningless for ssobol sequence.");

  RNG rng;
  // Generate first sample randomly.
  for (int d = 0; d < nd; d++) samples[d] = rng.GetUniformFloat();

  SampleGrid2D sample_grid;
  if (n_candidates > 1) sample_grid.AddSample({samples[0], samples[1]});
  // Used to temporarily store candidate positions for each new point.
  std::vector<SampleGrid2D::Point2D> candidates(n_candidates);

  // Used to temporarily store the strata for each dimension, for a new point.
  std::vector<int> strata(nd);
  for (int log_n = 0; (1 << log_n) < num_samples; log_n++) {
    int prev_len = 1 << log_n;
    int n_strata = prev_len * 2;
    double i_strata = 1.0 / n_strata;
    for (int i = 0; i < prev_len && (prev_len + i) < num_samples; i++) {
      for (int d = 0; d < nd; d++) {
        // It may be better if the sobol_xors were transposed for this, so that
        // the inner loop corresponded to contiguous array values.
        int prev_idx = i ^ sobol_xors[d][log_n];
        int prev_stratum = samples[prev_idx * nd + d] * n_strata;
        strata[d] = prev_stratum^1;
      }

      double* sample = &(samples[(prev_len + i) * nd]);
      GetSobolPoint(strata, i_strata, nd, &rng, sample);

      // Do best-candidate sampling.
      if (n_candidates > 1) {
        candidates[0] = {sample[0], sample[1]};
        for (int cand_idx = 1; cand_idx < n_candidates; cand_idx++) {
          double candidate[2];
          GetSobolPoint(strata, i_strata, 2, &rng, candidate);
          candidates[cand_idx] = {candidate[0], candidate[1]};
        }
        int best_candidate = sample_grid.GetBestCandidate(candidates);
        sample_grid.AddSample(candidates[best_candidate]);
        sample[0] = candidates[best_candidate].x;
        sample[1] = candidates[best_candidate].y;
      }
    }
  }

  if (shuffle) ProgressiveShuffleSamples<2>(num_samples, nd, samples);
}

namespace {
inline void GetPMJ02Point(const int x_stratum,
                          const int y_stratum,
                          const double i_strata,
                          RNG* rng,
                          double* sample) {
  sample[0] = (rng->GetUniformFloat() + x_stratum) * i_strata;
  sample[1] = (rng->GetUniformFloat() + y_stratum) * i_strata;
}
}  // namespace

void GetPMJ02Samples(const int num_samples,
                     const int nd,
                     const bool shuffle,
                     const int n_candidates,
                     const bool owen,
                     double* samples) {
  ASSERT(nd == 2, "PMJ02 only works for 2 dimensions. Use ssobol instead.");
  ERRIF(owen, "--owen flag is meaningless for pmj02 sequence.");

  RNG rng;
  // Generate first sample randomly.
  for (int d = 0; d < 2; d++) samples[d] = rng.GetUniformFloat();

  SampleGrid2D sample_grid;
  if (n_candidates > 1) sample_grid.AddSample({samples[0], samples[1]});
  // Used to temporarily store candidate positions for each new point.
  std::vector<SampleGrid2D::Point2D> candidates(n_candidates);

  for (int log_n = 0; (1 << log_n) < num_samples; log_n++) {
    int prev_len = 1 << log_n;
    int n_strata = prev_len * 2;
    double i_strata = 1.0 / n_strata;
    for (int i = 0; i < prev_len && (prev_len + i) < num_samples; i++) {
      const int prev_x_idx = i ^ pmj02_xors[0][log_n];
      const int prev_x_stratum = samples[prev_x_idx*2] * n_strata;
      const int x_stratum = prev_x_stratum^1;

      const int prev_y_idx = i ^ pmj02_xors[1][log_n];
      const int prev_y_stratum = samples[prev_y_idx*2 + 1] * n_strata;
      const int y_stratum = prev_y_stratum^1;

      double* sample = &(samples[(prev_len + i) * 2]);
      GetPMJ02Point(x_stratum, y_stratum, i_strata, &rng, sample);

      // Do best-candidate sampling.
      if (n_candidates > 1) {
        candidates[0] = {sample[0], sample[1]};
        for (int cand_idx = 1; cand_idx < n_candidates; cand_idx++) {
          double candidate[2];
          GetPMJ02Point(x_stratum, y_stratum, i_strata, &rng, candidate);
          candidates[cand_idx] = {candidate[0], candidate[1]};
        }
        int best_candidate = sample_grid.GetBestCandidate(candidates);
        sample_grid.AddSample(candidates[best_candidate]);
        sample[0] = candidates[best_candidate].x;
        sample[1] = candidates[best_candidate].y;
      }
    }
  }

  if (shuffle) ProgressiveShuffleSamples<2>(num_samples, nd, samples);
}

/*
 * From here on down are implementations of the "stateless" scrambled Sobol'
 * function from Section 5.3 of the paper.
 */
namespace {
/*
 * randfloat() from Correlated Multi-Jittered Sampling by Andrew Kensler (2013).
 */
inline uint32_t HashToRndInt(int index, unsigned seed) {
  unsigned result = index;
  result ^= seed;
  result ^= result >> 17;
  result ^= result >> 10;  result *= 0xb36534e5;
  result ^= result >> 12;
  result ^= result >> 21;  result *= 0x93fc4795;
  result ^= 0xdf6e307f;
  result ^= result >> 17;  result *= 1 | seed >> 18;
  return result;
}
inline double HashToRnd(int index, unsigned seed) {
  return static_cast< double >(HashToRndInt(index, seed)) / 4298115584.0;
}

/*
 * Gets the index of the most significant bit, or -1 if the value is zero.
 * Copied from https://stackoverflow.com/a/31718095/624250 with the addition of
 * a specific check for zero.
 */
inline int GetMSB(uint32_t v) {
  if (v == 0) return -1;

  static constexpr int MultiplyDeBruijnBitPosition[32] = {
      0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
      8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
  };

  // First round down to one less than a power of 2.
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;

  return MultiplyDeBruijnBitPosition[(uint32_t)( v * 0x07C4ACDDU ) >> 27] - 1;
}
}  // namespace

double GetSobolStateless(int idx, int dim, uint32_t seed, int nd) {
  assert(nd <= MAX_SOBOL_DIM);
  // Base case, return first randomly placed point.
  if (!idx) return HashToRnd(dim, seed);

  // Determine stratum size and place in previous strata.
  int log_n = GetMSB(idx);  // (Right-most bit is numbered zero)
  int prev_len = 1 << log_n;
  int n_strata = prev_len * nd;
  int i = idx - prev_len;

  // Recursively get stratum of previous sample.
  int prev_stratum =
      GetSobolStateless(i ^ sobol_xors[dim][log_n], dim, seed) * n_strata;
  // Generate new sample in adjacent stratum.
  return ((prev_stratum ^ 1) + HashToRnd(idx * nd + dim, seed)) / n_strata;
}

double GetSobolStatelessIter(int idx, int dim, uint32_t seed, int nd) {
  assert(nd <= MAX_SOBOL_DIM);
  uint32_t bits = HashToRndInt(idx * nd + dim, seed);
  int msb = GetMSB(idx);
  while (idx > 0) {
    int next_idx = idx ^ (1 << msb);
    next_idx ^= sobol_xors[dim][msb];
    int next_msb = GetMSB(next_idx);

    // The main key here is that we use some number of bits from the *next*
    // value, depending on the gap between this msb and the next msb.
    uint32_t rand_bits = HashToRndInt(next_idx * nd + dim, seed);
    int bits_to_set = (msb - next_msb);
    // The ^1 here corresponds to stratum swapping.
    rand_bits = (rand_bits >> (32 - bits_to_set)) ^ 1;
    bits = (rand_bits << (32 - bits_to_set)) ^ (bits >> bits_to_set);

    msb = next_msb;
    idx = next_idx;
  }

  return static_cast<double>(bits) / 4298115584.0;
}

}  // namespace sampling
