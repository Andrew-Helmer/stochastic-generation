/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * Implementations of stochastic/scrambled Halton sequences.
 */
#include "shalton.h"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>

#include "bn_utils.h"
#include "rng.h"
#include "utils.h"

namespace sampling {

using std::vector;

namespace {
/*
 * We statically define the prime numbers and the functions to add least
 * significant digits in those bases, which make them faster.
 */
static constexpr int MAX_HALTON_DIM = 10;
static constexpr uint32_t primes[MAX_HALTON_DIM] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29
};
typedef unsigned (*add_ls_fn) (unsigned, unsigned);
static constexpr add_ls_fn add_ls_digit_fns[MAX_HALTON_DIM] = {
  &AddLSDigit<primes[0]>,
  &AddLSDigit<primes[1]>,
  &AddLSDigit<primes[2]>,
  &AddLSDigit<primes[3]>,
  &AddLSDigit<primes[4]>,
  &AddLSDigit<primes[5]>,
  &AddLSDigit<primes[6]>,
  &AddLSDigit<primes[7]>,
  &AddLSDigit<primes[8]>,
  &AddLSDigit<primes[9]>,
};

// Computes a stratum offset from a previous pass to this pass. This is used
// for correlated swapping.
int GetStratumOffset(const int pass,
                     const int b,
                     const uint32_t seed) {
  int stratum_offset = (Permute(pass, b - 1, seed) + 1);
  // We actually want the offsets from a previous pass, not the previous
  // power, so we subtract the offset from the previous pass.
  if (pass > 0) stratum_offset -= (Permute(pass - 1, b - 1, seed) + 1);
  return (b + stratum_offset) % b;
}

// With correlated swapping, we calculate all the 1D strata for a new point.
inline void GetHaltonStrataCS(const int i,
                              const int nd,
                              const vector<int>& num_strata,
                              const vector<int>& strata_offsets,
                              const double* samples,
                              vector<int>* strata) {
  for (int d = 0; d < nd; d++) {
    const int base = primes[d];
    // The beautiful thing about the Halton sequence is that, although all the
    // dimensions are in different bases, the "xor_values" are all zero. This
    // means we can find the previous sample index using only a subtraction.
    const int prev_idx = i - num_strata[d]/base;
    const int prev_stratum = samples[prev_idx * nd + d] * num_strata[d];
    (*strata)[d] = add_ls_digit_fns[d](prev_stratum, strata_offsets[d]);
  }
}

// Given a set of strata, one for each dimension, generate a new point.
inline void GetHaltonPoint(const vector<int>& strata,
                           const vector<int>& num_strata,
                           const int nd,
                           RNG* rng,
                           double* sample) {
  for (int d = 0; d < nd; d++) {
    // Xor-values for the Halton sequence are always zero, so we just need to
    // subtract to get the index in previous pass.
    sample[d] = (rng->GetUniformFloat() + strata[d]) / num_strata[d];
  }
}

// Given a set of 1D strata, one in each dimension, generates a number of
// candidate points in those strata, finds the candidate with the highest
// minimum 2D distance using the sample_grid, and then writes that point into
// the "sample" value.
inline void GetBestHaltonPoint(const vector<int>& strata,
                               const vector<int>& num_strata,
                               const SampleGrid2D& sample_grid,
                               const int n_candidates,
                               const int nd,
                               RNG* rng,
                               double* sample) {
  // Generate a single n-dimensional point in the output.
  GetHaltonPoint(strata, num_strata, nd, rng, sample);

  if (n_candidates > 1) {
    // We actually only generate 2D candidates.
    std::vector<SampleGrid2D::Point2D> candidates(n_candidates);
    candidates[0] = {sample[0], sample[1]};
    for (int cand_idx = 1; cand_idx < n_candidates; cand_idx++) {
      double candidate[2];
      GetHaltonPoint(strata, num_strata, /*nd=*/2, rng, candidate);
      candidates[cand_idx] = {candidate[0], candidate[1]};
    }
    int best_candidate = sample_grid.GetBestCandidate(candidates);
    // Write the best two dimensions into the output sample.
    sample[0] = candidates[best_candidate].x;
    sample[1] = candidates[best_candidate].y;
  }
}

void GetStochasticHaltonCSSamples(const int num_samples,
                                  const int nd,
                                  const int n_candidates,
                                  double *samples) {
  ASSERT(nd <= MAX_HALTON_DIM,
         "Only " << MAX_HALTON_DIM << " dimensions allowed.");
  RNG rng;
  for (int d = 0; d < nd; d++)
    samples[d] = rng.GetUniformFloat();

  SampleGrid2D sample_grid;
  if (n_candidates > 1) sample_grid.AddSample({samples[0], samples[1]});

  // Because the Halton sequence has a different base for each dimension, we
  // need to keep track of a bunch of things separately for each dimension.
  vector<int> num_strata(nd, 1);
  vector<int> cur_pass(nd, 0);
  vector<uint32_t> strata_hash_seeds(nd);
  vector<int> strata_offsets(nd);
  vector<int> strata(nd);
  for (int i = 1; i < num_samples; i++) {
    for (int d = 0; d < nd; d++) {
      const int base = primes[d];
      if (i >= num_strata[d]) {
        // If all of the 1D strata are occupied, we subdivide the 1D strata,
        // reset the pass to zero.
        num_strata[d] *= base;
        cur_pass[d] = 0;
        strata_hash_seeds[d] = rng.GetUniformInt();
        strata_offsets[d] =
            GetStratumOffset(/*pass=*/0, base, strata_hash_seeds[d]);
      } else if (i*base >= num_strata[d]*(cur_pass[d]+2)) {
        // Check to see if we've reached the next "pass" for this dimension.
        cur_pass[d] += 1;
        strata_offsets[d] =
            GetStratumOffset(cur_pass[d], base, strata_hash_seeds[d]);
      }
    }

    GetHaltonStrataCS(
        i, nd, num_strata, strata_offsets, samples, &strata);

    double* sample = &(samples[i * nd]);
    GetBestHaltonPoint(strata, num_strata, sample_grid, n_candidates, nd, &rng,
                       sample);
    if (n_candidates > 1) sample_grid.AddSample({sample[0], sample[1]});
  }
}

int GetStratumOffset(const int pass,
                     const int sample_interval,
                     const int b,
                     const uint32_t base_seed) {
  uint32_t seed = CombineHashes(base_seed, Hash(sample_interval));
  return GetStratumOffset(pass, b, seed);
}

inline void GetHaltonStrataOwen(const int i,
                                const int nd,
                                const vector<int>& num_strata,
                                const vector<int>& cur_pass,
                                const vector<uint32_t>& strata_hash_seeds,
                                const double* samples,
                                vector<int>* strata) {
  for (int d = 0; d < nd; d++) {
    const int base = primes[d];
    // The beautiful thing about the Halton sequence is that, although all the
    // dimensions are in different bases, the "xor_values" are all zero. This
    // means we can find the previous sample index using only a subtraction.
    const int prev_idx = i - num_strata[d] / base;
    const int prev_stratum = samples[prev_idx * nd + d] * num_strata[d];
    const int interval = prev_stratum / base;
    // We independently calculate the strata offsets/deltas for each sample.
    const int stratum_offset =
        GetStratumOffset(cur_pass[d], interval, base, strata_hash_seeds[d]);
    (*strata)[d] = add_ls_digit_fns[d](prev_stratum, stratum_offset);
  }
}

void GetStochasticHaltonOwenSamples(const int num_samples,
                                    const int nd,
                                    const int n_candidates,
                                    double *samples) {
  ASSERT(nd <= MAX_HALTON_DIM,
         "Only " << MAX_HALTON_DIM << " dimensions allowed.");
  RNG rng;
  for (int d = 0; d < nd; d++)
    samples[d] = rng.GetUniformFloat();

  SampleGrid2D sample_grid;
  if (n_candidates > 1) sample_grid.AddSample({samples[0], samples[1]});

  // Because the Halton sequence has a different base for each dimension, we
  // need to keep track of a bunch of things separately for each dimension.
  vector<int> num_strata(nd, 1);
  vector<int> cur_pass(nd, 0);
  vector<uint32_t> strata_hash_seeds(nd);
  vector<int> strata(nd);
  for (int i = 1; i < num_samples; i++) {
    for (int d = 0; d < nd; d++) {
      const int base = primes[d];
      if (i >= num_strata[d]) {
        // If all of the 1D strata are occupied, we subdivide the 1D strata,
        // reset the pass to zero.
        num_strata[d] *= base;
        cur_pass[d] = 0;
        strata_hash_seeds[d] = rng.GetUniformInt();
      } else if (i*base >= num_strata[d]*(cur_pass[d]+2)) {
        // Check to see if we've reached the next "pass" for this dimension.
        cur_pass[d] += 1;
      }
    }

    GetHaltonStrataOwen(
        i, nd, num_strata, cur_pass, strata_hash_seeds, samples, &strata);

    double* sample = &(samples[i * nd]);
    GetBestHaltonPoint(strata, num_strata, sample_grid, n_candidates, nd, &rng,
                       sample);
    if (n_candidates > 1) sample_grid.AddSample({sample[0], sample[1]});
  }
}
}  // namespace

void GetStochasticHaltonSamples(const int num_samples,
                                const int nd,
                                const bool shuffle,
                                const int candidates,
                                const bool owen,
                                double *samples) {
  ASSERT(nd <= MAX_HALTON_DIM,
         "Only " << MAX_HALTON_DIM << " dimensions allowed.");
  ASSERT(!shuffle, "Cannot effectively shuffle Halton sequences.");
  if (owen) {
    GetStochasticHaltonOwenSamples(num_samples, nd, candidates, samples);
  } else {
    GetStochasticHaltonCSSamples(num_samples, nd, candidates, samples);
  }
}
}  // namespace sampling
