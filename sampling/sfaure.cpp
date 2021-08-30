/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * Implementations of stochastic/scrambled Faure (0,s)-sequences.
 */
#include "sfaure.h"

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

using std::vector;

// Computes the strata offsets from a previous pass to this pass. This is used
// for correlated swapping.
template <int b>
vector<int> GetStrataOffsets(const int pass, const int nd,
                             const uint32_t baseSeed) {
  vector<int> strata_offsets(nd);
  for (int d = 0; d < nd; d++) {
    const uint32_t seed = CombineHashes(baseSeed, Hash(d));
    int strataOffset = (Permute(pass, b - 1, seed) + 1);
    // We actually want the offsets from a previous pass, not the previous
    // power, so we subtract the offset from the previous pass.
    if (pass > 0) strataOffset -= (Permute(pass - 1, b - 1, seed) + 1);
    strata_offsets[d] = (b + strataOffset) % b;
  }
  return strata_offsets;
}

// For a given sample, calculates the strata that sample should be placed in,
// one for each dimension, using correlated swapping.
template <int b, int n_xor_vals, const uint32_t xor_vals[b][n_xor_vals]>
inline void GetFaureStrataCS(const int prev_pass_idx,
                             const int nd,
                             const int log_n,
                             const int num_strata,
                             const vector<int>& strata_offsets,
                             const double* samples,
                             vector<int>* strata) {
  for (int d = 0; d < nd; d++) {
    const int prev_idx =
        CarrylessAdd<b>(prev_pass_idx, xor_vals[d][log_n]);
    const int prev_stratum = samples[prev_idx * nd + d] * num_strata;
    (*strata)[d] = AddLSDigit<b>(prev_stratum, strata_offsets[d]);
  }
}

// Generates a sample point in a given set of 1D strata.
inline void GetFaurePoint(const vector<int>& strata,
                          const int num_strata,
                          const int nd,
                          RNG* rng,
                          double* sample) {
  for (int d = 0; d < nd; d++) {
    // Xor-values for the Halton sequence are always zero, so we just need to
    // subtract to get the index in previous pass.
    sample[d] = (rng->GetUniformFloat() + strata[d]) / num_strata;
  }
}

// Given a set of 1D strata, one in each dimension, generates a number of
// candidate points in those strata, finds the candidate with the highest
// minimum 2D distance using the sample_grid, and then writes that point into
// the "sample" value.
inline void GetBestFaurePoint(const vector<int>& strata,
                              const int num_strata,
                              const SampleGrid2D& sample_grid,
                              const int n_candidates,
                              const int nd,
                              RNG* rng,
                              double* sample) {
  // Generate a single n-dimensional point in the output.
  GetFaurePoint(strata, num_strata, nd, rng, sample);

  if (n_candidates > 1) {
    // We actually only generate 2D candidates.
    std::vector<SampleGrid2D::Point2D> candidates(n_candidates);
    candidates[0] = {sample[0], sample[1]};
    for (int cand_idx = 1; cand_idx < n_candidates; cand_idx++) {
      double candidate[2];
      GetFaurePoint(strata, num_strata, /*nd=*/2, rng, candidate);
      candidates[cand_idx] = {candidate[0], candidate[1]};
    }
    int best_candidate = sample_grid.GetBestCandidate(candidates);
    // Write the best two dimensions into the output sample.
    sample[0] = candidates[best_candidate].x;
    sample[1] = candidates[best_candidate].y;
  }
}

// Generalized implementation of stochastically generated Faure sequences,
// with correlated swapping.
template <int b, int n_xor_vals, const uint32_t xor_vals[b][n_xor_vals]>
void GetSFaureCSSamples(const int num_samples,
                        const int nd,
                        const int n_candidates,
                        double* samples) {
  ASSERT(nd <= b, "Only " << b << " dimensions allowed.");
  RNG rng;
  for (int d = 0; d < nd; d++) samples[d] = rng.GetUniformFloat();

  SampleGrid2D sample_grid;
  if (n_candidates > 1) sample_grid.AddSample({samples[0], samples[1]});

  int next_sample = 1;
  int num_1d_strata = 1;

  vector<int> strata(nd);
  for (int log_n = 0; next_sample < num_samples; log_n++) {
    num_1d_strata *= b;

    const uint32_t strata_offset_seed = rng.GetUniformInt();
    const int prevLen = next_sample;
    for (int pass = 0; pass < b - 1; pass++) {
      const vector<int> strata_offsets =
          GetStrataOffsets<b>(pass, nd, strata_offset_seed);
      for (int i = 0; i < prevLen; i++) {
        const int prev_pass_idx = i + (prevLen * pass);

        GetFaureStrataCS<b, n_xor_vals, xor_vals>(prev_pass_idx, nd, log_n,
                                                  num_1d_strata, strata_offsets,
                                                  samples, &strata);
        double* sample = &(samples[next_sample * nd]);
        GetBestFaurePoint(
            strata, num_1d_strata, sample_grid, n_candidates, nd, &rng, sample);

        if (n_candidates > 1) sample_grid.AddSample({sample[0], sample[1]});

        next_sample++;
        if (next_sample >= num_samples) break;
      }
      if (next_sample >= num_samples) break;
    }
  }
}

// The Owen-scrambled version of the stochastic Faure sequence is the same
// except that we use independently chosen strata offsets for each new sample.
template <int b>
int GetStrataOffset(const uint32_t base_seed, const int sample_interval,
                    const int pass) {
  uint32_t seed = CombineHashes(base_seed, Hash(sample_interval));
  int strata_offset = (Permute(pass, b - 1, seed) + 1);
  if (pass > 0) strata_offset -= (Permute(pass - 1, b - 1, seed) + 1);
  return (b + strata_offset) % b;
}

template <int b, int n_xor_vals, const uint32_t xor_vals[b][n_xor_vals]>
inline void GetFaureStrataOwen(const int prev_pass_idx,
                               const int nd,
                               const int log_n,
                               const int pass,
                               const int num_strata,
                               const vector<uint32_t>& strata_offset_seeds,
                               const double* samples,
                               vector<int>* strata) {
  for (int d = 0; d < nd; d++) {
    const int prev_idx =
        CarrylessAdd<b>(prev_pass_idx, xor_vals[d][log_n]);
    const int prev_stratum = samples[prev_idx * nd + d] * num_strata;

    // Get the strata offset for this interval and pass.
    const int strata_interval = prev_stratum / b;
    const int strata_offset =
        GetStrataOffset<b>(strata_offset_seeds[d], strata_interval, pass);

    (*strata)[d] = AddLSDigit<b>(prev_stratum, strata_offset);
  }
}

template <int b, int n_xor_vals, const uint32_t xor_vals[b][n_xor_vals]>
void GetSFaureOwenSamples(int num_samples, int nd, int n_candidates,
                          double* samples) {
  ASSERT(nd <= b, "Only " << b << " dimensions allowed.");
  RNG rng;
  for (int d = 0; d < nd; d++) samples[d] = rng.GetUniformFloat();

  SampleGrid2D sample_grid;
  if (n_candidates > 1) sample_grid.AddSample({samples[0], samples[1]});

  int next_sample = 1;
  int num_1d_strata = 1;

  vector<uint32_t> strata_offset_seeds(nd);
  vector<int> strata(nd);
  for (int log_n = 0; next_sample < num_samples; log_n++) {
    num_1d_strata *= b;

    // Get the base seeds for each dimension.
    for (int d = 0; d < nd; d++) strata_offset_seeds[d] = rng.GetUniformInt();

    const int prevLen = next_sample;
    for (int pass = 0; pass < b - 1; pass++) {
      for (int i = 0; i < prevLen; i++) {
        const int prev_pass_idx = i + (prevLen * pass);

        GetFaureStrataOwen<b, n_xor_vals, xor_vals>(
            prev_pass_idx, nd, log_n, pass, num_1d_strata, strata_offset_seeds,
            samples, &strata);

        double* sample = &(samples[next_sample * nd]);
        GetBestFaurePoint(
            strata, num_1d_strata, sample_grid, n_candidates, nd, &rng, sample);

        if (n_candidates > 1) sample_grid.AddSample({sample[0], sample[1]});

        next_sample++;
        if (next_sample >= num_samples) break;
      }
      if (next_sample >= num_samples) break;
    }
  }
}

template <int b, int n_xor_vals, const uint32_t xor_vals[b][n_xor_vals]>
void GetSFaureSamples(int num_samples, int nd, bool shuffle,
                      int candidates, bool owen, double* samples) {
  if (!owen) {
    GetSFaureCSSamples<b, n_xor_vals, xor_vals>(
        num_samples, nd, candidates, samples);
  } else {
    GetSFaureOwenSamples<b, n_xor_vals, xor_vals>(
        num_samples, nd, candidates, samples);
  }

  if (shuffle) ProgressiveShuffleSamples<b>(num_samples, nd, samples);
}
}  // namespace

/*
 * These functions are basically just wrappers around SFaureSamples<>.
 */
void GetStochasticFaure03Samples(int num_samples, int nd, bool shuffle,
                                 int candidates, bool owen, double* samples) {
  GetSFaureSamples<3, FAURE03_XOR_SIZE, faure03_xors>(
        num_samples, nd, shuffle, candidates, owen, samples);
}

void GetStochasticFaure05Samples(int num_samples, int nd, bool shuffle,
                                 int candidates, bool owen, double* samples) {
  GetSFaureSamples<5, FAURE05_XOR_SIZE, faure05_xors>(
        num_samples, nd, shuffle, candidates, owen, samples);
}

void GetStochasticFaure07Samples(int num_samples, int nd, bool shuffle,
                                 int candidates, bool owen, double* samples) {
  GetSFaureSamples<7, FAURE07_XOR_SIZE, faure07_xors>(
        num_samples, nd, shuffle, candidates, owen, samples);
}

void GetStochasticFaure011Samples(int num_samples, int nd, bool shuffle,
                                  int candidates, bool owen, double* samples) {
  GetSFaureSamples<11, FAURE011_XOR_SIZE, faure011_xors>(
        num_samples, nd, shuffle, candidates, owen, samples);
}
}  // namespace sampling
