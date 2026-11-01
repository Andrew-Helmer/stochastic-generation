/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * This header defines functions to stochastically generate the
 * high-dimensional Sobol' sequence, the pmj02 sequence from
 * "Progressive Multi-Jittered Sample Sequences"
 * by Christensen, Kensler, and Kilpatrick (2018), and the SZ sequence from
 * "SZ Sequences: Binary-Based (0, 2^q)-Sequences" by Ahmed et al. (2025).
 */
#ifndef SAMPLING_SSOBOL_H
#define SAMPLING_SSOBOL_H

#include <cstdint>

namespace sampling {

// Fills the nd-dimensional sample array with an Owen-scrambled
// Sobol sequence, using the stochastic generation algorithm. Each of the nd
// coordinates are stored for a sample are stored contiguously (i.e. in
// "sample major" order).
void GetStochasticSobolSamples(const int num_samples,
                               const int nd,
                               const bool shuffle,
                               const int candidates,
                               const bool owen,
                               double* samples);

// Fills the nd-dimensional sample array with a PMJ02 sequence, using the
// stochastic generation algorithm and the xor-values derived in the
// supplemental materials. ND must be 2, for higher dimensions use the
// ssobol sequence.
void GetPMJ02Samples(const int num_samples,
                     const int nd,
                     const bool shuffle,
                     const int candidates,
                     const bool owen,
                     double* samples);

// Fills the nd-dimensional sample array with an Owen-scrambled SZ sequence,
// using the stochastic generation algorithm, with generator matrices from a
// new direct construction (see xor_values.py). Up to 16 dimensions.
//
// Every aligned pair of dimensions (2k, 2k+1) is a base-2 (0,2)-sequence,
// every aligned group of four (4k, ..., 4k+3) is a base-4 (0,4)-sequence, and
// all 16 dimensions together are a base-16 (0,16)-sequence. For rendering, the
// first four dimensions are an excellent choice, e.g. combining a 2D pixel
// position and a 2D lens position.
void GetStochasticSZSamples(const int num_samples,
                            const int nd,
                            const bool shuffle,
                            const int candidates,
                            const bool owen,
                            double* samples);

// Query for an arbitrary coordinate from a scrambled Sobol (0,2)-sequence. Seed
// can be the same across all dimensions, or it can be different for each
// dimension.
//
// nd is the number of dimensions shared by the same seed value, to get separate
// hashed values for each sample index. If a different seed is provided for each
// dimension, this can be used with nd=1.
double GetSobolStateless(const int idx,
                         const int dim,
                         const uint32_t seed,
                         const int nd = 2);
// Same as above but implemented iteratively.
double GetSobolStatelessIter(const int idx,
                             const int dim,
                             const uint32_t seed,
                             const int nd = 2);

// Pretty much just for testing the stateless Sobol functions, but similar
// to the GetStochasticSobolSamples above. Note that best candidate sampling
// doesn't work, candidates must be equal to 1.
void GetStochasticSobolStatelessSamples(const int num_samples,
                                        const int nd,
                                        const bool shuffle,
                                        const int candidates,
                                        const bool owen,
                                        double* samples);

}  // namespace sampling

#endif  // SAMPLING_SSOBOL_H
