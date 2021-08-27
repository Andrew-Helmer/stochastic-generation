/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * This header defines functions to stochastically generate both the
 * high-dimensional Sobol' sequence, as well as the pmj02 sequence from
 * "Progressive Multi-Jittered Sample Sequences"
 * by Christensen, Kensler, and Kilpatrick (2018)
 */
#ifndef SAMPLING_SSOBOL_H
#define SAMPLING_SSOBOL_H

#include <cstdint>

namespace sampling {

// Fills the nd-dimensional sample array with an Owen-scrambled
// Sobol sequence, using the stochastic generation algorithm. Each of the nd
// coordinates are stored for a sample are stored contiguously (i.e. in
// "sample major" order).
void GetStochasticSobolSamples(int nSamples, int nd, bool shuffle,
                               int candidates, bool owen, double* samples);

// Fills the nd-dimensional sample array with a PMJ02 sequence, using the
// stochastic generation algorithm and the xor-values derived in the
// supplemental materials. ND must be 1 or 2, for higher dimensions use the
// ssobol sequence.
void GetPMJ02Samples(int nSamples, int nd, bool shuffle, int candidates,
                     bool owen, double* samples);

// Query for an arbitrary coordinate from a scrambled Sobol (0,2)-sequence. Seed
// can be the same across all dimensions, or it can be different for each
// dimension.
//
// nd is the number of dimensions shared by the same seed value, to get separate
// hashed values for each sample index. If a different seed is provided for each
// dimension, this can be used with nd=1.
double GetSobolStateless(int idx, int dim, uint32_t seed, int nd = 2);
// Same as above but implemented iteratively.
double GetSobolStatelessIter(int idx, int dim, uint32_t seed, int nd = 2);

}  // namespace sampling

#endif  // SAMPLING_SSOBOL_H
