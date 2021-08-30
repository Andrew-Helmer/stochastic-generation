/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * Defines a function to generate stochastic/scrambled Halton sequences.
 */
#ifndef SAMPLING_SHALTON_H
#define SAMPLING_SHALTON_H

#include <cstdint>

namespace sampling {

// Fills the nd-dimensional sample array with a stochastically generated
// Halton sequence.
void GetStochasticHaltonSamples(const int num_samples,
                                const int nd,
                                const bool shuffle,
                                const int candidates,
                                const bool owen,
                                double *samples);

}  // namespace sampling

#endif  // SAMPLING_SHALTON_H
