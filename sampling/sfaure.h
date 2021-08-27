/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * Defines functions to generate stochastic/scrambled Faure (0,s)-sequences.
 * Can generate (0,3), (0,5), (0,7), and (0,11) sequences, with or without
 * decorrelation (shuffling), with or without best-candidate samping for the
 * first two dimensions.
 */
#ifndef SAMPLING_SFAURE_H
#define SAMPLING_SFAURE_H

namespace sampling {

// Fills the nd-dimensional sample array with a stochastically generated
// Faure sequences. These are all (0,s)-sequences, and the chosen dimensionality
// must be <= to the base, which is value after the zero. I.e. the (0,5)
// can only generate up to 5 dimensions.
//
// If owen=false, then these sequences use Correlated Shuffling from the paper,
// otherwise they use "un"correlated shuffling, i.e. the sequences are
// equivalent to fully Owen-scrambled sequences. We haven't observed a
// noticeable difference in the quality of either sequence.
void GetStochasticFaure03Samples(int nSamples, int nd, bool shuffle,
                                 int candidates, bool owen, double *samples);
void GetStochasticFaure05Samples(int nSamples, int nd, bool shuffle,
                                 int candidates, bool owen, double *samples);
void GetStochasticFaure07Samples(int nSamples, int nd, bool shuffle,
                                 int candidates, bool owen, double *samples);
void GetStochasticFaure011Samples(int nSamples, int nd, bool shuffle,
                                  int candidates, bool owen, double *samples);

}  // namespace sampling

#endif  // SAMPLING_SFAURE_H
