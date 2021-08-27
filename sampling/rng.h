/*
 * Licensed under MIT Open-Source License: see LICENSE.
 * Note from Andrew Helmer: virtually everything in here was taken from another
 * source, either a paper or something with a permissive license, so I do not
 * claim copyright here.
 *
 * Implementations of pseudorandom number generations and hash functions.
 * These are implemented in the header file, because it makes the sampling code
 * much faster.
 */
#ifndef SAMPLING_RNG_H
#define SAMPLING_RNG_H

#include <cstdint>
#include <limits>
#include <random>

namespace sampling {

/*
 * C++ standard random number generation. Can be used to seed the LCD and PCG
 * classes below.
 */
static uint64_t GetDeviceUniformInt() {
  thread_local static std::random_device r;
  thread_local static std::default_random_engine gen(r());
  thread_local static std::uniform_int_distribution<uint64_t> uniform;

  static const std::uniform_int_distribution<uint64_t>::param_type param(
      0, std::numeric_limits<uint64_t>::max());

  return uniform(gen, param);
}


#define USE_PCG32 1
/*
 * Random number generator. Uses either PCG32 (Permuted Congruential Generator)
 * or a linear congruential generator, equivalent to drand48(), depending on
 * whether USE_PCG32 is set to 1.
 */
class RNG {
#if USE_PCG32
// PCG implementation based on the supplemental code provided in
//  * "Practical Hash-based Owen Scrambling" (Burley 2020), which is in turn
// from pcg-random.org.
 public:
  RNG() {
    pcg32_srandom_r(GetDeviceUniformInt(), 0);
  }
  explicit RNG(uint64_t seed) {
    pcg32_srandom_r(seed, 0);
  }
  uint32_t GetUniformInt() {
    uint64_t oldstate = state;
    // Advance internal state
    state = oldstate * 6364136223846793005ULL + (inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
  }
  double GetUniformFloat() {
    static constexpr double S = static_cast<double>(1.0/(1ul<<32));
    return GetUniformInt() * S;
  }

 private:
  void pcg32_srandom_r(uint64_t initstate,
                       uint64_t initseq) {
    state = 0U;
    inc = (initseq << 1u) | 1u;
    GetUniformInt();
    state += initstate;
    GetUniformInt();
  }

  uint64_t state;
  uint64_t inc;
#else 
// Linear congruential generator provided by Andrew Kensler. Should give
// results equivalent to drand48(), but it's faster to have it in a header file.
 public:
  RNG() {drand48_rng = GetDeviceUniformInt();}
  explicit RNG(uint64_t seed) : drand48_rng(seed) {}
  uint64_t GetUniformInt() {
    drand48_rng = 0x5deece66d * drand48_rng + 0xb;
    return drand48_rng >> 17 & 0x7fffffff;
  }
  double GetUniformFloat() {
    return static_cast<double>(GetUniformInt()) / 2147483648.0;
  }

 private:
  uint64_t drand48_rng;
#endif
 public:
  // Get an integer in the range [min, max), without bias.
  uint32_t GetUniformInt(uint32_t min, uint32_t max) {
    static constexpr uint32_t int_max = std::numeric_limits<uint32_t>::max();
    const uint32_t range = max-min;
    uint32_t base_int;
    // Reject a uniform int that's too high, otherwise the values are biased.
    const int max_value = int_max - (int_max % range);
    do {
      base_int = this->GetUniformInt();
    } while (base_int > max_value);
    return (base_int % range) + min;
  }
};

/*
 * permute() from Correlated Multi-Jittered Sampling (Kensler 2013).
 * For a given seed,
 * this will map every index to a different, unique, permuted index in [0, len].
 */
inline uint32_t Permute(uint32_t idx, uint32_t len, uint32_t seed) {
  uint32_t mask = len-1;
  mask |= mask >> 1;
  mask |= mask >> 2;
  mask |= mask >> 4;
  mask |= mask >> 8;
  mask |= mask >> 16;

  do {
    idx ^= seed; idx *= 0xe170893d;
    idx ^= seed >> 16;
    idx ^= (idx & mask) >> 4;
    idx ^= seed >> 8; idx *= 0x0929eb3f;
    idx ^= seed >> 23;
    idx ^= (idx & mask) >> 1; idx *= 1 | seed >> 27;
    idx *= 0x6935fa69;
    idx ^= (idx & mask) >> 11; idx *= 0x74dcb303;
    idx ^= (idx & mask) >> 2; idx *= 0x9e501cc3;
    idx ^= (idx & mask) >> 2; idx *= 0xc860a3df;
    idx &= mask;
    idx ^= idx >> 5;
  } while (idx >= len);
  return (idx + seed) % len;
}

/*
 * Hash functions from supplemental code of
 * "Practical Hash-based Owen Scrambling" (Burley 2020).
 */
inline uint32_t CombineHashes(uint32_t seed, uint32_t v) {
  return seed ^ (v + (seed << 6) + (seed >> 2));
}

inline uint32_t Hash(uint32_t x) {
    // finalizer from murmurhash3
    x ^= x >> 16;
    x *= 0x85ebca6bu;
    x ^= x >> 13;
    x *= 0xc2b2ae35u;
    x ^= x >> 16;
    return x;
}

}  // namespace sampling

#endif  // SAMPLING_RNG_H
