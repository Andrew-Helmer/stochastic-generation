/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * Just a few utilities used in various places in the code.
 */
#ifndef SAMPLING_UTILS_H
#define SAMPLING_UTILS_H

#include <cstdint>
#include <iostream>

#   define ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#   define ERRIF(condition, message) \
    do { \
        if ((condition)) { \
            std::cerr << message << std::endl; \
        } \
    } while (false)

namespace sampling {
/*
 * Generalization of bitwise xor to other prime bases. Adds the BASE digits of i
 * and j modulo the base.
 */
template<unsigned BASE>
unsigned CarrylessAdd(unsigned i, unsigned j) {
  unsigned sum = 0, bPow = 1;
  while (j > 0 && i > 0) {
    sum += ((i + j) % BASE) * bPow;
    i /= BASE; j /= BASE;
    bPow *= BASE;
  }
  return sum + bPow*(i+j);
}

/*
 * Adds a single digit in an arbitrary base to the least significant digit,
 * without carrying. Specialized (faster) case of the above function, used for
 * strata swapping. j should be a single digit, i.e. 0 <= j < BASE.
 */
template<unsigned BASE>
unsigned AddLSDigit(unsigned i, unsigned j) {
  int i_lsb = i % BASE;
  i = i - i_lsb;
  int lsb = i_lsb+j;
  return i + (lsb >= BASE ? lsb-BASE : lsb);
}

}  // namespace sampling

#endif  // SAMPLING_UTILS_H
