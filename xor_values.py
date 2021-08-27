#!/usr/bin/env python

"""Generates the Faure and pmj02 xor_values, used in xor_values.h."""

import numpy as np
import scipy.linalg

def print_hex(vals):
  array_elements = ", ".join([hex(val) for val in vals])
  string = "{{{0}}}".format(array_elements)
  print(string)

def get_stirling_matrix(size, base=2):
  m = np.identity(size).astype(int)
  for k in range(size):
    for n in range(k, size):
      if (n == 0) and (k == 0): m[k][n] = 1
      elif (k == 0): m[k][n] = ((n - 1) * m[k][n-1]) % base
      else:
        m[k][n] = (m[k-1][n-1] + (n - 1) * m[k][n-1]) % base
  return m

def get_xor_values(gen_matrix, base=2):
  # Invert the matrix.
  m_inv = np.linalg.inv(gen_matrix).astype(int) % base
  # Truncate the diagonal.
  m_inv -= np.identity(m_inv.shape[0], dtype=int)

  # Compute the xor-values from the negated columns.
  # For base 2, base_pow = (1,2,4,8,16...)
  base_pow = np.power(base, np.arange(0, gen_matrix.shape[0]))
  return [np.sum((-col % base)*base_pow) for col in m_inv.T]

def print_pmj02_xors(n):
  s = get_stirling_matrix(n+1)
  xor_vals = get_xor_values(s[:n, :n], base=2)
  print_hex(xor_vals)
  xor_vals = get_xor_values(s[1:n+1, 1:n+1], base=2)
  print_hex(xor_vals)

def print_faure_xors(base, n):
  print("Base {0} Faure (0,{0}) xor-values: ".format(base))
  m = np.identity(n)
  pascal_matrix = scipy.linalg.pascal(n, kind='upper') % base

  for i in range(base):
    xor_vals = get_xor_values(m, base)
    print_hex(xor_vals) 
    m = np.matmul(m, pascal_matrix) % base

print_pmj02_xors(32)
print_faure_xors(2, 32)
print_faure_xors(3, 21)
print_faure_xors(5, 14)
print_faure_xors(7, 12)
print_faure_xors(11, 10)