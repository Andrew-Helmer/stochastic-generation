#!/usr/bin/env python

"""Generates the Faure, pmj02, and SZ sequence xor_values, used in xor_values.h.

This works by constructing the generator matrices, and then "reversing" those
matrices as explained in Section 4.4 of the paper.
"""

from math import comb
import numpy as np
import scipy.linalg

def get_xor_values(gen_matrix, base=2):
  """Compute xor-values from a non-singular upper triangular matrix.

  Listing 3 from the paper.
  """
  # Invert the matrix.
  m_inv = np.linalg.inv(gen_matrix).astype(int) % base
  # Truncate the diagonal.
  m_inv -= np.identity(m_inv.shape[0], dtype=int)

  # Compute the xor-values from the negated columns.
  # For base 2, base_pow = (1,2,4,8,16...)
  base_pow = np.power(base, np.arange(0, gen_matrix.shape[0]))
  return [np.sum((-col % base)*base_pow) for col in m_inv.T]

def print_hex(vals):
  """For an array of numbers, print them as hexadecimal values."""
  array_elements = ", ".join([hex(val) for val in vals])
  string = "{{{0}}}".format(array_elements)
  print(string)

def get_stirling_matrix(size, base=2):
  """Construct a matrix of unsigned Stirling numbers of the first kind."""
  m = np.identity(size).astype(int)
  for k in range(size):
    for n in range(k, size):
      if (n == 0) and (k == 0): m[k][n] = 1
      elif (k == 0): m[k][n] = ((n - 1) * m[k][n-1]) % base
      else:
        m[k][n] = (m[k-1][n-1] + (n - 1) * m[k][n-1]) % base
  return m

def print_pmj02_xors(n):
  """Print pmj02 xor-values using Stirling numbers.

  See section 1.1 of the supplemental materials.
  """
  s = get_stirling_matrix(n+1)
  xor_vals = get_xor_values(s[:n, :n], base=2)
  print_hex(xor_vals)
  # The second matrix is shifted up and to the left by one.
  xor_vals = get_xor_values(s[1:n+1, 1:n+1], base=2)
  print_hex(xor_vals)

def print_faure_xors(base, n):
  """Print Faure sequences xor-values"""
  print("Base {0} Faure (0,{0}) xor-values: ".format(base))
  m = np.identity(n)
  # Faure sequences are constructed from exponentiated Pascal matrices.
  pascal_matrix = scipy.linalg.pascal(n, kind='upper') % base

  for i in range(base):
    xor_vals = get_xor_values(m, base)
    print_hex(xor_vals) 
    m = np.matmul(m, pascal_matrix) % base

def repeat_diagonal(block, target_size):
  """Repeat a square block matrix along the main diagonal to target_size."""
  block_size = len(block)
  num_repeats = target_size // block_size
  zero_block = np.zeros((block_size, block_size), dtype=int)
  return np.block([
      [block if row == col else zero_block for col in range(num_repeats)]
      for row in range(num_repeats)
  ])

def get_block_pascal_matrix(block, matrix_size):
  """Construct an upper-triangular block-Pascal matrix from a square block."""
  block_size = len(block)
  num_blocks = matrix_size // block_size
  block_powers = [np.identity(block_size, dtype=int)]
  for _ in range(num_blocks - 1):
    block_powers.append(np.matmul(block_powers[-1], block) % 2)

  pascal_matrix = np.zeros((matrix_size, matrix_size), dtype=int)
  for row in range(num_blocks):
    for col in range(row, num_blocks):
      if comb(col, row) % 2 == 1:
        row_start, row_end = row * block_size, (row + 1) * block_size
        col_start, col_end = col * block_size, (col + 1) * block_size
        pascal_matrix[row_start:row_end, col_start:col_end] = (
            block_powers[col - row]
        )
  return pascal_matrix

def generate_sz_basis_matrices(num_matrices):
  """Generate the basis block matrices used to construct SZ sequence matrices.

  Each doubling step constructs a quadratic finite-field extension over GF(2)
  in matrix form, using the companion matrix [[0, I], [beta, I]] where beta is
  the last basis matrix from the previous extension level.
  """
  basis_matrices = [np.array([[1]], dtype=int)]
  subblock_size = 1
  while subblock_size < num_matrices:
    beta = basis_matrices[subblock_size - 1]
    zero_subblock = np.zeros((subblock_size, subblock_size), dtype=int)
    for idx in range(subblock_size):
      if subblock_size + idx >= num_matrices:
        break
      enlarged_basis = repeat_diagonal(basis_matrices[idx], subblock_size)
      next_basis = np.block([
          [zero_subblock,                           enlarged_basis],
          [np.matmul(beta, enlarged_basis) % 2,     enlarged_basis]
      ])
      basis_matrices.append(next_basis)
    subblock_size *= 2
  return basis_matrices

def generate_sz_matrices(num_dims, matrix_size=32):
  """Generate non-singular upper-triangular SZ generator matrices.

  Uses a novel direct construction for SZ sequences ("SZ Sequences:
  Binary-Based (0, 2^q)-Sequences" by Ahmed et al., 2025) that builds each
  dimension's generator matrix directly from block-Pascal matrices of basis
  blocks without any brute-force search or sorting.
  """
  matrices = [np.identity(matrix_size, dtype=int)]
  if num_dims <= 1:
    return matrices[:num_dims]

  num_basis = (num_dims - 1).bit_length()
  basis_matrices = generate_sz_basis_matrices(num_basis)
  pascal_matrices = [
      get_block_pascal_matrix(basis, matrix_size) for basis in basis_matrices
  ]

  # Each block-Pascal matrix doubles the number of generated dimensions by
  # right-multiplying the previously generated matrices.
  for pascal_matrix in pascal_matrices:
    num_remaining = num_dims - len(matrices)
    matrices += [
        np.matmul(prev_matrix, pascal_matrix) % 2
        for prev_matrix in matrices[:num_remaining]
    ]
  return matrices

def print_sz_xors(num_dims, matrix_size=32):
  """Print SZ sequence xor-values."""
  print("SZ sequence ({0} dimensions) xor-values: ".format(num_dims))
  for gen_matrix in generate_sz_matrices(num_dims, matrix_size):
    xor_vals = get_xor_values(gen_matrix, base=2)
    print_hex(xor_vals)

print_pmj02_xors(32)
print_faure_xors(2, 32)
print_faure_xors(3, 21)
print_faure_xors(5, 14)
print_faure_xors(7, 12)
print_faure_xors(11, 10)
print_sz_xors(16, 32)
