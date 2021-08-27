/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * This headers defines utilities for performing best-candidate sampling in
 * to improve minimum distances between points.
 */
#ifndef SAMPLING_BN_UTILS_H
#define SAMPLING_BN_UTILS_H

#include <deque>
#include <vector>

namespace sampling {

static constexpr int kBestCandidateSamples = 100;

/*
 * SampleGrid2D is a simple utility class to keep track of points in a 2D
 * uniform grid, and query minimum distances within that grid. When points are
 * added to the grid, if there is another point in the same grid cell, the grid
 * is subdivided until there are an equal number of points in each cell.
 *
 * We use this rather than a tree-structure like a kd-tree or quadtree, because
 * the sample sequences we're working on all have 2D stratification guarantees,
 * so the sample grid will never get too much bigger than the number of samples
 * (although it can be not great with a higher-base sequence,
 * like the SFaure011).
 */
class SampleGrid2D {
 public:
  SampleGrid2D() : sample_grid_(1, nullptr) {}

  struct Point2D {
    double x;
    double y;
  };

  // Add a new sample point into the grid. This may automatically subdivide the
  // grid, if the new point is in the same cell as a previous point.
  void AddSample(const Point2D sample);

  // From a list of candidates, return the index of the candidate with the
  // furthest minimum toroidal distance from points in the grid.
  int GetBestCandidate(const std::vector<Point2D>& candidates) const;

 private:
  double GetMinDistSqFast(const Point2D sample,
                          const double max_min_dist_sq,
                          Point2D* nearest) const;
  void SubdivideGrid(const int subdivisions);

  // The actual grid is just stored as a vector of pointers.
  std::vector<Point2D*> sample_grid_;
  // We keep track of all our points individually. We use a deque to make sure
  // that pointers are persistent. A memory arena would probably be better, but
  // this is simple.
  std::deque<Point2D> points_;
  int width_ = 1;
};

}  // namespace sampling

#endif  // SAMPLING_BN_UTILS_H
