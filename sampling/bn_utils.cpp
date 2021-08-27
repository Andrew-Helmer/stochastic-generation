/*
 * Copyright (C) Andrew Helmer 2021.
 * Licensed under MIT Open-Source License: see LICENSE.
 *
 * This file implements the SampleGrid2D class, which is used to store 2D points
 * in a uniform grid to accelerate nearest-neighbor searches.
 */
#include "bn_utils.h"

#include <algorithm>
#include <cmath>
#include <utility>

#include "utils.h"

namespace sampling {
namespace {
inline int GetGridIndex(const SampleGrid2D::Point2D& point,
                        const int grid_width) {
  const int x_pos = point.x * grid_width;
  const int y_pos = point.y * grid_width;
  return y_pos * grid_width + x_pos;
}
}

void SampleGrid2D::AddSample(const Point2D sample) {
  int grid_idx = GetGridIndex(sample, width_);
  // If there's already a point in this grid cell, we need to repeatedly
  // subdivide the grid until there that point and the new point won't be in the
  // same grid cell.
  if (sample_grid_[grid_idx] != nullptr) {
    // Figure out how many times we need to subdivide the grid.
    Point2D conflicting_point = *sample_grid_[grid_idx];
    int subdivisions = 0;
    int temp_width = width_;
    do {
      subdivisions++;
      temp_width <<= 1;
    } while (static_cast<int>(conflicting_point.x * temp_width)
             == static_cast<int>(sample.x * temp_width) &&
             static_cast<int>(conflicting_point.y * temp_width)
             == static_cast<int>(sample.y * temp_width));

    // Subdivide the grid and recalculate the grid index for our new point.
    SubdivideGrid(subdivisions);
    grid_idx = GetGridIndex(sample, width_);
  }
  points_.push_back(sample);
  sample_grid_[grid_idx] = &(points_.back());
}

void SampleGrid2D::SubdivideGrid(const int subdivisions) {
  width_ <<= subdivisions;
  // We don't use resize because there's no reason to do the copy. Instead we
  // just allocate a new vector, which is fine because it's still a geometric
  // growth.
  sample_grid_ = std::vector<Point2D*>(width_*width_, nullptr);
  const int num_points = points_.size();
  // Add the points back into the grid.
  for (int i = 0; i < num_points; i++) {
    int grid_idx = GetGridIndex(points_[i], width_);
    sample_grid_[grid_idx] = &points_[i];
  }
}

namespace {
// Get the distance between two points, with wrap-around.
inline double GetToroidalDistanceSq(const SampleGrid2D::Point2D& p0,
                                    const SampleGrid2D::Point2D& p1) {
  double x_diff = abs(p1.x-p0.x);
  if (x_diff > 0.5) x_diff = 1.0 - x_diff;
  double y_diff = abs(p1.y-p0.y);
  if (y_diff > 0.5) y_diff = 1.0 - y_diff;

  return (x_diff*x_diff)+(y_diff*y_diff);
}
// Wrap an integer around a fixed length (i.e. into the [0, limit) range).
inline int WrapInt(const int index,
                     const int limit) {
  if (index < 0) return index+limit;
  if (index >= limit) return index-limit;
  return index;
}

}  // namespace

double SampleGrid2D::GetMinDistSqFast(const Point2D sample,
                                      const double max_min_dist_sq,
                                      Point2D* nearest) const {
  double min_dist_sq = GetToroidalDistanceSq(sample, points_[0]);
  const Point2D* nearest_ptr = &points_[0];

  const double cell_width = 1.0 / width_;

  // Each iteration of this loop checks the cells on the boundary of the
  // outer_width * outer_width grid.
  int outer_width = 1;
  int num_cells_to_check = 1;
  int starting_x_pos = sample.x * width_;
  int starting_y_pos = sample.y * width_;
  while (outer_width <= (width_ + 1)) {
    // This next loop starts at the starting_x_pos and starting_y_pos, and
    // iterates over the boundary cells, checking the distances of any point
    // in the cells. It starts out going from left to right, then turns 90
    // degrees once it reaches the corner, and so on until it gets back to the
    // original cell.
    int x_pos = starting_x_pos, y_pos = starting_y_pos;
    int x_dir = 1, y_dir = 0;
    for (int i = 0; i < num_cells_to_check; i++) {
      const int lookup_x_pos = WrapInt(x_pos, width_);
      const int lookup_y_pos = WrapInt(y_pos, width_);
      const int grid_idx = lookup_y_pos * width_ + lookup_x_pos;
      Point2D* point = sample_grid_[grid_idx];
      if (point != nullptr) {
        double dist_sq = GetToroidalDistanceSq(sample, *point);
        if (dist_sq < min_dist_sq) {
          min_dist_sq = dist_sq;
          nearest_ptr = point;
          if (min_dist_sq < max_min_dist_sq) break;
        }
      }

      // Turn 90 degrees.
      if (i > 0 && (i % (outer_width-1) == 0)) {
        y_dir = -y_dir;
        std::swap(x_dir, y_dir);
      }
      x_pos += x_dir;
      y_pos += y_dir;
    }

    num_cells_to_check = 2*(2*outer_width+2);  // Outer square - inner square.
    outer_width += 2;
    starting_x_pos--;
    starting_y_pos--;

    double min_grid_dist = std::max((outer_width / 2) - 1, 0) * cell_width;
    if (min_dist_sq < (min_grid_dist*min_grid_dist) ||
        min_dist_sq < max_min_dist_sq) {
      break;
    }
  }

  if (nearest != nullptr) *nearest = *nearest_ptr;
  return min_dist_sq;
}

int SampleGrid2D::GetBestCandidate(
    const std::vector<Point2D>& candidates) const {
  int best_candidate = 0;
  Point2D prev_nearest;
  double max_min_dist_sq = GetMinDistSqFast(candidates[0], 0.0, &prev_nearest);
  const int n_candidates = candidates.size();
  for (int i = 1; i < n_candidates; i++) {
    if (GetToroidalDistanceSq(candidates[i], prev_nearest) < max_min_dist_sq)
      continue;

    double min_dist_sq =
        GetMinDistSqFast(candidates[i], max_min_dist_sq, &prev_nearest);
    if (min_dist_sq > max_min_dist_sq) {
      max_min_dist_sq = min_dist_sq;
      best_candidate = i;
    }
  }
  return best_candidate;
}

}  // namespace sampling
