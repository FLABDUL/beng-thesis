/*
Copyright (c) 2016 Ravi Peters

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "compute_ma_processing.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

#ifdef VERBOSEPRINT
#include <chrono>
#include <iostream>
using Clock = std::chrono::high_resolution_clock;
#endif

namespace {

const Point nan_point(
   std::numeric_limits<Scalar>::quiet_NaN(),
   std::numeric_limits<Scalar>::quiet_NaN(),
   std::numeric_limits<Scalar>::quiet_NaN());

Scalar compute_radius(const Vector3 &p, const Vector3 &n, const Vector3 &q) {
   const Scalar distance = (p - q).norm();
   if (distance <= std::numeric_limits<Scalar>::epsilon()) {
      return std::numeric_limits<Scalar>::quiet_NaN();
   }

   const Scalar cos_theta = n.dot(p - q) / distance;
   if (std::abs(cos_theta) <= std::numeric_limits<Scalar>::epsilon()) {
      return std::numeric_limits<Scalar>::quiet_NaN();
   }

   return distance / (Scalar(2) * cos_theta);
}

Scalar cos_angle(const Vector3 &lhs, const Vector3 &rhs) {
   const Scalar denominator = lhs.norm() * rhs.norm();
   if (denominator <= std::numeric_limits<Scalar>::epsilon()) {
      return Scalar(1);
   }
   return std::clamp(lhs.dot(rhs) / denominator, Scalar(-1), Scalar(1));
}

ma_result shrinking_ball_point(
   const ma_parameters &parameters,
   const Vector3 &p,
   const Vector3 &n,
   const pcl::search::KdTree<Point>::Ptr &kd_tree) {
   Scalar radius = parameters.initial_radius;
   Point center;
   center.getVector3fMap() = p - n * radius;

   if (!center.getVector3fMap().allFinite()) {
      return {nan_point, -1, std::numeric_limits<Scalar>::quiet_NaN()};
   }

   int q_index = -1;
   unsigned int iteration = 0;
   std::vector<int> nearest_indices(1);
   std::vector<Scalar> nearest_distances(1);

   while (iteration < parameters.iteration_limit) {
      if (kd_tree->nearestKSearch(center, 1, nearest_indices, nearest_distances) <= 0) {
         break;
      }

      const int next_q_index = nearest_indices.front();
      const Vector3 q = kd_tree->getInputCloud()->at(next_q_index).getVector3fMap();
      const Scalar convergence_radius = std::max(Scalar(0), radius - parameters.convergence_delta);

      // PCL returns squared distances. A ball has converged when its nearest
      // surface point lies on (or outside) the current radius.
      if (p == q || nearest_distances.front() >= convergence_radius * convergence_radius) {
         break;
      }

      const Scalar next_radius = compute_radius(p, n, q);
      if (!std::isfinite(next_radius) || next_radius <= Scalar(0) || next_radius >= radius) {
         break;
      }

      const Vector3 next_center = p - n * next_radius;
      if (!next_center.allFinite()) {
         break;
      }

      if (parameters.denoise_preserve > 0 || parameters.denoise_planar > 0) {
         const Scalar separation_angle = std::acos(cos_angle(p - next_center, q - next_center));
         if (iteration == 0 && parameters.denoise_planar > 0 &&
             separation_angle < parameters.denoise_planar) {
            break;
         }
         if (iteration > 0 && parameters.denoise_preserve > 0 &&
             separation_angle < parameters.denoise_preserve && next_radius > (q - p).norm()) {
            break;
         }
      }

      center.getVector3fMap() = next_center;
      radius = next_radius;
      q_index = next_q_index;
      ++iteration;
   }

   if (iteration == 0 && parameters.nan_for_initr) {
      return {nan_point, -1, std::numeric_limits<Scalar>::quiet_NaN()};
   }

   // Derive the reported radius from the final centre so both outputs always
   // describe the same medial ball.
   const Scalar final_radius = (p - center.getVector3fMap()).norm();
   return {center, q_index, final_radius};
}

void shrinking_ball_points(
   ma_parameters &parameters,
   ma_data &madata,
   bool inner,
   const progress_callback &callback) {
   const size_t count = madata.coords->size();
   const size_t offset = inner ? 0 : count;

#ifdef WITH_OPENMP
#pragma omp parallel for
#endif
   for (std::int64_t i = 0; i < static_cast<std::int64_t>(count); ++i) {
      const Vector3 p = (*madata.coords)[static_cast<size_t>(i)].getVector3fMap();
      const Vector3 n = inner
         ? (*madata.normals)[static_cast<size_t>(i)].getNormalVector3fMap()
         : -(*madata.normals)[static_cast<size_t>(i)].getNormalVector3fMap();

      const ma_result result = shrinking_ball_point(parameters, p, n, madata.kd_tree);
      const size_t output_index = static_cast<size_t>(i) + offset;
      (*madata.ma_coords)[output_index] = result.c;
      madata.ma_qidx[output_index] = result.qidx;
      madata.ma_rs[output_index] = result.r;
   }

   if (callback) {
      callback(offset + count);
   }
}

}  // namespace

void compute_masb_points(ma_parameters &parameters, ma_data &madata, progress_callback callback) {
#ifdef VERBOSEPRINT
   auto start_time = Clock::now();
#endif

   if (!madata.kd_tree) {
      madata.kd_tree.reset(new pcl::search::KdTree<Point>());
      madata.kd_tree->setInputCloud(madata.coords);
#ifdef VERBOSEPRINT
      auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
      std::cout << "Constructed kd-tree in " << elapsed.count() << " ms\n";
      start_time = Clock::now();
#endif
   }

   shrinking_ball_points(parameters, madata, true, callback);
#ifdef VERBOSEPRINT
   auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
   std::cout << "Shrank interior balls in " << elapsed.count() << " ms\n";
   start_time = Clock::now();
#endif

   shrinking_ball_points(parameters, madata, false, callback);
#ifdef VERBOSEPRINT
   elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);
   std::cout << "Shrank exterior balls in " << elapsed.count() << " ms\n";
#endif
}
