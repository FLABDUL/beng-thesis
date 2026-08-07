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

#include "io.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <cnpy.h>

namespace {

cnpy::NpyArray read_npy_array(const std::filesystem::path &path) {
   if (!std::filesystem::is_regular_file(path)) {
      throw std::runtime_error("Missing NumPy array: " + path.string());
   }
   cnpy::NpyArray array = cnpy::npy_load(path.string());
   if (array.fortran_order) {
      throw std::runtime_error("Fortran-order arrays are not supported: " + path.string());
   }
   return array;
}

void require_shape(const cnpy::NpyArray &array, size_t rows, size_t columns, const std::string &name) {
   if (array.shape.size() != 2 || array.shape[0] != rows || array.shape[1] != columns) {
      throw std::runtime_error(name + " must have shape " + std::to_string(rows) + "x" + std::to_string(columns));
   }
}

void require_vector_length(const cnpy::NpyArray &array, size_t length, const std::string &name) {
   const bool one_dimensional = array.shape.size() == 1 && array.shape[0] == length;
   const bool column_vector = array.shape.size() == 2 && array.shape[0] == length && array.shape[1] == 1;
   if (!one_dimensional && !column_vector) {
      throw std::runtime_error(name + " must contain " + std::to_string(length) + " values");
   }
}

Scalar floating_value(const cnpy::NpyArray &array, size_t index, const std::string &name) {
   if (array.word_size == sizeof(float)) {
      return static_cast<Scalar>(array.data<float>()[index]);
   }
   if (array.word_size == sizeof(double)) {
      return static_cast<Scalar>(array.data<double>()[index]);
   }
   throw std::runtime_error(name + " must use float32 or float64 values");
}

std::vector<float> flatten_points(const PointCloud &points, size_t offset, size_t count) {
   std::vector<float> values(count * 3);
   for (size_t i = 0; i < count; ++i) {
      const Point &point = points.at(i + offset);
      values[i * 3] = point.x;
      values[i * 3 + 1] = point.y;
      values[i * 3 + 2] = point.z;
   }
   return values;
}

}  // namespace

void npy2madata(const std::string &input_dir_path, ma_data &madata, const io_parameters &params) {
   const std::filesystem::path input_dir(input_dir_path);

   if (params.coords) {
      cnpy::NpyArray coords = read_npy_array(input_dir / "coords.npy");
      if (coords.shape.size() != 2 || coords.shape[1] != 3) {
         throw std::runtime_error("coords.npy must have shape Nx3");
      }
      madata.coords.reset(new PointCloud);
      madata.coords->reserve(coords.shape[0]);
      for (size_t i = 0; i < coords.shape[0]; ++i) {
         madata.coords->push_back(Point(
            floating_value(coords, i * 3, "coords.npy"),
            floating_value(coords, i * 3 + 1, "coords.npy"),
            floating_value(coords, i * 3 + 2, "coords.npy")));
      }
   }

   if ((params.normals || params.ma_coords || params.ma_qidx || params.lfs) && !madata.coords) {
      throw std::runtime_error("Coordinate data must be loaded before dependent arrays");
   }

   const size_t point_count = madata.coords ? madata.coords->size() : 0;

   if (params.normals) {
      cnpy::NpyArray normals = read_npy_array(input_dir / "normals.npy");
      require_shape(normals, point_count, 3, "normals.npy");
      madata.normals.reset(new NormalCloud);
      madata.normals->reserve(point_count);
      for (size_t i = 0; i < point_count; ++i) {
         madata.normals->push_back(Normal(
            floating_value(normals, i * 3, "normals.npy"),
            floating_value(normals, i * 3 + 1, "normals.npy"),
            floating_value(normals, i * 3 + 2, "normals.npy")));
      }
   }

   if (params.ma_coords) {
      cnpy::NpyArray inner = read_npy_array(input_dir / "ma_coords_in.npy");
      cnpy::NpyArray outer = read_npy_array(input_dir / "ma_coords_out.npy");
      require_shape(inner, point_count, 3, "ma_coords_in.npy");
      require_shape(outer, point_count, 3, "ma_coords_out.npy");

      madata.ma_coords.reset(new PointCloud);
      madata.ma_coords->reserve(2 * point_count);
      for (const auto *array : {&inner, &outer}) {
         for (size_t i = 0; i < point_count; ++i) {
            madata.ma_coords->push_back(Point(
               floating_value(*array, i * 3, "medial-axis coordinates"),
               floating_value(*array, i * 3 + 1, "medial-axis coordinates"),
               floating_value(*array, i * 3 + 2, "medial-axis coordinates")));
         }
      }
   }

   if (params.ma_qidx) {
      cnpy::NpyArray inner = read_npy_array(input_dir / "ma_qidx_in.npy");
      cnpy::NpyArray outer = read_npy_array(input_dir / "ma_qidx_out.npy");
      require_vector_length(inner, point_count, "ma_qidx_in.npy");
      require_vector_length(outer, point_count, "ma_qidx_out.npy");
      if (inner.word_size != sizeof(int) || outer.word_size != sizeof(int)) {
         throw std::runtime_error("Medial-axis index arrays must use 32-bit integers");
      }
      madata.ma_qidx.reserve(2 * point_count);
      madata.ma_qidx.insert(madata.ma_qidx.end(), inner.data<int>(), inner.data<int>() + point_count);
      madata.ma_qidx.insert(madata.ma_qidx.end(), outer.data<int>(), outer.data<int>() + point_count);
   }

   if (params.lfs) {
      cnpy::NpyArray lfs = read_npy_array(input_dir / "lfs.npy");
      require_vector_length(lfs, point_count, "lfs.npy");
      madata.lfs.reserve(point_count);
      for (size_t i = 0; i < point_count; ++i) {
         madata.lfs.push_back(floating_value(lfs, i, "lfs.npy"));
      }
   }
}

void madata2npy(const std::string &npy_path, const ma_data &madata, const io_parameters &params) {
   const std::filesystem::path output_dir(npy_path);
   std::filesystem::create_directories(output_dir);
   const size_t point_count = madata.coords ? madata.coords->size() : 0;
   if (point_count == 0) {
      throw std::runtime_error("Cannot write output for an empty point cloud");
   }

   const std::vector<size_t> point_shape{point_count, 3};
   const std::vector<size_t> vector_shape{point_count};

   if (params.coords) {
      std::vector<float> values = flatten_points(*madata.coords, 0, point_count);
      cnpy::npy_save((output_dir / "coords.npy").string(), values.data(), point_shape, "w");
   }

   if (params.normals) {
      if (!madata.normals || madata.normals->size() != point_count) {
         throw std::runtime_error("Normals are missing or do not match the point count");
      }
      std::vector<float> values(point_count * 3);
      for (size_t i = 0; i < point_count; ++i) {
         const Normal &normal = madata.normals->at(i);
         values[i * 3] = normal.normal_x;
         values[i * 3 + 1] = normal.normal_y;
         values[i * 3 + 2] = normal.normal_z;
      }
      cnpy::npy_save((output_dir / "normals.npy").string(), values.data(), point_shape, "w");
   }

   if (params.ma_coords) {
      if (!madata.ma_coords || madata.ma_coords->size() != 2 * point_count) {
         throw std::runtime_error("Medial-axis coordinates are missing or incomplete");
      }
      std::vector<float> inner = flatten_points(*madata.ma_coords, 0, point_count);
      std::vector<float> outer = flatten_points(*madata.ma_coords, point_count, point_count);
      cnpy::npy_save((output_dir / "ma_coords_in.npy").string(), inner.data(), point_shape, "w");
      cnpy::npy_save((output_dir / "ma_coords_out.npy").string(), outer.data(), point_shape, "w");
   }

   if (params.ma_qidx) {
      if (madata.ma_qidx.size() != 2 * point_count) {
         throw std::runtime_error("Medial-axis indices are missing or incomplete");
      }
      cnpy::npy_save((output_dir / "ma_qidx_in.npy").string(), madata.ma_qidx.data(), vector_shape, "w");
      cnpy::npy_save((output_dir / "ma_qidx_out.npy").string(), madata.ma_qidx.data() + point_count, vector_shape, "w");
   }

   if (params.ma_rs) {
      if (madata.ma_rs.size() != 2 * point_count) {
         throw std::runtime_error("Medial-ball radii are missing or incomplete");
      }
      cnpy::npy_save((output_dir / "ma_rad_in.npy").string(), madata.ma_rs.data(), vector_shape, "w");
      cnpy::npy_save((output_dir / "ma_rad_out.npy").string(), madata.ma_rs.data() + point_count, vector_shape, "w");
   }

   if (params.lfs) {
      if (madata.lfs.size() != point_count) {
         throw std::runtime_error("Local feature-size values are missing or incomplete");
      }
      cnpy::npy_save((output_dir / "lfs.npy").string(), madata.lfs.data(), vector_shape, "w");
   }

   if (params.mask) {
      if (madata.mask.size() != point_count) {
         throw std::runtime_error("Simplification mask is missing or incomplete");
      }
      std::unique_ptr<bool[]> mask(new bool[point_count]);
      for (size_t i = 0; i < point_count; ++i) {
         mask[i] = madata.mask[i];
      }
      cnpy::npy_save((output_dir / "decimate_lfs.npy").string(), mask.get(), vector_shape, "w");
   }
}

void convertNPYtoXYZ(const std::string &input_dir_path) {
   const std::filesystem::path input_dir(input_dir_path);
   cnpy::NpyArray coords = read_npy_array(input_dir / "coords.npy");
   if (coords.shape.size() != 2 || coords.shape[1] != 3) {
      throw std::runtime_error("coords.npy must have shape Nx3");
   }

   const std::filesystem::path output_path = input_dir / "coords.xyz";
   std::ofstream output(output_path);
   if (!output) {
      throw std::runtime_error("Could not write " + output_path.string());
   }

   output << "x y z\n";
   for (size_t i = 0; i < coords.shape[0]; ++i) {
      output
         << floating_value(coords, i * 3, "coords.npy") << ' '
         << floating_value(coords, i * 3 + 1, "coords.npy") << ' '
         << floating_value(coords, i * 3 + 2, "coords.npy") << '\n';
   }
}
