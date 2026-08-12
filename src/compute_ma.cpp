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

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <tclap/CmdLine.h>

#include "compute_ma_processing.h"
#include "io.h"
#include "madata.h"
#include "version.h"

int main(int argc, char **argv) {
   try {
      TCLAP::CmdLine cmd(
         "Approximate interior and exterior medial-axis balls using the shrinking-ball algorithm.",
         ' ', MASBCPP_VERSION);

      TCLAP::UnlabeledValueArg<std::string> input_arg(
         "input", "Directory containing Nx3 coords.npy and normals.npy arrays.",
         true, "", "input directory", cmd);
      TCLAP::UnlabeledValueArg<std::string> output_arg(
         "output", "Output directory; defaults to the input directory.",
         false, "", "output directory", cmd);
      TCLAP::ValueArg<double> preserve_arg(
         "d", "preserve", "Stable-ball denoising angle in degrees.", false, 20, "degrees", cmd);
      TCLAP::ValueArg<double> planar_arg(
         "p", "planar", "Planar denoising angle in degrees.", false, 32, "degrees", cmd);
      TCLAP::ValueArg<double> radius_arg(
         "r", "radius", "Initial ball radius in model units.", false, 200, "number", cmd);
      TCLAP::ValueArg<double> convergence_arg(
         "c", "convergence", "Minimum radius change before convergence.", false, 1E-7, "number", cmd);
      TCLAP::ValueArg<unsigned int> iterations_arg(
         "i", "iterations", "Maximum shrinking-ball iterations per sample.", false, 200, "integer", cmd);
      TCLAP::SwitchArg nan_arg(
         "a", "nan", "Write NaN when a ball remains at its initial radius.", cmd, false);

      cmd.parse(argc, argv);

      if (radius_arg.getValue() <= 0 || convergence_arg.getValue() <= 0 || iterations_arg.getValue() == 0) {
         throw TCLAP::ArgParseException(
            "radius, convergence and iterations must be greater than zero", "parameters");
      }

      constexpr double pi = 3.14159265358979323846;
      ma_parameters parameters{};
      parameters.initial_radius = static_cast<Scalar>(radius_arg.getValue());
      parameters.denoise_preserve = (pi / 180.0) * preserve_arg.getValue();
      parameters.denoise_planar = (pi / 180.0) * planar_arg.getValue();
      parameters.nan_for_initr = nan_arg.getValue();
      parameters.convergence_delta = static_cast<Scalar>(convergence_arg.getValue());
      parameters.iteration_limit = iterations_arg.getValue();

      const std::string output_path = output_arg.isSet() ? output_arg.getValue() : input_arg.getValue();
      std::filesystem::create_directories(output_path);

      io_parameters io_params{};
      io_params.coords = true;
      io_params.normals = true;

      ma_data madata{};
      npy2madata(input_arg.getValue(), madata, io_params);
      madata.ma_coords.reset(new PointCloud);
      madata.ma_coords->resize(2 * madata.coords->size());
      madata.ma_qidx.resize(2 * madata.coords->size());
      madata.ma_rs.resize(2 * madata.coords->size());
      compute_masb_points(parameters, madata);

      io_params = {};
      io_params.ma_coords = true;
      io_params.ma_qidx = true;
      io_params.ma_rs = true;
      madata2npy(output_path, madata, io_params);

      std::filesystem::path metadata_path = std::filesystem::path(output_path) / "compute_ma.txt";
      std::ofstream metadata(metadata_path);
      if (!metadata) {
         throw std::runtime_error("Could not write metadata to " + metadata_path.string());
      }
      metadata
         << "version " << MASBCPP_VERSION << '\n'
         << "initial_radius " << parameters.initial_radius << '\n'
         << "nan_for_initial_radius " << parameters.nan_for_initr << '\n'
         << "denoise_preserve_degrees " << preserve_arg.getValue() << '\n'
         << "denoise_planar_degrees " << planar_arg.getValue() << '\n'
         << "convergence_delta " << parameters.convergence_delta << '\n'
         << "iteration_limit " << parameters.iteration_limit << '\n';
   } catch (const TCLAP::ArgException &error) {
      std::cerr << "Error: " << error.error() << " for " << error.argId() << '\n';
      return 2;
   } catch (const std::exception &error) {
      std::cerr << "Error: " << error.what() << '\n';
      return 1;
   }

   return 0;
}
