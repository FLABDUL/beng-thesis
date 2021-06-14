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

#include <iostream>
#include <fstream>
#include <string>

#include <cnpy.h>

#include "madata.h"
#include "types.h"

inline cnpy::NpyArray read_npyarray(std::string input_file_path) {
   // windows fix
   std::replace(input_file_path.begin(), input_file_path.end(), '\\', '/');
   // check if file exists
   {
      std::ifstream infile(input_file_path.c_str());
      if (!infile) {
         std::cerr << "Invalid file path " << input_file_path << std::endl;
         exit(1);
      }
   }
   //std::cout << "Reading array from " << input_file_path <<std::endl;
   cnpy::NpyArray npy_array = cnpy::npy_load(input_file_path.c_str());
   return npy_array;
}

void npy2madata(std::string input_dir_path, ma_data &madata, io_parameters &params) {
   if (params.coords) {
      //std::cout << "Reading coords array..." << std::endl;

      cnpy::NpyArray npy_array = read_npyarray(input_dir_path + "/coords.npy");//read data loads byte-by-byte
      //float* coords_carray = npy_array.data<float>();
      double* coords_carray = npy_array.data<double>();//
	  
      madata.coords.reset(new PointCloud);
      madata.coords->reserve(npy_array.shape[0]);

      //std::cout << "shape 0 = " << npy_array.shape[0] << std::endl;//8 vertices
      //std::cout << "shape 1 = " << npy_array.shape[1] << std::endl;//3 per
	  
      for (size_t i = 0; i < npy_array.shape[0]; i++){
          //std::cout << "DEBUG: x =  " << std::to_string(coords_carray[i * 3 + 0]) << std::endl << std::flush;
         //std::cout << "DEBUG: y =  " << std::to_string(coords_carray[i * 3 + 1]) << std::endl << std::flush;
         //std::cout << "DEBUG: z =  " << std::to_string(coords_carray[i * 3 + 2]) << std::endl << std::flush;
         madata.coords->push_back(Point(
            (float)coords_carray[i * 3 + 0],
            (float)coords_carray[i * 3 + 1],
            (float)coords_carray[i * 3 + 2]
         ));
	  }
      //npy_array.destruct();
   }

   if (params.normals) {
      //std::cout << "Reading normals array..." << std::endl;

      cnpy::NpyArray npy_array_norm = read_npyarray(input_dir_path + "/normals.npy");
      double* normals_carray = npy_array_norm.data<double>();

      if (npy_array_norm.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and normals" << std::endl;
         exit(1);
      }

      madata.normals.reset(new NormalCloud);
      madata.normals->reserve(madata.coords->size());

      for (size_t i = 0; i < madata.coords->size(); i++){
          //std::cout << "DEBUG: nx =  " << std::to_string(normals_carray[i * 3 + 0]) << std::endl << std::flush;
         //std::cout << "DEBUG: ny =  " << std::to_string(normals_carray[i * 3 + 1]) << std::endl << std::flush;
         //std::cout << "DEBUG: nz =  " << std::to_string(normals_carray[i * 3 + 2]) << std::endl << std::flush;
         madata.normals->push_back(Normal(
            (float)normals_carray[i * 3 + 0],
            (float)normals_carray[i * 3 + 1],
            (float)normals_carray[i * 3 + 2]
         ));
	  }
      //npy_array.destruct();
   }

   if (params.ma_coords) {
      //std::cout << "Reading ma coords arrays..." << std::endl;

      cnpy::NpyArray in_npy_array = read_npyarray(input_dir_path + "/ma_coords_in.npy");
      float* in_ma_coords_carray = in_npy_array.data<float>();

      if (in_npy_array.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and inner ma coords" << std::endl;
         exit(1);
      }

      cnpy::NpyArray out_npy_array = read_npyarray(input_dir_path + "/ma_coords_out.npy");
      float* out_ma_coords_carray = out_npy_array.data<float>();

      if (out_npy_array.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and outer ma coords" << std::endl;
         exit(1);
      }

      madata.ma_coords.reset(new PointCloud);
      madata.ma_coords->reserve(2 * madata.coords->size());

      for (size_t i = 0; i < madata.coords->size(); i++){
         madata.ma_coords->push_back(Point(
            in_ma_coords_carray[i * 3 + 0],
            in_ma_coords_carray[i * 3 + 1],
            in_ma_coords_carray[i * 3 + 2]
         ));
	  }
      //in_npy_array.destruct();

      for (size_t i = 0; i < madata.coords->size(); i++)
         madata.ma_coords->push_back(Point(
            out_ma_coords_carray[i * 3 + 0],
            out_ma_coords_carray[i * 3 + 1],
            out_ma_coords_carray[i * 3 + 2]
         ));
      //out_npy_array.destruct();
   }

   if (params.ma_qidx) {
      //std::cout << "Reading q index arrays..." << std::endl;

      cnpy::NpyArray in_npy_array = read_npyarray(input_dir_path + "/ma_qidx_in.npy");
      int* in_qidx_carray = in_npy_array.data<int>();

      if (in_npy_array.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and inner q indices" << std::endl;
         exit(1);
      }

      cnpy::NpyArray out_npy_array = read_npyarray(input_dir_path + "/ma_qidx_out.npy");
      int* out_qidx_carray = out_npy_array.data<int>();

      if (out_npy_array.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and outer q indices" << std::endl;
         exit(1);
      }

      madata.ma_qidx.reserve(2 * madata.coords->size());

      for (size_t i = 0; i < madata.coords->size(); i++)
         madata.ma_qidx.push_back(in_qidx_carray[i]);
      //in_npy_array.destruct();

      for (size_t i = 0; i < madata.coords->size(); i++)
         madata.ma_qidx.push_back(out_qidx_carray[i]);
      //out_npy_array.destruct();
   }

   if (params.lfs) {
      //std::cout << "Reading lfs array..." << std::endl;

      cnpy::NpyArray npy_array = read_npyarray(input_dir_path + "/lfs.npy");
      float* lfs_carray = npy_array.data<float>();

      if (npy_array.shape[0] != madata.coords->size()) {
         std::cerr << "Mismatched number of coords and lfs" << std::endl;
         exit(1);
      }

      madata.lfs.reserve(madata.coords->size());

      for (size_t i = 0; i < madata.coords->size(); i++)
         madata.lfs.push_back(lfs_carray[i]);
      //npy_array.destruct();
   }
}

void madata2npy(std::string npy_path, ma_data &madata, io_parameters &params) {
   if (params.coords) {
      //std::cout << "Writing coords array..." << std::endl;

      //std::cout << "coords size = ." << madata.coords->size() << std::endl;
	  
      //const unsigned int shape[] = { static_cast<unsigned int>(madata.coords->size()), 3 };
	  const std::vector<size_t> shape{ static_cast<size_t>(madata.coords->size()), 3 };
      float* coords_carray = new float[madata.coords->size() * 3];
      for (size_t i = 0; i < madata.coords->size(); i++) {
         coords_carray[i * 3 + 0] = madata.coords->at(i).x;
         coords_carray[i * 3 + 1] = madata.coords->at(i).y;
         coords_carray[i * 3 + 2] = madata.coords->at(i).z;
      }
      //cnpy::npy_save(npy_path + "/coords.npy", coords_carray, shape, 2, "w");
	  cnpy::npy_save(npy_path + "/coords.npy", coords_carray, shape);
      delete[] coords_carray; coords_carray = nullptr;
   }

   if (params.normals) {
      //std::cout << "Writing normals array..." << std::endl;

      //const unsigned int shape[] = { static_cast<unsigned int>(madata.coords->size()), 3 };
      const std::vector<size_t> shape{ static_cast<size_t>(madata.coords->size()), 3 };
      float* normals_carray = new float[madata.coords->size() * 3];
      for (size_t i = 0; i < madata.coords->size(); i++) {
         normals_carray[i * 3 + 0] = madata.normals->at(i).normal_x;
         normals_carray[i * 3 + 1] = madata.normals->at(i).normal_y;
         normals_carray[i * 3 + 2] = madata.normals->at(i).normal_z;
      }
      //cnpy::npy_save(npy_path + "/normals.npy", normals_carray, shape, 2, "w");
	  cnpy::npy_save(npy_path + "/normals.npy", normals_carray, shape);
      delete[] normals_carray; normals_carray = nullptr;
   }

   if (params.ma_coords) {
      //std::cout << "Writing ma coords arrays..." << std::endl;

      //const unsigned int shape[] = { static_cast<unsigned int>(madata.coords->size()), 3 };
      const std::vector<size_t> shape{ static_cast<size_t>(madata.coords->size()), 3 };

      float* in_ma_coords_carray = new float[madata.coords->size() * 3];
      for (size_t i = 0; i < madata.coords->size(); i++) {
         in_ma_coords_carray[i * 3 + 0] = madata.ma_coords->at(i).x;
         in_ma_coords_carray[i * 3 + 1] = madata.ma_coords->at(i).y;
         in_ma_coords_carray[i * 3 + 2] = madata.ma_coords->at(i).z;
      }
      //cnpy::npy_save(npy_path + "/ma_coords_in.npy", in_ma_coords_carray, shape, 2, "w");
	  cnpy::npy_save(npy_path + "/ma_coords_in.npy", in_ma_coords_carray, shape);
      delete[] in_ma_coords_carray; in_ma_coords_carray = nullptr;

      float* out_ma_coords_carray = new float[madata.coords->size() * 3];
      for (size_t i = 0; i < madata.coords->size(); i++) {
         out_ma_coords_carray[i * 3 + 0] = madata.ma_coords->at(i + madata.coords->size()).x;
         out_ma_coords_carray[i * 3 + 1] = madata.ma_coords->at(i + madata.coords->size()).y;
         out_ma_coords_carray[i * 3 + 2] = madata.ma_coords->at(i + madata.coords->size()).z;
      }
      //cnpy::npy_save(npy_path + "/ma_coords_out.npy", out_ma_coords_carray, shape, 2, "w");
	  cnpy::npy_save(npy_path + "/ma_coords_out.npy", out_ma_coords_carray, shape);
      delete[] out_ma_coords_carray; out_ma_coords_carray = nullptr;
   }

   if (params.ma_qidx) {
      //std::cout << "Writing q index arrays..." << std::endl;

      //const unsigned int shape[] = { static_cast<unsigned int>(madata.coords->size()) };
      const std::vector<size_t> shape{ static_cast<size_t>(madata.coords->size()), 1 };

      //cnpy::npy_save(npy_path + "/ma_qidx_in.npy", &madata.ma_qidx[0], shape, 1, "w");
	  cnpy::npy_save(npy_path + "/ma_qidx_in.npy", &madata.ma_qidx[0], shape, "w");
      //cnpy::npy_save(npy_path + "/ma_qidx_out.npy", &madata.ma_qidx[madata.coords->size()], shape, 1, "w");
	  cnpy::npy_save(npy_path + "/ma_qidx_out.npy", &madata.ma_qidx[madata.coords->size()], shape);
   }
   
      if (params.ma_rs) {
      //std::cout << "Writing radius arrays..." << std::endl;

      const std::vector<size_t> shape{ static_cast<size_t>(madata.coords->size()) };

      //cnpy::npy_save(npy_path + "/ma_qidx_in.npy", &madata.ma_qidx[0], shape, 1, "w");
	  cnpy::npy_save(npy_path + "/ma_rad_in.npy", &madata.ma_rs[0], shape, "w");
      //cnpy::npy_save(npy_path + "/ma_qidx_out.npy", &madata.ma_qidx[madata.coords->size()], shape, 1, "w");
	  cnpy::npy_save(npy_path + "/ma_rad_out.npy", &madata.ma_rs[madata.coords->size()], shape);
   }

   if (params.lfs) {
      //std::cout << "Writing lfs array..." << std::endl;

      //const unsigned int shape[] = { static_cast<unsigned int>(madata.coords->size()) };
	  const std::vector<size_t> shape{ static_cast<size_t>(madata.coords->size()) };
      //cnpy::npy_save(npy_path + "/lfs.npy", &madata.lfs[0], shape, 1, "w");
	  cnpy::npy_save(npy_path + "/lfs.npy", &madata.lfs[0], shape);
   }

   if (params.mask) {
      //std::cout << "Writing mask array..." << std::endl;

	  const std::vector<size_t> shape{ static_cast<size_t>(madata.coords->size()) };
	  //const unsigned int shape[] = { static_cast<unsigned int>(madata.coords->size()) };
      bool* out_mask_carray = new bool[madata.coords->size()];
      for (size_t i = 0; i < madata.coords->size(); i++) {
         out_mask_carray[i] = madata.mask[i];
      }
      //cnpy::npy_save(npy_path + "/decimate_lfs.npy", out_mask_carray, shape, 1, "w");
	  cnpy::npy_save(npy_path + "/decimate_lfs.npy", out_mask_carray, shape);
      delete[] out_mask_carray; out_mask_carray = nullptr;
   }
}

// Just a convenience function, to call when necessary.
void convertNPYtoXYZ(std::string input_dir_path)
{
   // Read in the data:
   cnpy::NpyArray coords_npy = read_npyarray(input_dir_path + "/coords.npy");
   float* coords_carray = coords_npy.data<float>();

   unsigned int num_points = coords_npy.shape[0];
   unsigned int dim = coords_npy.shape[1];

   // Write this out to a pointcloudxyz file:
   std::string outFile(input_dir_path + "/coords.xyz");
   std::ofstream out_pointcloudxyz(outFile);
   if (!out_pointcloudxyz)
   {
      std::cerr << "Invalid file path " << outFile << std::endl;
      exit(1);
   }

   // Header
   out_pointcloudxyz << "x y z\n";

   // coords
   for (int i = 0; i < num_points; i++)
   {
      for (int j = 0; j < 3; j++)
      {
         if (j > 0) out_pointcloudxyz << " ";
         out_pointcloudxyz << coords_carray[i * 3 + j];
         //std::cout << "DEBUG: i =  " << std::to_string(i) << std::endl << std::flush;
         //std::cout << "DEBUG: j =  " << std::to_string(j) << std::endl << std::flush;
         //std::cout << "DEBUG: coords =  " << std::to_string(coords_carray[i * 3 + j]) << std::endl << std::flush;
      }
      out_pointcloudxyz << "\n";
   }
   //coords_npy.destruct();
}
