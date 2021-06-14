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

#include <iostream>
#include <fstream>
#include <string>

#include <tclap/CmdLine.h>

#include "compute_ma_processing.h"
#include "io.h"
#include "madata.h"
#include "types.h"

int main(int argc, char **argv) {//?
   // parse command line arguments
   
   try {//all the outputs here can be seen in the Jupyter Notebook
	   
      //std::cout << "DEBUG: Starting compute_ma.cpp. " << std::endl;//indicate executable starting
	   
      TCLAP::CmdLine cmd("Computes a MAT point approximation, see also https://github.com/tudelft3d/masbcpp", ' ', "0.1");//command line info

      TCLAP::UnlabeledValueArg<std::string> inputArg("input", "path to directory with inside it a 'coords.npy' and a 'normals.npy' file. Both should be Nx3 float arrays where N is the number of input points.", true, "", "input dir", cmd);//input help/error info
      TCLAP::UnlabeledValueArg<std::string> outputArg("output", "path to output directory", false, "", "output dir", cmd);//output help/error info

      TCLAP::ValueArg<double> denoise_preserveArg("d", "preserve", "denoise preserve threshold", false, 20, "double", cmd);//dp var
      TCLAP::ValueArg<double> denoise_planarArg("p", "planar", "denoise planar threshold", false, 32, "double", cmd);//dpt var
      TCLAP::ValueArg<double> initial_radiusArg("r", "radius", "initial ball radius", false, 200, "double", cmd);//ibr var

      TCLAP::SwitchArg nan_for_initrSwitch("a", "nan", "write nan for points with radius equal to initial radius", cmd, false);//a help/error info

      cmd.parse(argc, argv);//?

      ma_parameters input_parameters;//state vars

      input_parameters.initial_radius = float(initial_radiusArg.getValue());//get ir var
      input_parameters.denoise_preserve = (M_PI / 180.0) * denoise_preserveArg.getValue();//get dps var
      input_parameters.denoise_planar = (M_PI / 180.0) * denoise_planarArg.getValue();//get dpn var
      input_parameters.nan_for_initr = nan_for_initrSwitch.getValue();//get nan var

      std::string output_path = outputArg.isSet() ? outputArg.getValue() : inputArg.getValue();//check path i/o

      //std::cout << "Parameters: denoise_preserve=" << denoise_preserveArg.getValue() << ", denoise_planar=" << denoise_planarArg.getValue() << ", initial_radius=" << input_parameters.initial_radius << "\n";//output set vars by cmd line

      //std::cout << "DEBUG: Initial parameters created/set " << std::endl;//indicate pars given/set

      io_parameters io_params = {};//empty i/p par var?
      io_params.coords = true;//take coords
      io_params.normals = true;//take normals

      ma_data madata = {};//empty data output?
      npy2madata(inputArg.getValue(), madata, io_params);//?

      // Perform the actual processing
      madata.ma_coords.reset(new PointCloud);//create point cloud
      madata.ma_coords->resize(2 * madata.coords->size());//resize coords
      madata.ma_qidx.resize(2 * madata.coords->size());//resize index
	  madata.ma_rs.resize(2 * madata.coords->size());//resize index
      compute_masb_points(input_parameters, madata);//compute ma using ip and resized data

      io_params.coords = false;//check coords taken?
      io_params.normals = false;//check normals taken?
      io_params.ma_coords = true;//check output coords created?
      io_params.ma_qidx = true;//check output array created?
      io_params.ma_rs = true;//check output array created?

      io_params.c_p = true;//HAKIM FINDING C-P

      madata2npy(output_path, madata, io_params);//?

      {
          //std::cout << "DEBUG: Errors are printed after this line. " << std::endl;//potential error warnings
		 
         std::string output_path_metadata = output_path + "/compute_ma";//?
         std::replace(output_path_metadata.begin(), output_path_metadata.end(), '\\', '/');//?

         std::ofstream metadata(output_path_metadata.c_str());//?
         if (!metadata) {//check for invalid data file path
            throw TCLAP::ArgParseException("invalid filepath", output_path);//invalid file path warning
         }

         metadata//data from cmd line input
            << "initial_radius " << input_parameters.initial_radius << std::endl//ir ip data
            << "nan_for_initr " << input_parameters.nan_for_initr << std::endl//ir nan ip data
            << "denoise_preserve " << denoise_preserveArg.getValue() << std::endl//dps input data
            << "denoise_planar " << denoise_planarArg.getValue() << std::endl;//dpn input data
         metadata.close();
      }
   }
   catch (TCLAP::ArgException &e) { std::cerr << "Error: " << e.error() << " for " << e.argId() << std::endl; }//check for error and state what it is

   return 0;//compulsory return value
}
