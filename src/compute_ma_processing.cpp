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

#include <limits>

#ifdef VERBOSEPRINT
#include <chrono>
#include <iostream>
#endif

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#ifdef VERBOSEPRINT
typedef std::chrono::high_resolution_clock Clock;//clock for computation time measurement
#endif

//==============================
//   COMPUTE MA
//==============================

const Scalar delta_convergance = 1E-7f;
const unsigned int iteration_limit = 200;
const Point nanPoint(std::numeric_limits<Scalar>::quiet_NaN(), std::numeric_limits<Scalar>::quiet_NaN(), std::numeric_limits<Scalar>::quiet_NaN());

inline Scalar compute_radius(const Vector3 &p, const Vector3 &n, const Vector3 &q) {
   // compute radius of the ball that touches points p and q and whose center falls on the normal n from p
   //std::cout << "DEBUG: Compute radius. " << std::endl; 
   Scalar d = (p - q).norm();
   //std::cout << "DEBUG: Compute radius.  d=" << d << std::endl;
   Scalar cos_theta = n.dot(p - q) / d;
   //std::cout << "DEBUG: Compute radius.  cos_theta=" << cos_theta << std::endl;
   return Scalar(d / (2 * cos_theta));
}

inline Scalar cos_angle(const Vector3 p, const Vector3 q) {//check for threshold for ball at point
   // calculate the cosine of angle between vector p and q, see http://en.wikipedia.org/wiki/Law_of_cosines#Vector_formulation
   //std::cout << "DEBUG: Find cos angle. " << std::endl;
   Scalar result = p.dot(q) / (p.norm() * q.norm());

   if (result > 1) 
   {
	   //std::cout << "DEBUG: Cos angle results > 1. " << std::endl;
	  return 1;//positive angle
   }
   else if (result < -1)
   {
	   //std::cout << "DEBUG: Cos angle results < -1. " << std::endl; 
	   return -1;//negative angle
   }
   return result;//result is scalar cos angle
}

ma_result sb_point(const ma_parameters &input_parameters, const Vector3 &p, const Vector3 &n, pcl::search::KdTree<Point>::Ptr kd_tree) {
   // calculate a medial ball for a given oriented point using the shrinking ball algorithm,
   // see https://3d.bk.tudelft.nl/rypeters/pdfs/16candg.pdf section 3.2 for details
   unsigned int j = 0;//j is the iteration counter equivalent to i in the paper
   Scalar r = input_parameters.initial_radius, d;
   Vector3 q, last_q, c_next;//point index and next centre
   int qidx = -1, qidx_next;
   Point c; c.getVector3fMap() = p - n * r;// c, paper co = -r0np + p, is center of ball

   // we can't continue if we have bad input, we won't be able to perform nearest neighbour searches
   if (!c.getVector3fMap().allFinite())//check input is finite
   {
	   //std::cout << "DEBUG: Bad input. " << std::endl;
   
      return{ nanPoint, -1 };
   }
   // Results from our search
   std::vector<int> k_indices(1);//indexes
   std::vector<Scalar> k_distances(1);//distances between points

   //std::cout << "============================================================================ " << std::endl << std::flush;

   while (true) {
      // find closest point to c
	  //std::cout << "DEBUG: Finding closest point to centre. " << std::endl << std::flush;
      kd_tree->nearestKSearch(c, 1, k_indices, k_distances);//does a k search
      qidx_next = k_indices[0];//start the index at the first k search index
      q = kd_tree->getInputCloud()->at(qidx_next).getVector3fMap();//get point cloud points from k search
      d = k_distances[0];//get distances from k search

      // this should handle all (special) cases where we want to break the loop
      // - normal case when ball no longer shrinks
      // - the case where q==p, equivalent to pi+1 = pi or p in paper, implying that Bi is empty, cleaarly touches 2 points, and medial ball found, Ci is medial axis point
      // - any duplicate point cases
	  float rdel2 = (r-delta_convergance)*(r-delta_convergance);//convergence threshold check

	  //std::cout << "DEBUG: r = " << std::to_string(r) << std::endl << std::flush;
	  //std::cout << "DEBUG: (r-del)^2 = " << std::to_string(rdel2) << std::endl << std::flush;
	  //std::cout << "DEBUG: d =  " << std::to_string(d) << std::endl << std::flush;
	  //std::cout << "DEBUG: n[0] =  " << std::to_string(n[0]) << std::endl << std::flush;
	  //std::cout << "DEBUG: n[1] =  " << std::to_string(n[1]) << std::endl << std::flush;
	  //std::cout << "DEBUG: n[2] =  " << std::to_string(n[2]) << std::endl << std::flush;
	  //std::cout << "DEBUG: p[0] =  " << std::to_string(p[0]) << std::endl << std::flush;
	  //std::cout << "DEBUG: p[1] =  " << std::to_string(p[1]) << std::endl << std::flush;
	  //std::cout << "DEBUG: p[2] =  " << std::to_string(p[2]) << std::endl << std::flush;
	  //std::cout << "DEBUG: q[0] =  " << std::to_string(q[0]) << std::endl << std::flush;
	  //std::cout << "DEBUG: q[1] =  " << std::to_string(q[1]) << std::endl << std::flush;
	  //std::cout << "DEBUG: q[2] =  " << std::to_string(q[2]) << std::endl << std::flush;
		  
      if (p==q)//check threshold and break accordingly
	  {
		  //std::cout << "BREAK: p==q" << std::endl << std::flush;//r too small, always drop out
		  if (j == 0)
			  r = -1.0;
		  else
			  //r = compute_radius(p, n, last_q);//find next radius
			  r = 0.5*(p - q).norm();
		  //std::cout << "DEBUG: compute_rad = " << r << std::endl << std::flush;
		  break;
	  }
	  	  
	  if (d >= (r-delta_convergance)*(r-delta_convergance))//check threshold and break accordingly
	  {
		  //std::cout << "BREAK: d converged" << std::endl << std::flush;//distance is same as r put in
		  if (j == 0)//if d bigger than r, shrinking closer to start point
			//r = compute_radius(p, n, q);
			  r = 0.5 * (p - q).norm();
		  else
			  //r = compute_radius(p, n, last_q);//find next radius
			  r = 0.5 * (p - last_q).norm();
		  //std::cout << "DEBUG: compute_rad = " << r << std::endl << std::flush;
		  break;
	  }
	  // Compute next ball center
      r = compute_radius(p, n, q);//find next radius
	  //std::cout << "DEBUG: compute_rad = " << r << std::endl << std::flush;
      c_next = p - n * r;//find next centre

      if (!c_next.allFinite())//check if next centre finite
	  { 
		  //std::cout << "BREAK: Next centre infinite. " << std::endl;
		break;
	  }
      // denoising
	  //see https://3d.bk.tudelft.nl/rypeters/pdfs/16candg.pdf section 3.3 for details
      if (input_parameters.denoise_preserve || input_parameters.denoise_planar) 
	  {//create denoising thresholds
         Scalar a = cos_angle(p - c_next, q - c_next);
         Scalar separation_angle = std::acos(a);//scale invariant separation angle
		 //std::cout << " denoising. " << std::endl;
         if (j == 0 && input_parameters.denoise_planar > 0 && separation_angle < input_parameters.denoise_planar)//plane detection for slightly perturbed points
		 {
			 //std::cout << "BREAK: Denoise planar. " << std::endl;
		    break;
         }
         if (j > 0 && input_parameters.denoise_preserve > 0 && (separation_angle < input_parameters.denoise_preserve && r > (q - p).norm()))//stable ball preservation to ignore noisy points
		 {
			 //std::cout << "BREAK: Denoise preserve check. " << std::endl;
		    break;
         }//break according to thresholds dpn and dps
      }
      // Stop iteration if this looks like an infinite loop
      if (j > iteration_limit)//iteration limit set above
	  {
		  //std::cout << "BREAK: Passed iteration limit. " << std::endl;
		 break;
	  }
      c.getVector3fMap() = c_next;//go to next centre
      qidx = qidx_next;//go to next index
      j++;//start next counter and point, i+1
	  last_q = q;//take the last index instead of the current one
   }

   if (j == 0 && input_parameters.nan_for_initr)//if initial point is nan
   {	 
	   //std::cout << "DEBUG: Initial point nan. " << std::endl;
     return{ nanPoint, -1, -1.0 };//return nan and -1
   }
   else//if initial point is not nan
   {
      //std::cout << "DEBUG: Not initial point nan. r= " << std::to_string(r) << std::endl;
      return{ c, qidx, r };//return the centre point and index
   }
}

void sb_points(ma_parameters &input_parameters, ma_data &madata, bool inner, progress_callback callback) {
   // outer mat should be written to second half of ma_coords/ma_qidx
   size_t offset = 0;
   if (inner == false)
      offset = madata.coords->size();

   size_t progress = offset;
   size_t accum = 0;
//#pragma omp parallel for firstprivate(accum)
   for (int i = 0; i < madata.coords->size(); i++)
   {
      Vector3 p = (*madata.coords)[i].getVector3fMap();
      Vector3 n;
      if (inner)
	  {
		  //std::cout << "DEBUG: inner " << std::endl << std::flush;
         n = (*madata.normals)[i].getNormalVector3fMap();
	  }
      else
	  {
		  //std::cout << "DEBUG: outer " << std::endl << std::flush; 
         n = -(*madata.normals)[i].getNormalVector3fMap();
	  }

	  //std::cout << "*****************************" << std::endl << std::flush;
	  //std::cout << "DEBUG: n[0] =  " << std::to_string(n[0]) << std::endl << std::flush;
	  //std::cout << "DEBUG: n[1] =  " << std::to_string(n[1]) << std::endl << std::flush;
	  //std::cout << "DEBUG: n[2] =  " << std::to_string(n[2]) << std::endl << std::flush;
	  //std::cout << "DEBUG: p[0] =  " << std::to_string(p[0]) << std::endl << std::flush;
	  //std::cout << "DEBUG: p[1] =  " << std::to_string(p[1]) << std::endl << std::flush;
	  //std::cout << "DEBUG: p[2] =  " << std::to_string(p[2]) << std::endl << std::flush;
	  
      ma_result res = sb_point(input_parameters, p, n, madata.kd_tree);

      (*madata.ma_coords)[i + offset] = res.c;
      madata.ma_qidx[i + offset] = res.qidx;
	  madata.ma_rs[i + offset] = res.r;	  

      accum++;
      if (accum == 5000)
      {
//#pragma omp critical
         {
            progress += accum;
            if (callback)
               callback(progress);
         }
         accum = 0;
      }
   }
}

void compute_masb_points(ma_parameters &input_parameters, ma_data &madata, progress_callback callback) {
#ifdef VERBOSEPRINT
   auto start_time = Clock::now();
#endif

   if (!madata.kd_tree) {
      madata.kd_tree.reset(new pcl::search::KdTree<Point>());
      madata.kd_tree->setInputCloud(madata.coords);
#ifdef VERBOSEPRINT
      auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);//measure duration
      std::cout << "Constructed kd-tree in " << elapsed_time.count() << " ms" << std::endl;//output time taken for kd-tree construction
      start_time = Clock::now();
#endif
   }

   // Inside processing
   sb_points(input_parameters, madata, 1, callback);
#ifdef VERBOSEPRINT
   auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);//measure duration
   std::cout << "Done shrinking interior balls, took " << elapsed_time.count() << " ms" << std::endl;//print time for interior ball
   start_time = Clock::now();
#endif

   // Outside processing
   sb_points(input_parameters, madata, 0, callback);
#ifdef VERBOSEPRINT
   elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start_time);//measure duration
   std::cout << "Done shrinking exterior balls, took " << elapsed_time.count() << " ms" << std::endl;//time for exterior ball
#endif
}

