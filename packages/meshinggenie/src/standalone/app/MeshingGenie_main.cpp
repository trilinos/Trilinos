//@HEADER
// ************************************************************************
// 
//               MeshingGenie: Fracture Meshing Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/*******************************************************************************
 * MeshingGenie_v1p0.cpp
 * Author: Mohamed S. Ebeida (msebeid@sandia.gov)  
 * Desctiption: Main entry point for the standlone version of MeshingGenie
   This software is based on a series of publication for solving the maximal 
   Poisson-disk sampling problem and utilizing the solution in meshing. For 
   more details, please read (cite) the following articles:
   1. M. S. Ebeida, A. Patney, S. A. Mitchell, A. A. Davidson, P. M. Knupp, 
      and J. D. Owens, "Efficient maximal Poisson-disk Sampling", ACM 
	  Transactions on Graphics (SIGGRAPH 2011), 30(4), August 2011.
   2. M. S. Ebeida, S. A. Mitchell, A. A. Davidson, A. Patney, P. M. Knupp, 
      and J. D. Owens, "Efficient and Good Delaunay meshes from random points",
	  Computer Aided Design, 43(11): 1506-1515.
   3. M. S. Ebeida, P. M. Knupp, V. J. Leung, J. E. Bishop, and M. J. Martinez,
      "Mesh generation for modeling and simulation of carbon sequestration 
	  processes", DOE Scientific Discovery through Advanced Computing (SciDAC) 
	  conference, July 10-14, Denver, 2011. 
   4. M. S. Ebeida, S. A. Mitchell, "Uniform random Voronoi meshes", 20th 
      International Meshing Roundtable, Paris October 23-26, 2011.
   5. M. S. Ebeida, S. A. Mitchell, A. Patney, A. A. Davidson, and J. D. Owens, 
      "A simple algorithm for maximal Poisson-disk sampling in high dimensions", 
	  Computer Graphics Forum (Eurographics 2012), 31(2), May 2012.
   6. S. A. Mitchell, A. Rand, M. S. Ebeida and C. Bajaj, "Variable radii 
      Poisson-disk sampling", 24th Canadian Conference on Computational Geometry
	  (CCCG'12), Charlottetown, August 2012.
 * Last modified: 11/21/2012
********************************************************************************/

#include "MeshingGenie_defs.h"
#include "MeshingGenie_mps_nd.h"

int read_input(std::string file_name);

double string_to_double(const std::string &s);

// Input variables
std::vector< std::vector<double> > _boundary_points; // coordinates of boundary points    
std::vector< std::vector<size_t> > _boundary_faces;  // planar boundary faces
size_t _ndim;                                        // Number of dimensions
double _r;                                           // distribution radius
size_t _random_seed;                                 // random seed
bool   _mps;                                         // Solve MPS
bool   _delaunay_tessellation;                       // Generate Delaunay Tessellation
bool   _voronoi_tessellation;                        // Generate Voronoi Tessellation

// Ourtput variables
std::vector< std::vector<double> > _sample_points;  // coordinates of output sample points
std::vector< std::vector<double> > _v_corners;      // coordinates of output voronoi corners
std::vector< std::vector<size_t> > _dtessellation;  // Delaunay tessellations
std::vector< std::vector< std::vector< size_t> > > _vtessellation;  // Voronoi tessellations

int main(int argc, char *argv[])
{
	_r = 0.0; _ndim = 0;
	_mps = false;
	_delaunay_tessellation = false;
	_voronoi_tessellation = false;

	if (argc == 1)
	{
		std::cout<< "Meshing Genie: Please specify input file!" << std::endl;
		char dummy_cc;
		std::cout<< "\n*** enter q to quit ...";
		std::cin >> dummy_cc;
		return 1;
	}

	if(read_input(argv[1])) 
	{
		char dummy_cc;
		std::cout<< "\n*** Error reading input file, enter q to quit ...";
		std::cin >> dummy_cc;
		return 1;
	}

	if (_mps)
	{
		MeshingGenie_mps_nd genie_mps;
		genie_mps.solve_mps(_ndim, _r, _random_seed, _boundary_points, _boundary_faces, _sample_points);
		genie_mps.save_maximal_sample("maximal_sample.dat");
	}

	char dummy_cc;
	std::cout<< "\n*** enter q to quit ...";
    std::cin >> dummy_cc;

	return 0;
}

double string_to_double(const std::string &s)
{			 
	std::istringstream i(s); double x;
	if (!(i >> x)) return 0.0; return x;			
}

int read_input(std::string file_name)
{	
	#pragma region Read input from file:
	size_t iline(0); size_t num_bp(0), num_bf(0), i_bp(0), i_bf(0);
	std::string line;
	std::ifstream myfile(file_name.c_str());
	if (myfile.is_open())
	{    				
		while (!myfile.eof())
		{					
			getline(myfile, line); iline++;	
			if (line.size() == 0) break;
           					
			std::istringstream iss(line);			
			std::vector<std::string> tokens;					
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
				 std::back_inserter<std::vector<std::string> >(tokens));
					
			if (iline == 1)
			{		
				std::string str = tokens[0];
				if (str.find("p")!=str.npos) _mps = true;
				if (str.find("d")!=str.npos) _delaunay_tessellation = true;
				if (str.find("v")!=str.npos) _voronoi_tessellation = true;
			}
			else if (iline == 2)
			{
				_ndim = (size_t) string_to_double(tokens[0]);
				_r = string_to_double(tokens[1]);
				_random_seed = (size_t) string_to_double(tokens[2]);
				num_bp = (size_t) string_to_double(tokens[3]);
				num_bf = (size_t) string_to_double(tokens[4]);
				_boundary_points.resize(num_bp);
				_boundary_faces.resize(num_bf);
			}
			else if (i_bp < num_bp)
			{
				for (size_t i = 0; i < _ndim; i++)
				{
					_boundary_points[i_bp].push_back(string_to_double(tokens[i]));
				}
				i_bp++;
			}
			else if (i_bf < num_bf)
			{
				for (size_t i = 0; i < _ndim; i++)
				{
					_boundary_faces[i_bf].push_back((size_t)string_to_double(tokens[i]));
				}
				i_bf++;
			}
		}
		return 0;
	}
	std::cout<< "Meshing Genie: Error reading input file!" << std::endl;
	return 1;
	#pragma endregion
}
