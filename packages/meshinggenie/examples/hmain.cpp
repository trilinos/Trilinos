
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

// R5.0

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>


// Only use headers that are required in the testing functions
#include "HybridMesher_2d.h"

double string_to_double(const std::string &s)
{			 
	std::istringstream i(s); double x;
	if (!(i >> x)) return 0; return x;			
}

int plot_tessellation(std::string file_name, std::vector<double> &x, std::vector<double> &y,
		      std::vector< std::vector<size_t> > &elements, double dm,
		      std::vector<double> &VorBound, std::vector< std::vector<double> > &Holes, std::vector< std::vector<double> > &Cracks)
{
#pragma region Plot Voronoi mesh:
  std::fstream file(file_name.c_str(), std::ios::out);
  file << "%!PS-Adobe-3.0" << std::endl;
  file << "72 72 scale     % one unit = one inch" << std::endl;

#pragma region Retrieve bounding box, scale, and translate:

  double xmin(x[0]), ymin(y[0]), xmax(x[0]), ymax(y[0]);
  size_t num_points(x.size());
  for (size_t ii = 1; ii < num_points; ii++)
    {
      if (x[ii] < xmin) xmin = x[ii];
      if (x[ii] > xmax) xmax = x[ii];
      if (y[ii] < ymin) ymin = y[ii];
      if (y[ii] > ymax) ymax = y[ii];
    }

  double Lx(xmax - xmin), Ly(ymax - ymin);

  double scale_x, scale_y, scale;
  double shift_x, shift_y;

  scale_x = 6.5 / Lx;
  scale_y = 9.0 / Ly;

  if (scale_x < scale_y)
    {
      scale = scale_x;
      shift_x = 1.0 - xmin * scale;
      shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
    }
  else
    {
      scale = scale_y;
      shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
      shift_y = 1.0 - ymin * scale;
    }
  file << shift_x << " " << shift_y << " translate" << std::endl;
    #pragma endregion

#pragma region Definitions of Shapes:

  file << "/circ    % stack: x y r" << std::endl;
  file << "{0 360 arc" << std::endl;
  file << " closepath" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/fcirc    % stack: x y r" << std::endl;
  file << "{0 360 arc" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.6 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/seg2      % stack: x1 y1 x2 y2" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/v2_bold      % stack: x1 y1 x2 y2" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " 0.01 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/v3      % stack: x1 y1 x2 y2 x3 y3" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  for (int i = 0; i < 2; i++) file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.9 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/v4      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  for (int i = 0; i < 3; i++) file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.9 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/v5      % stack: x1 y1 x2 y2 x3 y3 x4 y4 x5 y5" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  for (int i = 0; i < 4; i++) file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.9 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/v6      % stack: x1 y1 x2 y2 x3 y3 x4 y4 x5 y5 x6 y6" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  for (int i = 0; i < 5; i++) file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.9 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/v7      % stack: x1 y1 x2 y2 x3 y3 x4 y4 x5 y5 x6 y6 x7 y7" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  for (int i = 0; i < 6; i++) file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.9 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/v8      % stack: x1 y1 x2 y2 x3 y3 x4 y4 x5 y5 x6 y6 x7 y7 x8 y8" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  for (int i = 0; i < 7; i++) file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.9 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/v9      % stack: x1 y1 x2 y2 x3 y3 x4 y4 x5 y5 x6 y6 x7 y7 x8 y8 x9 y9" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  for (int i = 0; i < 8; i++) file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.9 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/v10      % stack: x1 y1 x2 y2 x3 y3 x4 y4 x5 y5 x6 y6 x7 y7 x8 y8 x9 y9 x10 y10" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  for (int i = 0; i < 9; i++) file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.9 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/v11      % stack: x1 y1 x2 y2 x3 y3 x4 y4 x5 y5 x6 y6 x7 y7 x8 y8 x9 y9 x10 y10 x11 y11" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  for (int i = 0; i < 10; i++) file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.9 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;

  file << "/v12      % stack: x1 y1 x2 y2 x3 y3 x4 y4 x5 y5 x6 y6 x7 y7 x8 y8 x9 y9 x10 y10 x11 y11 x12 y12" << std::endl;
  file << "{newpath" << std::endl;
  file << " moveto" << std::endl;
  for (int i = 0; i < 11; i++) file << " lineto" << std::endl;
  file << " closepath" << std::endl;
  file << " gsave" << std::endl;
  file << " 0.9 setgray fill" << std::endl;
  file << " grestore" << std::endl;
  file << " 0.002 setlinewidth" << std::endl;
  file << " stroke" << std::endl;
  file << "} def" << std::endl;
   #pragma endregion

  size_t num_elements(elements.size());
  for (size_t ie = 0; ie < num_elements; ie++)
    {
#pragma region Plot Element:
      size_t num_nodes(elements[ie].size());
      if (num_nodes < 3 || num_nodes > 12)
	{
	  std::cout<< "Error!! Number of sides of Voronoi Cell = " << num_nodes << std::endl;
	  continue;
	}

      for (size_t ii = 0; ii < num_nodes; ii++)
	{
	  file << x[elements[ie][ii]] * scale << "  " << y[elements[ie][ii]] * scale << "  ";
	}
      if       (num_nodes == 3)       file << "v3"       << std::endl;
      else if  (num_nodes == 4)       file << "v4"       << std::endl;
      else if  (num_nodes == 5)       file << "v5"       << std::endl;
      else if  (num_nodes == 6)       file << "v6"       << std::endl;
      else if  (num_nodes == 7)       file << "v7"       << std::endl;
      else if  (num_nodes == 8)       file << "v8"       << std::endl;
      else if  (num_nodes == 9)       file << "v9"       << std::endl;
      else if  (num_nodes == 10)      file << "v10"      << std::endl;
      else if  (num_nodes == 11)      file << "v11"      << std::endl;
      else if  (num_nodes == 12)      file << "v12"      << std::endl;
                #pragma endregion
    }

#pragma region Boundary Edges:
  if (true)
    {
      size_t num = VorBound.size();
      for (size_t i = 0; i < num; i+=2)
	{
	  size_t ip(i + 1), ipp(ip + 1), ippp(ipp + 1);
	  if (ipp == num) {ipp = 0; ippp = 1;}
	  file << VorBound[i] * scale << "  " << VorBound[ip] * scale << "  ";
	  file << VorBound[ipp] * scale << "  " << VorBound[ippp] * scale << " v2_bold" << std::endl;
	}
    }
  size_t num_H = Holes.size();
  for (size_t ih = 0; ih < num_H; ih++)
    {
      size_t num = Holes[ih].size();
      for (size_t i = 0; i < num; i+=2)
	{
	  size_t ip(i + 1), ipp(ip + 1), ippp(ipp + 1);
	  if (ipp == num) {ipp = 0; ippp = 1;}
	  file << Holes[ih][i] * scale << "  " << Holes[ih][ip] * scale << "  ";
	  file << Holes[ih][ipp] * scale << "  " << Holes[ih][ippp] * scale << " v2_bold" << std::endl;
	}
    }
  size_t num_C = Cracks.size();
  for (size_t ic = 0; ic < num_C; ic++)
    {
      size_t num = Cracks[ic].size() - 2;
      for (size_t i = 0; i < num; i+=2)
	{
	  size_t ip(i + 1), ipp(ip + 1), ippp(ipp + 1);
	  file << Cracks[ic][i] * scale << "  " << Cracks[ic][ip] * scale << "  ";
	  file << Cracks[ic][ipp] * scale << "  " << Cracks[ic][ippp] * scale << " v2_bold" << std::endl;
	}
    }
        #pragma endregion

  file << "showpage" << std::endl;
  return 0;
    #pragma endregion
}


int save_mesh(std::string file_name, std::vector<double> &x, std::vector<double> &y, std::vector< std::vector<size_t> > &elements)	          
{
	#pragma region Save Boundary Points:
	std::fstream file(file_name.c_str(), std::ios::out);
	size_t num_nodes(x.size());
	file << num_nodes << std::endl;		
	for (size_t i = 0; i < num_nodes; i++)
	{
		file << x[i] << " " << y[i] << std::endl;
	}
	size_t num_elements(elements.size());
	file << num_elements << std::endl;
	for (size_t ie = 0; ie < num_elements; ie++)
	{
		size_t num_element_nodes(elements[ie].size());
		file << num_element_nodes;
		for (size_t ii = 0; ii <  num_element_nodes; ii++)
		{
			file << " " << elements[ie][ii];		
		}
		file << std::endl;	
	}
	file << "0" << std::endl;
	return 0;		
	#pragma endregion
}

int read_input(std::string file_name, double &dm, std::vector<double> &VorBound,
	                                  std::vector< std::vector<double> > &Holes, 
	       std::vector< std::vector<double> > &Cracks, std::vector<double> &StrBound, int &NumStr)
{	
	#pragma region Read input from file:
	std::string line; size_t num_points(0);	
	std::ifstream myfile(file_name.c_str());
	size_t iline(0); size_t num_holes(0), num_cracks(0);
	VorBound.clear(); Holes.clear(); Cracks.clear();
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
				dm = string_to_double(tokens[0]);
				num_holes = size_t(string_to_double(tokens[1]));
				Holes.resize(num_holes);
				num_cracks = size_t(string_to_double(tokens[2]));
				Cracks.resize(num_cracks);
				continue;
			}
			if (iline == 2)
			{
				int num_tokens(tokens.size());
				VorBound.resize(num_tokens);
				for (int i = 0; i < num_tokens; i++)
				{
					VorBound[i] = string_to_double(tokens[i]);
				}
				continue;
			}
			if (iline >= 3 && iline < 3 + num_holes)
			{
				int num_tokens(tokens.size());
				Holes[iline - 3].resize(num_tokens);
				for (int i = 0; i < num_tokens; i++)
				{
					Holes[iline - 3][i] = string_to_double(tokens[i]);
				}
				continue;				
			}
			if (iline >= 3 + num_holes && iline < 3 + num_holes + num_cracks)
			{
				int num_tokens(tokens.size());
				Cracks[iline - 3 - num_holes].resize(num_tokens);
				for (int i = 0; i < num_tokens; i++)
				{
					Cracks[iline - 3 - num_holes][i] = string_to_double(tokens[i]);
				}
				continue;
			}
			if (iline == 3 + num_holes + num_cracks) {
			  int num_tokens(tokens.size());
			  StrBound.resize(num_tokens);
			  for (int i = 0; i < num_tokens; i++) {
			    StrBound[i] = string_to_double(tokens[i]);
			  }
			  continue;
			}
			if (iline == 3 + num_holes + num_cracks + 1) {
			  NumStr = (int)string_to_double(tokens[0]);
			  continue;
			}
		}
	}
	return 0;
	#pragma endregion
}


int main(int argc, char *argv[])
{   	
	double dm(0.1); 
	std::vector<double> VorBound, StrBound;
	std::vector< std::vector<double> > Holes;
	std::vector< std::vector<double> > Cracks;
	int NumStr;

	read_input(argv[1], dm, VorBound, Holes, Cracks, StrBound, NumStr);

	HybridMesher_2d hybrid = HybridMesher_2d(dm, VorBound, Holes, Cracks, StrBound, NumStr, atoi(argv[2]));

	hybrid.execute();

	std::vector<double> x; std::vector<double> y;
	std::vector< std::vector<size_t> > elements;
	hybrid.get_Tessellation(x, y, elements);
	plot_tessellation("final_hybrid.ps", x, y, elements, dm, VorBound, Holes, Cracks);
	save_mesh("tessellation.dat", x, y, elements);
    std::cout<< "\n*** Mission Accomplished ***." << std::endl;
    
    return 0;
}



