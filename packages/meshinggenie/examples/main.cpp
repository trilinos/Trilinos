
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

#include "MeshingGenie_3d.h"

double string_to_double(const std::string &s)
{			 
	std::istringstream i(s); double x;
	if (!(i >> x)) return 0; return x;			
};

void save_polygonal_mesh(std::string file_name, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
	                   std::vector<std::vector<size_t> > &faces, double dm, double tol)
{
	#pragma region Save polgonal mesh:
	std::fstream file(file_name.c_str(), std::ios::out);
	file << dm << " " << tol << " "  << x.size() << " " << faces.size() << std::endl;

	// Points
	size_t num_points(x.size());
	for (size_t i = 0; i < num_points; i++)
	{
		file << x[i] << " " << y[i] << " " << z[i] << std::endl;	
	}
	// faces
	size_t num_faces(faces.size());
	for (size_t i = 0; i < num_faces; i++)
	{
		size_t num_corners(faces[i].size());
		file << num_corners;
		for (size_t j = 0; j < num_corners; j++)
		{
			file << " " << faces[i][j];
		}						
		file << std::endl;
	}
	#pragma endregion
};

void save_face_normals(std::string file_name, std::vector<double> nx, std::vector<double> &ny, std::vector<double> &nz)
{
	#pragma region Save polgonal mesh:
	std::fstream file(file_name.c_str(), std::ios::out);
	// Points
	size_t num_faces(nx.size());
	for (size_t i = 0; i < num_faces; i++)
	{
		file << nx[i] << " " << ny[i] << " " << nz[i] << std::endl;	
	}
	#pragma endregion
};

int read_input(std::string file_name, double &dm, size_t &random_seed, double &tol, std::vector<double> &xb, std::vector<double> &yb, std::vector<double> &zb,
	                                  std::vector< std::vector<size_t> > &bfaces, 
									  std::vector< std::vector<double> > &xf, std::vector< std::vector<double> > &yf, std::vector< std::vector<double> > &zf)
{	
	#pragma region Read input from file:
	std::string line;
	std::ifstream myfile(file_name.c_str());
	size_t iline(0), num_points(0), num_faces(0), num_internal_faces(0), iface(0);
	xb.clear();
	yb.clear();
	zb.clear();
	bfaces.clear();
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
				random_seed = size_t(string_to_double(tokens[1]));
				tol = string_to_double(tokens[2]);
				num_points = size_t(string_to_double(tokens[3]));
				num_faces = size_t(string_to_double(tokens[4]));
				num_internal_faces = size_t(string_to_double(tokens[5]));
				xb.reserve(num_points);
				yb.reserve(num_points);
				zb.reserve(num_points);
				bfaces.reserve(num_faces);
				xf.resize(num_internal_faces);
				yf.resize(num_internal_faces);
				zf.resize(num_internal_faces);
				continue;
			}
			else if (iline >= 2 && iline - 2 < num_points)
			{
				xb.push_back(string_to_double(tokens[0]));
				yb.push_back(string_to_double(tokens[1]));
				zb.push_back(string_to_double(tokens[2]));				
				continue;
			}
			else
			{
				if (bfaces.size() < num_faces)
				{
					size_t num_face_nodes(size_t(string_to_double(tokens[0])));
					std::vector<size_t> face(num_face_nodes);
					for (size_t i = 1; i <= num_face_nodes; i++)
					{
						face[i - 1] = (size_t) string_to_double(tokens[i]);
					}
					bfaces.push_back(face);
				}
				else
				{
					size_t num_face_nodes(size_t(string_to_double(tokens[0])));
					std::vector<double> face(num_face_nodes * 3);
					for (size_t i = 0; i < num_face_nodes; i++)
					{
						xf[iface].push_back(string_to_double(tokens[i * 3 + 1]));
						yf[iface].push_back(string_to_double(tokens[i * 3 + 2]));
						zf[iface].push_back(string_to_double(tokens[i * 3 + 3]));
					}
					iface++;
				}
				continue;
			}
		}
	}
	return 0;
	#pragma endregion
};

int read_input_stl(std::string file_name, std::vector<double> &xb, std::vector<double> &yb, std::vector<double> &zb,
	                                  std::vector< std::vector<size_t> > &bfaces)
{	
	#pragma region Read input from stl file:
	std::string line;
	std::ifstream myfile(file_name.c_str());
	size_t iline(0), num_points(0), num_faces(0);
	xb.clear();
	yb.clear();
	zb.clear();
	bfaces.clear();
	if (myfile.is_open())
	{    				
		size_t inode(0), i1(0), i2(0), i3(0);
		while (!myfile.eof())
		{					
			getline(myfile, line); iline++;	
			if (line.size() == 0) break;
           					
			std::istringstream iss(line);			
			std::vector<std::string> tokens;					
			copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
						std::back_inserter<std::vector<std::string> >(tokens));
					
			if (tokens[0] == "vertex")
			{	
				double x(string_to_double(tokens[1]));
				double y(string_to_double(tokens[2]));
				double z(string_to_double(tokens[3]));

				size_t num_nodes(xb.size()); bool found(false);
				for (size_t i = 0; i < num_nodes; i++)
				{
					if (fabs(xb[i] - x) < 1E-10 && fabs(yb[i] - y) < 1E-10 && fabs(zb[i] - z) < 1E-10) {found = true; break;}
					
					if (inode == 0) i1++;
					if (inode == 1) i2++;
					if (inode == 2) i3++;
				}
				if (!found)
				{
					xb.push_back(x);
					yb.push_back(y);
					zb.push_back(z);
				}

				inode++;
				if (inode == 3)
				{
					std::vector<size_t> face(3);
					face[0] = i1; face[1] = i2; face[2] = i3;
					bfaces.push_back(face);
					inode = 0; i1 = 0; i2 = 0; i3 = 0;
				}
			}
		}
	}
	return 0;
	#pragma endregion
};


int main(int argc, char *argv[])
{
	
	std::cout<<"This is Release I of Meshing Genie:" << std::endl;

	std::vector<double> xb;
	std::vector<double> yb;
	std::vector<double> zb;
	std::vector< std::vector<size_t> > bfaces;

	std::vector< std::vector<double> > xf;
	std::vector< std::vector<double> > yf;
	std::vector< std::vector<double> > zf;

	double r(2.0), tol(1E-6);

	//read_input("topmod-test-ascii.stl.dat", r, tol, xb, yb, zb, bfaces);
	//read_input_stl("topmod-test-ascii.stl", xb, yb, zb, bfaces); // r = 0.5
	//save_polygonal_mesh("topmod-test-ascii.dat", xb, yb, zb, bfaces, r, tol);

	size_t rand_seed(1);
	read_input("Box.dat", r, rand_seed, tol, xb, yb, zb, bfaces, xf, yf, zf);
	//read_input("castle.dat", r, rand_seed, tol, xb, yb, zb, bfaces);
	MeshingGenie_3d genie = MeshingGenie_3d(xb, yb, zb, bfaces, r, tol, xf, yf, zf);

	genie.execute(rand_seed);

	char dummy_cc;/*
	std::cout<< "\n*** enter q to quit ...";
	std::cin >> dummy_cc;*/
    return 0;
};

