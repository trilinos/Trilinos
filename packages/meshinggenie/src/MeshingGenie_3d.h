// 09/01/2011 5:41 pm

// R1.0: Non-convex domains (sampling + Voronoi meshing + output and testing planar faces + non-convex domains + internal faces)

// MeshingGenie_3d.h

// Author: Mohamed S. Ebeida 

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

#ifndef MESHING_GENIE_3D_H
#define MESHING_GENIE_3D_H

#include <cmath>
#include <vector>
#include <list>
#include <utility>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <cstddef>
#include <cstdlib>
#include <limits>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>

#include <time.h>

#include "MeshingGenie_2d.h"


class MeshingGenie_3d   
{
	protected:
		
		struct EdgePoint
		{
			EdgePoint() {};
			double x;
			double y;
			double z;
			EdgePoint* next;
			double h;
			double v;
			size_t inode;
		};

		struct Polyhedron;

		struct Point 
		{
			/// Default constructor does nothing (for performance).
			Point() {}
			double x; double y; double z;
			Point* next;
			Polyhedron* poly;
			bool _on_reflex_corner;
			bool _on_internal_boundaries;
		};

		struct VoronoiPoint 
		{
			/// Default constructor does nothing (for performance).
			VoronoiPoint() {}
			double x; double y; double z;
			std::set<Polyhedron*> polys;
			VoronoiPoint* nextVpoint;
			bool on_boundary;
		};

		struct FacePoint
		{
			FacePoint(){};
			bool invalid;
			FacePoint* nextcorner;
			FacePoint* prevcorner;
			VoronoiPoint* Vpoint;	
		};

		struct PolyFace
		{
			PolyFace(){};
			FacePoint* corner;
			PolyFace* nextface;
			size_t    num_corners;
			double nx; double ny; double nz;
			bool on_boundary;
			bool inverted;
			bool preserved;
			Point* q;
		};

		struct Polyhedron
		{
			/// Default constructor does nothing (for performance).
			Polyhedron() {};
			PolyFace* face;
			Polyhedron* nextpoly;
		};

		public:
        
			//! constructor        
			MeshingGenie_3d(std::vector<double> &xb, std::vector<double> &yb, std::vector<double> &zb, 
					std::vector< std::vector<size_t> > &bfaces, double r, double tol, 
					std::vector< std::vector<double> > &xf, std::vector< std::vector<double> > &yf, std::vector< std::vector<double> > &zf);
        
			//! Destructor                        
			~MeshingGenie_3d(){ };


		public:
			
			// old code

			inline void scale_input(std::vector<double> &xb, std::vector<double> &yb, std::vector<double> &zb, 
				                    double &r, double &tol, 
									std::vector< std::vector<double> > &xf, std::vector< std::vector<double> > &yf, std::vector< std::vector<double> > &zf)
			{
				#pragma region Scale Input:
				_SF = 1.0 / r;
				double xmin(xb[0]), xmax(xb[0]), ymin(yb[0]), ymax(yb[0]), zmin(zb[0]), zmax(zb[0]);
				size_t num_input_points(xb.size());
				for (size_t i = 0; i < num_input_points; i++)
				{
					if (xb[i] < xmin) xmin = xb[i];
					if (xb[i] > xmax) xmax = xb[i];
					if (yb[i] < ymin) ymin = yb[i];
					if (yb[i] > ymax) ymax = yb[i];
					if (zb[i] < zmin) zmin = zb[i];
					if (zb[i] > zmax) zmax = zb[i];
				}
				_DX = 0.5 * (xmin + xmax);
				_DY = 0.5 * (ymin + ymax);
				_DZ = 0.5 * (zmin + zmax);
				r *= _SF;
				tol *= _SF;
				for (size_t i = 0; i < num_input_points; i++)
				{
					xb[i] -= _DX; xb[i] *= _SF;
					yb[i] -= _DY; yb[i] *= _SF;
					zb[i] -= _DZ; zb[i] *= _SF;
				}
				size_t num_internal_faces(xf.size());
				for (size_t i = 0; i < num_internal_faces; i++)
				{
					size_t num_face_corners(xf[i].size());
					for (size_t j = 0; j < num_face_corners; j++)
					{
						xf[i][j] -= _DX; xf[i][j] *= _SF;
						yf[i][j] -= _DY; yf[i][j] *= _SF;
						zf[i][j] -= _DZ; zf[i][j] *= _SF;
					}
				}
				_SF_inv = 1.0 / _SF;
				#pragma endregion
			};

			inline void scale_output(double &x, double &y, double &z)
			{
				x*= _SF_inv; x+= _DX;
				y*= _SF_inv; y+= _DY;
				z*= _SF_inv; z+= _DZ;
			};

			inline void get_b_edge_faces(size_t i, size_t j, std::vector<size_t> &edge_faces)
			{
				#pragma region Retrieve boundaary edge face:
				edge_faces.clear();
				size_t num_node_faces(_b_node_faces[i].size());
				for (size_t iface = 0; iface < num_node_faces; iface++)
				{
					size_t face_index = _b_node_faces[i][iface];
					size_t num_face_node(_bfaces[face_index].size());
					for (size_t inode = 0; inode < num_face_node; inode++)
					{
						if (_bfaces[face_index][inode] == j)
						{
							edge_faces.push_back(face_index);
							break;
						}
					}
				}
				return;
				#pragma endregion
			};

			inline Point* closest_point(double xx, double yy, double zz, double tol)
			{
				#pragma region Get Closest Point:
				size_t i, j, k, icell, jcell;
				i = size_t((xx - _xo) / _s);
				j = size_t((yy - _yo) / _s);
				k = size_t((zz - _zo) / _s);
				get_cell_index(i, j, k, icell);
				Point* pnt(_cell_points[icell]);

				while (pnt != 0)
				{
					if (fabs(pnt->x - xx) < tol && fabs(pnt->y - yy) < tol && fabs(pnt->z - zz) < tol)	
					{
						return pnt;
					}
					pnt = pnt->next;
				}

				// check points in the neighbor if point lies close to the cell boundaries
				bool ip(false), jp(false), kp(false), im(false), jm(false), km(false);
				double xmin(_xo + i * _s), ymin(_yo + j * _s), zmin(_zo + k * _s), xmax(xmin + _s), ymax(ymin + _s), zmax(zmin + _s);
				if (xmax - xx < tol) ip = true;
				else if (xx - xmin < tol) im = true;
				if (ymax - yy < tol) jp = true;
				else if (yy - ymin < tol) jm = true;
				if (zmax - zz < tol) kp = true;
				else if (zz - zmin < tol) km = true;

				if (!ip && !jp && !kp && !im && !jm && !km) return 0;

				// search neighbors
				// check points in the neighbor 27 cells
				for (size_t index = 0; index < 27; index++)
				{
					jcell = icell + _neighbors[index];
					// check all point in that cell
					pnt = _cell_points[jcell];
					while (pnt != 0)
					{
						if (fabs(pnt->x - xx) < tol && fabs(pnt->y - yy) < tol && fabs(pnt->z - zz) < tol) 
						{
							return pnt;
						}
						pnt = pnt->next;
					}
				}
				return 0;
				#pragma endregion
			};

			inline void plot_polyhedron(Polyhedron* poly)
			{
				#pragma region plot Voronoi cell Faces:			
				PolyFace* f = poly->face;size_t iface(0);
				std::vector<double> x;
				std::vector<double> y;					
				std::vector<double> z;
				std::vector< std::vector<size_t> > faces;
				while (f != 0)
				{				
					size_t num_face_corners(f->num_corners);
					if (num_face_corners < 3)
					{
						f = f->nextface;
						continue;
					}
					std::vector<size_t> face;
					FacePoint* pnt = f->corner;
					for (size_t inode = 0; inode < num_face_corners; inode++)
					{
						double xx(pnt->Vpoint->x);
						double yy(pnt->Vpoint->y);
						double zz(pnt->Vpoint->z);
						size_t ii(0);
						for (size_t i = 0; i < x.size(); i++)
						{
							if (fabs(xx - x[i]) < 1E-10 && fabs(yy - y[i]) < 1E-10 && fabs(zz - z[i]) < 1E-10) break;							
							ii++;							
						}
						if (ii == x.size())
						{
							x.push_back(xx); y.push_back(yy); z.push_back(zz);
						}
						face.push_back(ii);
						pnt = pnt->nextcorner;
					}	
					faces.push_back(face);
					f = f->nextface;
				}
				save_vtk_polygons(x, y, z, faces);
				f = f;
				#pragma endregion
			};

			inline void plot_cell_faces(Polyhedron* poly)
			{
				#pragma region plot Voronoi cell Faces:			
				PolyFace* f = poly->face;size_t iface(0);
				while (f != 0)
				{
					std::vector<double> xx;
					std::vector<double> yy;					
					std::vector<double> zz;
					std::vector<size_t> face;
			
					FacePoint* pnt = f->corner;
					for (size_t inode = 0; inode < f->num_corners; inode++)
					{
						xx.push_back(pnt->Vpoint->x);
						yy.push_back(pnt->Vpoint->y);
						zz.push_back(pnt->Vpoint->z);
						face.push_back(inode);
						pnt = pnt->nextcorner;
					}
					if (f->num_corners > 2) save_vtk_face(xx, yy, zz, face, iface);
					f = f->nextface;				
					iface++;
				}
				f = f;
				#pragma endregion
			};

			inline void plot_cell_face(PolyFace* f, size_t iface)
			{
				#pragma region plot Voronoi cell Face:			
				if (f == 0) return;
				if (f->num_corners == 0) return;
				std::vector<double> xx;
				std::vector<double> yy;					
				std::vector<double> zz;
				std::vector<size_t> face;
			
				FacePoint* pnt = f->corner;
				for (size_t inode = 0; inode < f->num_corners; inode++)
				{
					xx.push_back(pnt->Vpoint->x);
					yy.push_back(pnt->Vpoint->y);
					zz.push_back(pnt->Vpoint->z);
					face.push_back(inode);
					pnt = pnt->nextcorner;
				}
				if (f->num_corners > 2) save_vtk_face(xx, yy, zz, face, iface);
									
				#pragma endregion
			};

			inline void plot_grid_cell_faces(size_t icell)
			{
				#pragma region Plot grid cell:
				std::vector<double> x;
				std::vector<double> y;
				std::vector<double> z;
				std::vector<std::vector<size_t> > faces;

				size_t i, j, k;
				get_cell_indices(icell, i, j, k);
				double xo(_xo + i * _s), yo(_yo + j * _s), zo(_zo + k * _s);
					
				// first point
				double x1(xo), y1(yo), z1(zo);
				x.push_back(x1); y.push_back(y1); z.push_back(z1);
					
				// second point
				double x2(xo + _s), y2(yo), z2(zo);
				x.push_back(x2); y.push_back(y2); z.push_back(z2);
				
				// third point
				double x3(xo + _s), y3(yo + _s), z3(zo);
				x.push_back(x3); y.push_back(y3); z.push_back(z3);
						
				// fourth point
				double x4(xo), y4(yo + _s), z4(zo);
			    x.push_back(x4); y.push_back(y4); z.push_back(z4);
								
				// fifth point				
				double x5(xo), y5(yo), z5(zo + _s);
				x.push_back(x5); y.push_back(y5); z.push_back(z5);

				// sixth point
				double x6(xo + _s), y6(yo), z6(zo + _s);
				x.push_back(x6); y.push_back(y6); z.push_back(z6);
				
				// seventh point
				double x7(xo + _s), y7(yo + _s), z7(zo + _s);
				x.push_back(x7); y.push_back(y7); z.push_back(z7);
						
				// eighth point
				double x8(xo), y8(yo + _s), z8(zo + _s);
				x.push_back(x8); y.push_back(y8); z.push_back(z8);

				size_t i1(0), i2(1), i3(2), i4(3), i5(4), i6(5), i7(6), i8(7);

				std::vector<size_t> f1(4);
				f1[0] = i1; f1[1] = i2; f1[2] = i3; f1[3] = i4;
				faces.push_back(f1);

				std::vector<size_t> f2(4);
				f2[0] = i1; f2[1] = i4; f2[2] = i8; f2[3] = i5;
				faces.push_back(f2);

				std::vector<size_t> f3(4);
				f3[0] = i1; f3[1] = i5; f3[2] = i6; f3[3] = i2;
				faces.push_back(f3);

				std::vector<size_t> f4(4);
				f4[0] = i7; f4[1] = i8; f4[2] = i4; f4[3] = i3;
				faces.push_back(f4);

				std::vector<size_t> f5(4);
				f5[0] = i7; f5[1] = i3; f5[2] = i2; f5[3] = i6;
				faces.push_back(f5);

				std::vector<size_t> f6(4);
				f6[0] = i7; f6[1] = i6; f6[2] = i5; f6[3] = i8;
				faces.push_back(f6);

				save_vtk_polygons(x, y, z, faces);
				#pragma endregion
			};

			
			inline bool edge_cutting_face(FacePoint* p, PolyFace* f)
			{
				#pragma region Edge in Face Check:
				double xo(p->Vpoint->x), yo(p->Vpoint->y), zo(p->Vpoint->z);
				double x1(p->nextcorner->Vpoint->x), y1(p->nextcorner->Vpoint->y), z1(p->nextcorner->Vpoint->z);
				FacePoint* pnt(f->corner);
				FacePoint* pntp;
				FacePoint* fpnt(pnt);
				while (true)
				{
					pntp = pnt->nextcorner;
					if (fabs(x1 - pnt->Vpoint->x) < 1E-10 && fabs(y1 - pnt->Vpoint->y) < 1E-10 && fabs(z1 - pnt->Vpoint->z) < 1E-10 &&
						fabs(xo - pntp->Vpoint->x) < 1E-10 && fabs(yo - pntp->Vpoint->y) < 1E-10 && fabs(zo - pntp->Vpoint->z) < 1E-10)
					{
						return false; // edge bounds face
					}
					pnt = pntp;
					if (pnt == fpnt) break;
				}
				return true;
				#pragma endregion
			};

			// output
			inline bool segment_crossing_face(double px, double py, double pz, double qx, double qy, double qz, PolyFace* f, 
				                              double &xx, double &yy, double &zz)
			{
				#pragma region Segment Crossing Face Check:
				double xo(f->corner->Vpoint->x);
				double yo(f->corner->Vpoint->y);
				double zo(f->corner->Vpoint->z);
				double nx(f->nx), ny(f->ny), nz(f->nz);

				double ax(xo - px), ay(yo - py), az(zo - pz);
				double bx(qx - px), by(qy - py), bz(qz - pz);

				double dot_a = ax * nx + ay * ny + az * nz;
				double dot_b = bx * nx + by * ny + bz * nz;

				if (fabs(dot_a) < 1E-10 || fabs(dot_a) > fabs(dot_b)) return false;

				double u, v;

				u = dot_a / dot_b;
				if (u < 1E-10 || u > 1.0 - 1E-10) return false;

				xx = px + u * (qx - px);
				yy = py + u * (qy - py);
				zz = pz + u * (qz - pz);				
				
				// check of intersection point lies in that face
				double cx, cy, cz;
				FacePoint* fpnt = f->corner;
				FacePoint* p1 = fpnt->nextcorner;
				FacePoint* p2;
				while (true)
				{
					p2 = p1->nextcorner;
					if (p2 == fpnt) break;

					ax = (p1->Vpoint->x - xo);
					ay = (p1->Vpoint->y - yo);
					az = (p1->Vpoint->z - zo);

					bx = (p2->Vpoint->x - xo);
					by = (p2->Vpoint->y - yo);
					bz = (p2->Vpoint->z - zo);

					cx = (xx - xo);
					cy = (yy - yo);
					cz = (zz - zo);

					double aa = ax * ax + ay * ay + az * az;
					double bb = bx * bx + by * by + az * az;
					double ab = ax * bx + ay * by + az * bz;
					double ac = ax * cx + ay * cy + az * cz;
					double bc = bx * cx + by * cy + bz * cz;

					double inv_den = 1.0 / (aa * bb - ab * ab);
					u = (bb * ac - ab * bc) * inv_den;
					v = (aa * bc - ab * ac) * inv_den;

					if (u > 0.0 && v > 0.0 && u + v < 1.0) return true;

					p1 = p2;
				}
				return false;
				#pragma endregion
			};

			void save_vcells(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
					 std::vector<std::vector<size_t> > &faces, std::vector<size_t> &elements)
			{
				#pragma region Save output Voronoi Cells:
				std::fstream file("VCells.dat", std::ios::out);
				// Points
				size_t num_points(x.size());
				size_t num_elements(elements.size() - 1);
				file << num_points << " " << num_elements << std::endl;           // number of points
				for (size_t i = 0; i < num_points; i++)
				{
					double xx(x[i]), yy(y[i]), zz(z[i]);
					scale_output(xx, yy, zz);
					file << xx << " " << yy << " " << zz << std::endl;           // number of points		
				}

				// VCells
				for (size_t ie = 0; ie < num_elements; ie++)
				{
					size_t num_faces(elements[ie + 1] - elements[ie]);	
					file << num_faces;
					for (size_t iface = elements[ie]; iface < elements[ie + 1]; iface++)
					{
						size_t num_face_corners(faces[iface].size());
						file << " " << num_face_corners;
						for (size_t i = 0; i < num_face_corners; i++) file << " " <<  faces[iface][i];
					}
					file << std::endl;
				}
				#pragma endregion
			};

			void save_face_normals(std::vector<double> &fn)
			{
				#pragma region Save output normals:
				std::fstream file("normals.dat", std::ios::out);
				// Points
				size_t num_faces(fn.size() / 3); size_t ii(0);
				for (size_t i = 0; i < num_faces; i++)
				{
					file << fn[ii]; ii++;
					file << " " << fn[ii]; ii++;
					file << " " << fn[ii]; ii++;
					file << std::endl;           // number of points		
				}
				#pragma endregion
			};

			void save_vtk_polygons(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
					       std::vector<std::vector<size_t> > &faces)
			{
				#pragma region Save vtk Polygons:
				std::fstream file("polygons.vtk", std::ios::out);
				file << "# vtk DataFile Version 3.0" << std::endl; // header
				file << "# vtk output" << std::endl;               // title
				file << "ASCII" << std::endl;                      // data type
				file << "DATASET POLYDATA" << std::endl;           // Topology
				// Points
				size_t num_points(x.size());
				file << "POINTS " << num_points << " float" << std::endl;           // number of points
				for (size_t i = 0; i < num_points; i++)
				{
					file << x[i] << " " << y[i] << " " << z[i] << std::endl;           // number of points		
				}
				// Polygons
				size_t num_polygons(faces.size());
				size_t num_data(0);
				num_data += num_polygons;
				//num_data += 1;
				for (size_t ie = 0; ie < num_polygons; ie++) 
				{
					//if (ie != 212) continue;
					num_data += faces[ie].size();
				}

				file << "POLYGONS " << num_polygons << " " << num_data << std::endl;           // number of points
				//file << "POLYGONS " << 1 << " " << num_data << std::endl;                         // single face only - debug
				for (size_t ie = 0; ie < num_polygons; ie++) 
				{
					//if (ie != 212) continue;
					size_t num_poly_points(faces[ie].size());
					file << num_poly_points;
					for (size_t i = 0; i < num_poly_points; i++) file << " " <<  faces[ie][i];
					file << std::endl;
				}
				#pragma endregion
			};

			void save_ply_mesh(std::string file_name, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
					   std::vector<std::vector<size_t> > &faces)
			{
				#pragma region Save PLY mesh:
				std::fstream file(file_name.c_str(), std::ios::out);
				file << "ply" << std::endl; // header
				file << "comment MeshingGenie generated" << std::endl;            
				file << "element vertex" << x.size() << std::endl;
				file << "property float x" << std::endl;          
				file << "property float y" << std::endl;          
				file << "property float z" << std::endl;          
				file << "element face" << faces.size() << std::endl;
				file << "property list uchar int vertex_indices" << std::endl;            
				file << "end_header" << std::endl;            

				// Points
				size_t num_points(x.size());
				for (size_t i = 0; i < num_points; i++)
				{
					file << x[i] << " " << y[i] << " " << z[i] << std::endl;           // number of points		
				}
				// faces
				size_t num_faces(faces.size());
				for (size_t i = 0; i < num_points; i++)
				{
					size_t num_corners(faces[i].size());
					file << num_corners;
					for (size_t j = 0; j < num_corners; j++)
					{
						file << " " << faces[i][j];
					}						
					file << std::endl;           // number of points		
				}
				#pragma endregion
			};

			void save_vtk_face(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
	                           std::vector<size_t> &face, size_t iface)
			{
				#pragma region Save vtk Face:

				std::stringstream ss;
				ss << "face_" << iface <<".vtk";
				std::fstream file(ss.str().c_str(), std::ios::out);
				file << "# vtk DataFile Version 3.0" << std::endl; // header
				file << "# vtk output" << std::endl;               // title
				file << "ASCII" << std::endl;                      // data type
				file << "DATASET POLYDATA" << std::endl;           // Topology
				// Points
				size_t num_points(x.size());
				file << "POINTS " << num_points << " float" << std::endl;           // number of points
				for (size_t i = 0; i < num_points; i++)
				{
					file << x[i] << " " << y[i] << " " << z[i] << std::endl;           // number of points		
				}
				// Polygons
				
				size_t num_data(0);
				num_data += 1;
				num_data += face.size();				

				file << "POLYGONS " << 1 << " " << num_data << std::endl;
				size_t num_poly_points(face.size());
				file << num_poly_points;
				for (size_t i = 0; i < num_poly_points; i++) file << " " <<  face[i];
				file << std::endl;			
				#pragma endregion
			};

			void plot_boundary_cells()
			{
				#pragma region Plot Boundary Cells:
				std::vector<double> x;
				std::vector<double> y;
				std::vector<double> z;
				std::vector<std::vector<size_t> > faces;

				size_t num_pnts(0);
				for (size_t icell = 0; icell < _n; icell++)
				{
					if (_invalid_cells[icell]) continue;

					if (_boundary_cells[icell] || _external_cells[icell]) continue;

					size_t i, j, k;
					get_cell_indices(icell, i, j, k);
					double xo(_xo + i * _s), yo(_yo + j * _s), zo(_zo + k * _s);
					
					// first point
					double x1(xo), y1(yo), z1(zo); size_t i1(num_pnts);
					for (size_t ipnt = 0; ipnt < num_pnts; ipnt++)
					{
						if (fabs(x[ipnt] - x1) < 1E-10 && fabs(y[ipnt] - y1) < 1E-10 && fabs(z[ipnt] - z1) < 1E-10)
						{
							i1 = ipnt; break;
						}
					}
					if (i1 == num_pnts)
					{
						x.push_back(x1);
						y.push_back(y1);
						z.push_back(z1);
						num_pnts++;
					}

					// second point
					double x2(xo + _s), y2(yo), z2(zo); size_t i2(num_pnts);
					for (size_t ipnt = 0; ipnt < num_pnts; ipnt++)
					{
						if (fabs(x[ipnt] - x2) < 1E-10 && fabs(y[ipnt] - y2) < 1E-10 && fabs(z[ipnt] - z2) < 1E-10)
						{
							i2 = ipnt; break;
						}
					}
					if (i2 == num_pnts)
					{
						x.push_back(x2);
						y.push_back(y2);
						z.push_back(z2);
						num_pnts++;
					}

					// third point
					double x3(xo + _s), y3(yo + _s), z3(zo); size_t i3(num_pnts);
					for (size_t ipnt = 0; ipnt < num_pnts; ipnt++)
					{
						if (fabs(x[ipnt] - x3) < 1E-10 && fabs(y[ipnt] - y3) < 1E-10 && fabs(z[ipnt] - z3) < 1E-10)
						{
							i3 = ipnt; break;
						}
					}
					if (i3 == num_pnts)
					{
						x.push_back(x3);
						y.push_back(y3);
						z.push_back(z3);
						num_pnts++;
					}

					// fourth point
					double x4(xo), y4(yo + _s), z4(zo); size_t i4(num_pnts);
					for (size_t ipnt = 0; ipnt < num_pnts; ipnt++)
					{
						if (fabs(x[ipnt] - x4) < 1E-10 && fabs(y[ipnt] - y4) < 1E-10 && fabs(z[ipnt] - z4) < 1E-10)
						{
							i4 = ipnt; break;
						}
					}
					if (i4 == num_pnts)
					{
						x.push_back(x4);
						y.push_back(y4);
						z.push_back(z4);
						num_pnts++;
					}

					// fifth point
					double x5(xo), y5(yo), z5(zo + _s); size_t i5(num_pnts);
					for (size_t ipnt = 0; ipnt < num_pnts; ipnt++)
					{
						if (fabs(x[ipnt] - x5) < 1E-10 && fabs(y[ipnt] - y5) < 1E-10 && fabs(z[ipnt] - z5) < 1E-10)
						{
							i5 = ipnt; break;
						}
					}
					if (i5 == num_pnts)
					{
						x.push_back(x5);
						y.push_back(y5);
						z.push_back(z5);
						num_pnts++;
					}

					// sixth point
					double x6(xo + _s), y6(yo), z6(zo + _s); size_t i6(num_pnts);
					for (size_t ipnt = 0; ipnt < num_pnts; ipnt++)
					{
						if (fabs(x[ipnt] - x6) < 1E-10 && fabs(y[ipnt] - y6) < 1E-10 && fabs(z[ipnt] - z6) < 1E-10)
						{
							i6 = ipnt; break;
						}
					}
					if (i6 == num_pnts)
					{
						x.push_back(x6);
						y.push_back(y6);
						z.push_back(z6);
						num_pnts++;
					}

					// seventh point
					double x7(xo + _s), y7(yo + _s), z7(zo + _s); size_t i7(num_pnts);
					for (size_t ipnt = 0; ipnt < num_pnts; ipnt++)
					{
						if (fabs(x[ipnt] - x7) < 1E-10 && fabs(y[ipnt] - y7) < 1E-10 && fabs(z[ipnt] - z7) < 1E-10)
						{
							i7 = ipnt; break;
						}
					}
					if (i7 == num_pnts)
					{
						x.push_back(x7);
						y.push_back(y7);
						z.push_back(z7);
						num_pnts++;
					}

					// eighth point
					double x8(xo), y8(yo + _s), z8(zo + _s); size_t i8(num_pnts);
					for (size_t ipnt = 0; ipnt < num_pnts; ipnt++)
					{
						if (fabs(x[ipnt] - x8) < 1E-10 && fabs(y[ipnt] - y8) < 1E-10 && fabs(z[ipnt] - z8) < 1E-10)
						{
							i8 = ipnt; break;
						}
					}
					if (i8 == num_pnts)
					{
						x.push_back(x8);
						y.push_back(y8);
						z.push_back(z8);
						num_pnts++;
					}

					std::vector<size_t> f1(4);
					f1[0] = i1; f1[1] = i2; f1[2] = i3; f1[3] = i4;
					faces.push_back(f1);

					std::vector<size_t> f2(4);
					f2[0] = i1; f2[1] = i4; f2[2] = i8; f2[3] = i5;
					faces.push_back(f2);

					std::vector<size_t> f3(4);
					f3[0] = i1; f3[1] = i5; f3[2] = i6; f3[3] = i2;
					faces.push_back(f3);

					std::vector<size_t> f4(4);
					f4[0] = i7; f4[1] = i8; f4[2] = i4; f4[3] = i3;
					faces.push_back(f4);

					std::vector<size_t> f5(4);
					f5[0] = i7; f5[1] = i3; f5[2] = i2; f5[3] = i6;
					faces.push_back(f5);

					std::vector<size_t> f6(4);
					f6[0] = i7; f6[1] = i6; f6[2] = i5; f6[3] = i8;
					faces.push_back(f6);
				}

				save_vtk_polygons(x, y, z, faces);

				#pragma endregion
			}

			void test_voronoi_cells(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
						std::vector<std::vector<size_t> > &faces, std::vector<double> faces_normals)
			{
				#pragma region Test Voronoi faces:
				std::cout<< "*** Testing Voronoi faces ... " << std::endl;

				// Polygons
				size_t num_polygons(faces.size());

				double x1, y1, z1;
				double nx, ny, nz;				
				for (size_t ie = 0; ie < num_polygons; ie++) 
				{
					size_t num_poly_points(faces[ie].size());
					if (num_poly_points < 3)
					{
						std::cout<<" ** Error: Face " << ie << " is degenerated - number of corners = " << num_poly_points << std::endl;
						continue;
					}
					x1 = x[faces[ie][0]]; y1 = y[faces[ie][0]]; z1 = z[faces[ie][0]];

					nx = faces_normals[ie * 3];
					ny = faces_normals[ie * 3 + 1];
					nz = faces_normals[ie * 3 + 2];

					if (num_poly_points == 3) continue;

					for (size_t i = 1; i < num_poly_points; i++) 
					{
						double d  = nx * (x[faces[ie][i]] - x1) + ny * (y[faces[ie][i]] - y1) + nz * (z[faces[ie][i]] - z1);
						if (d > 0.05 * _dm)
						{					
							std::cout << "    ^ Error: Face " << ie << " is not coplanar, deviataion = " << d << "! " << std::endl;
							save_vtk_face(x, y, z, faces[ie], ie);
							break;						
						}
					}
				}
				std::cout<< "*** Testing Voronoi faces ... done!" << std::endl;
				#pragma endregion
			};

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Sampling Methods ::::::::::::::::
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			int execute(size_t fixed_seed);

			void generate_background_grid();

			void identify_external_cell();

			int sprinkle_points_along_sharp_edges(double d);

			int sample_edge(double xo, double yo, double zo, double xn, double yn, double zn, double d,
	                             std::vector<double> &ex, std::vector<double> &ey, std::vector<double> &ez);


			// used a lot
			inline double generate_a_random_number();

			inline double distance_squared(double dx, double dy, double dz);

			inline void get_cell_index(size_t i, size_t j, size_t k, size_t &icell){icell = _njk * i + _nk* j + k;};

			inline void get_cell_indices(size_t icell, size_t &i, size_t &j, size_t &k){i = icell / _njk; j = (icell - i * _njk) / _nk; k = icell - i * _njk - j * _nk;};

			inline void get_neighbor_indices(size_t index, size_t io, size_t jo, size_t ko, size_t &ii, size_t &jj, size_t &kk);

			inline void get_active_children(size_t i, size_t j, size_t k, size_t icell, 
											bool &c1, bool &c2, bool &c3, bool &c4, bool &c5, bool &c6, bool &c7, bool &c8);

			void use_grid_level(size_t refLevel);

			void test_output();

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Voronoi Meshing Methods ::::::::::::::::
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			void generate_CVM();

			inline bool construct_internal_boundary_voronoi_cell(size_t icell, size_t &v_cell_index);

			inline bool construct_reflex_voronoi_cell(size_t icell, size_t &v_cell_index);
			
			inline bool construct_regular_voronoi_cell(size_t icell, size_t &v_cell_index);

			inline bool trim_face(PolyFace* f, double x1, double y1, double z1,
								               double nx, double ny, double nz, 
								               double &px, double &py, double &pz, 
								               double &qx, double &qy, double &qz);

			inline bool get_trimming_points(PolyFace* f, double x1, double y1, double z1,
								               double nx, double ny, double nz, 
								               double &px, double &py, double &pz, 
								               double &qx, double &qy, double &qz);

			inline void weld_Voronoi_Point(VoronoiPoint* &vpoint, Polyhedron* poly, double tol);

			inline Polyhedron* remove_external_polys(Polyhedron* poly);

			inline void remove_degenerate_faces(Polyhedron* poly);

			inline void weld_poly_corners(Polyhedron* poly);

			inline bool add_face_to_poly(Polyhedron* poly, PolyFace* face);
			
			inline Polyhedron* create_a_cell(size_t icell);

			inline bool is_external_poly(Polyhedron* poly);

			inline void delete_poly_faces(Polyhedron* poly);

			inline Polyhedron* copy_polyhedron(Polyhedron* poly);	

			inline PolyFace* copy_polyface(PolyFace* f, bool invert_normals);

			inline PolyFace* construct_super_triangle(double xmid, double ymid, double zmid, double nx, double ny, double nz,
				                                      double xc, double yc, double zc, Point* q);

			inline bool point_in_boundary_face(size_t iface, double px, double py, double pz);

			inline void delete_face_corners(PolyFace* f);

			inline PolyFace* create_boundary_face(size_t iface);

			inline PolyFace* create_internal_face(size_t iface);
			
			void sample_face(size_t iface);

			void detect_internal_cracks_on_internal_faces();

		private:
		
			size_t _num_neighbors;
			int* _neighbors;

			size_t _n, _ni, _nj, _nk, _njk, _no, _pool_size;
			double _dm, _dm_squared, _xo, _yo, _zo, _s, _ss;
			double _tol;

			Point** _cell_points;
			bool* _invalid_cells;
			size_t* _active_pool_i;
			size_t* _active_pool_j;
			size_t* _active_pool_k;
			size_t* _active_pool_icell;

			bool* _boundary_cells;
			bool* _external_cells;

			Polyhedron** _v_cells;

			size_t _num_valid, _num_inserted_points;

			// Constants
			double _fixed_seed;
			size_t _MY_RAND_MAX;
			double _SQRT_3;
			double _SQRT_3_INV;
			size_t _TWO_POW_REF_LEVEL;
			double _RF;

			VoronoiPoint** _cell_vpoints;

			std::vector<double> _xb;
			std::vector<double> _yb;
			std::vector<double> _zb;
			std::vector< std::vector<size_t> > _bfaces;
			std::vector< std::vector<size_t> > _b_node_faces;
			std::vector<double> _bfaces_nx;
			std::vector<double> _bfaces_ny;
			std::vector<double> _bfaces_nz;

			// internal faces
			std::vector< std::vector<double> > _xf;
			std::vector< std::vector<double> > _yf;
			std::vector< std::vector<double> > _zf;
			std::vector< std::vector< std::vector<double> > > _xfc;
			std::vector< std::vector< std::vector<double> > > _yfc;
			std::vector< std::vector< std::vector<double> > > _zfc;
			std::vector< double > _nfx;
			std::vector< double > _nfy;
			std::vector< double > _nfz;

			std::vector<bool> _nonconvex_bfaces;

			std::map<size_t, std::set<size_t> > _cell_faces;
			std::map<size_t, std::set<size_t> > _cell_internal_faces;
			std::map<size_t, std::set<size_t> >::iterator _cell_faces_iter;
			std::vector<size_t> _sharp_edges;
			std::vector<size_t> _sharp_corners;
		
			double _box_xmin, _box_xmax, _box_ymin, _box_ymax, _box_zmin, _box_zmax;

			double _DX, _DY, _DZ, _SF, _SF_inv;

};
                                
#endif	
