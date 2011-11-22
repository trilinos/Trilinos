// 09/01/2011 5:41 pm

// R1.0: Non-convex domains (sampling + Voronoi meshing + output and testing planar faces + non-convex domains + internal faces)

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

// R5.0

#ifndef MESHING_GENIE_2D_H
#define MESHING_GENIE_2D_H

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

class MeshingGenie_2d   
{
	protected:
		struct PointEdge;

		struct Point 
		{
			/// Default constructor does nothing (for performance).
			Point() {}
			double x; double y;
			Point*     next; // A point could be a list of points
			PointEdge* edges;
			size_t num_edges;
			bool       done; // if true then edge connectivity is constructed
			bool       on_boundary; 
		};

		struct PointEdge 
		{
			PointEdge(){}			
			Point* edge_other_end;
			size_t edge_id;  // dirctionality 1 bit then edge index + 1
			PointEdge* next; // pointer to the next constrained edge
		};

		struct PolyPoint 
		{
			/// Default constructor does nothing (for performance).
			PolyPoint() {}
			double x;
			double y;
			double u;
			PolyPoint* next;
			PolyPoint* prev;
			Point* center_L;
			Point* center_R;
		};

		struct Poly
		{
			/// Default constructor does nothing (for performance).
			Poly() {}
			size_t icell;
			double area;
			PolyPoint* pp;		
		};

		struct LoopPoint 
		{
			/// Default constructor does nothing (for performance).
			LoopPoint() {}
			Point* pp;
			LoopPoint* next;
			LoopPoint* prev;
			double  Lsq;
		};

    public:
        //! constructor
        MeshingGenie_2d(double dm, std::vector<double> &ExternalBoundaries,
						           std::vector< std::vector<double> > &Holes, 
						           std::vector< std::vector<double> > &Cracks,
								   int output, // 0 : point cloud - 1: Voronoi Cells - 2: CDT - 3: All-quad
								   bool plot_results);
   

        //! Destructor                
        ~MeshingGenie_2d(){ };

		void use_fixed_seed(size_t fixed_seed){_fixed_seed = fixed_seed;};

        int execute();    

		void get_point_cloud(std::vector<double> &x, std::vector<double> &y);

		void get_CDT_Tessellation(std::vector<double> &x, std::vector<double> &y, std::vector< std::vector<size_t> > &elements);

		void get_Voronoi_Tessellation(std::vector<double> &x, std::vector<double> &y, 
			                          std::vector< std::vector<size_t> > &elements, double min_edge_length, 
									  bool generate_side_sets, bool split_cracktipcells);

		void get_side_sets(std::vector< std::vector<size_t> > &side_sets){
			side_sets = _side_sets;
		};

	private:

        int apply_grid_based_method();

        int extract_boudary_edges();	

        int generate_backgorund_grid(double ds);

        int mark_cells_crossing_boundaries( );

		int sprinkle_points_along_boundaries(double d);

		int fill_neighbors_data();

		int refill_neighbors_data();

		int apply_maximal_poisson_disk_sampler();

		int clean_crowded_cells(bool eliminate_isolated_voids, bool speak_to_me);
		
		int eliminate_voids_via_polys();

		inline void pick_a_random_point_from_a_convex_hull(PolyPoint* pp, double &xx, double &yy);

		inline bool is_valid_internal_point(size_t icell, double xx, double yy);

		inline bool cut_polygon(PolyPoint* &pp, size_t &num_poly_pnts, Point** circles, size_t num_circles);

		inline bool split_polygon(PolyPoint* pp);

		int CDT(bool Voronoi_Cells);

        int plot_vertices(std::string file_name, bool plot_circles, bool plot_grid);

		int plot_cells(std::string file_name, bool blue_cells);

		int plot_cell(std::string file_name, size_t icell, PolyPoint* pp, size_t kcell);

		int plot_point(std::string file_name, size_t icell, LoopPoint* po);

		int plot_star(std::string file_name, double xo, double yo, size_t icell, Point** neighbors, bool* invalid, int num);

		int plot_tessellation(std::string file_name, bool Voroni);

		// IO Methods
		void print_message(std::string s)    {std::cout << "     * " << s << std::endl;};                
		void print_message_end(std::string s){std::cout << "       " << s << std::endl;};
		std::string double_to_string(const double &d)
		{
			std::ostringstream i; i << d; return i.str();
		};

		double string_to_double(const std::string &s)
		{			 
			std::istringstream i(s); double x;
			if (!(i >> x)) return 0; return x;			
		};

		inline double generate_a_random_number()
        {           
            size_t num(rand());
            num = num << 15;
            num += rand(); // another random number
            return static_cast<double>(num * 1.0 / _MY_RAND_MAX);
        };
		
		inline void get_cell_indices(size_t handle, size_t &i, size_t &j)
        {
            j = handle % (_nr - 1); i = (handle - j) / (_nr - 1);            
        };

        inline void get_cell_handle(size_t i, size_t j, size_t &handle)
        {
            handle = i * (_nr - 1) + j;            
        };

		inline void get_cell_handle(double x, double y, size_t &handle)
        {						
			size_t i = size_t(floor((x - _bg_xo) / _bg_s));
			if (_bg_xo + i * _bg_s - x > _bg_s - 1E-10) i++;
			size_t j = size_t(floor((y - _bg_yo) / _bg_s));
			if (_bg_yo + j * _bg_s - y > _bg_s - 1E-10) j++;
            handle = i * (_nr - 1) + j;
        };

		int input_translator(std::vector<double> &ExtBound, std::vector< std::vector<double> > &Cracks)
		{			
			#pragma region check for intersections between cracks and external boundaries
			std::vector<double> new_boundaries;    
			std::vector<double> pxv;
			std::vector<double> pyv;
			std::vector<double> pdv;
			double xo, yo, xn, yn;
			int num_cracks = Cracks.size();
			int iimax = ExtBound.size();
			int num_boundary_edges(ExtBound.size()/2);
			int ii(0); double px, py;
			for (int ied = 0; ied < num_boundary_edges; ied++)
			{        							
				#pragma region Working On External Boundaries:
				pxv.clear(); pyv.clear();
				xo = ExtBound[ii]; ii++;
				yo = ExtBound[ii]; ii++;
				if (ii == iimax) ii = 0;
				xn = ExtBound[ii]; ii++;
				yn = ExtBound[ii]; ii++;

				new_boundaries.push_back(xo);            
				new_boundaries.push_back(yo);
		
				for (int ic = 0; ic < num_cracks; ic++)
				{
					int num_crack_edges(Cracks[ic].size() / 2 - 1);
					int jj = 0;
					for (int jed = 0; jed < num_crack_edges; jed++)
					{
						double x1 = Cracks[ic][jj]; jj++;
						double y1 = Cracks[ic][jj]; jj++;

						double x2 = Cracks[ic][jj]; jj++;
						double y2 = Cracks[ic][jj]; jj++;
                
						if (intersecting_segments(xo, yo, xn, yn, x1, y1, x2, y2, px, py))
						{
							pxv.push_back(px);                            
							pyv.push_back(py);
						}
						jj -= 2;
					}            
				}
				int num_intersections(pxv.size());
				pdv.clear(); pdv.resize(num_intersections);
				for (int i = 0; i < num_intersections; i++)
				{
					pdv[i] = distance_squared(xo - pxv[i], yo - pyv[i]);
				}

				// sort intersections based on pdv
				// TODO:: use more effecient sorting method!!!
				for (int i = 0; i < num_intersections; i++)
				{
					for (int j = i + 1; j < num_intersections; j++)
					{
						if (pdv[i] > pdv[j])
						{
							double tmp = pdv[i]; pdv[i] = pdv[j]; pdv[j] = tmp;
							tmp = pxv[i]; pxv[i] = pxv[j]; pxv[j] = tmp;
							tmp = pyv[i]; pyv[i] = pyv[j]; pyv[j] = tmp;
						}
					}
				}

				// add intersection points to external boundaries
				for (int i = 0; i < num_intersections; i++)
				{
					if (i >= 1)
					{
						if (fabs(pxv[i]- pxv[i - 1]) < 1E-10 && fabs(pyv[i]- pyv[i - 1]) < 1E-10) continue;			
					}
					new_boundaries.push_back(pxv[i]);
					new_boundaries.push_back(pyv[i]);            
				}

				ii -= 2;
				#pragma endregion
			}
			ExtBound = new_boundaries;

			std::vector< std::vector<double> > new_cracks(num_cracks);    
			for (int icrack = 0; icrack < num_cracks; icrack++)
			{        
				#pragma region Working On an Internal crack:

				ii = 0; int num_crack_edges(Cracks[icrack].size() / 2 - 1);

				for (int ied = 0; ied < num_crack_edges; ied++)
				{
					pxv.clear(); pyv.clear();
					xo = Cracks[icrack][ii]; ii++;
					yo = Cracks[icrack][ii]; ii++;        
					xn = Cracks[icrack][ii]; ii++;
					yn = Cracks[icrack][ii]; ii++;

					new_cracks[icrack].push_back(xo); new_cracks[icrack].push_back(yo);
					for (int ic = 0; ic < num_cracks; ic++)
					{
						if (ic == icrack) continue;

						int num_crack_edges(Cracks[ic].size() / 2 - 1);
						int jj = 0;
						for (int jed = 0; jed < num_crack_edges; jed++)
						{
							double x1 = Cracks[ic][jj]; jj++;
							double y1 = Cracks[ic][jj]; jj++;

							double x2 = Cracks[ic][jj]; jj++;
							double y2 = Cracks[ic][jj]; jj++;							
							if (intersecting_segments(xo, yo, xn, yn, x1, y1, x2, y2, px, py))
							{								
								pxv.push_back(px);                            
								pyv.push_back(py);								
							}
							jj -= 2;
						}            
					}
					int num_intersections(pxv.size());
					pdv.clear(); pdv.resize(num_intersections);
					for (int i = 0; i < num_intersections; i++)
					{
						pdv[i] = distance_squared(xo - pxv[i], yo - pyv[i]);
					}

					// sort intersections based on pdv
					// TODO:: use more effecient sorting method!!!
					for (int i = 0; i < num_intersections; i++)
					{
						for (int j = i + 1; j < num_intersections; j++)
						{
							if (pdv[i] > pdv[j])
							{
								double tmp = pdv[i]; pdv[i] = pdv[j]; pdv[j] = tmp;
								tmp = pxv[i]; pxv[i] = pxv[j]; pxv[j] = tmp;
								tmp = pyv[i]; pyv[i] = pyv[j]; pyv[j] = tmp;
							}
						}
					}

					// add intersection points to external boundaries
					for (int i = 0; i < num_intersections; i++)
					{
						if (i >= 1)
						{
							if (fabs(pxv[i]- pxv[i - 1]) < 1E-10 && fabs(pyv[i]- pyv[i - 1]) < 1E-10) continue;			
						}
						new_cracks[icrack].push_back(pxv[i]);
						new_cracks[icrack].push_back(pyv[i]);            
					}

					ii -= 2;
				}
				new_cracks[icrack].push_back(xn); new_cracks[icrack].push_back(yn);        				
				#pragma endregion
			}

			for (int icrack = 0; icrack < num_cracks; icrack++)
			{
				Cracks[icrack] = new_cracks[icrack];
			}
			return 0;
			#pragma endregion
		};

		inline void add_constraint_edge(Point* pi, Point* pj)
		{
			#pragma region Add Contained directional edge:
			if (pi->done || pj->done) return;

			if (pi->edges == 0) 
			{
				pi->edges = new PointEdge();
				pi->edges->edge_other_end = pj;
				pi->edges->edge_id = 0; // unclassified edge
				pi->edges->next = 0; pi->num_edges++;			
			}
			else
			{
				PointEdge* ed = pi->edges;
				while (true)
				{
					if (ed->edge_other_end == pj) return; // edge already exists
					if (ed->next == 0 || ed->next == pi->edges) break;
					ed = ed->next;
				}
				ed->next = new PointEdge();
				ed->next->edge_other_end = pj;
				ed->next->edge_id = 0; // unclassified edge
				ed->next->next = 0; pi->num_edges++;
			}

			if (pj->edges == 0) 
			{
				pj->edges = new PointEdge();
				pj->edges->edge_other_end = pi;
				pj->edges->edge_id = 0; // unclassified edge
				pj->edges->next = 0; pj->num_edges++;
			}
			else
			{
				PointEdge* ed = pj->edges;
				while (true)
				{
					if (ed->edge_other_end == pi) return; // edge already exists
					if (ed->next == 0 || ed->next == pj->edges) break;
					ed = ed->next;
				}
				ed->next = new PointEdge();
				ed->next->edge_other_end = pi;
				ed->next->edge_id = 0; // unclassified edge
				ed->next->next = 0; pj->num_edges++;
			}
			#pragma endregion
		};
		
		inline void add_constraint_edge(Point* pi, Point* pj, size_t iedge)
		{
			#pragma region Add Contained directional edge:
			
			if (pi->done || pj->done) return;

			if (pi->edges == 0) 
			{
				size_t forward_id = generate_edge_id(iedge, true);
				pi->edges = new PointEdge();
				pi->edges->edge_other_end = pj;
				pi->edges->edge_id = forward_id;
				pi->edges->next = 0; pi->num_edges++;
			}
			else
			{
				size_t forward_id = generate_edge_id(iedge, true);
				PointEdge* ed = pi->edges;
				while (true)
				{
					if (ed->edge_other_end == pj) return; // edge already exists
					if (ed->next == 0) break;
					ed = ed->next;
				}
				ed->next = new PointEdge();
				ed->next->edge_other_end = pj;
				ed->next->edge_id = forward_id;
				ed->next->next = 0; pi->num_edges++;
			}
			
			if (pj->edges == 0) 
			{
				size_t backward_id = generate_edge_id(iedge, false);
				pj->edges = new PointEdge();
				pj->edges->edge_other_end = pi;
				pj->edges->edge_id = backward_id;
				pj->edges->next = 0; pj->num_edges++;
			}
			else
			{
				size_t backward_id = generate_edge_id(iedge, false);
				PointEdge* ed = pj->edges;
				while (true)
				{
					if (ed->edge_other_end == pi) return; // edge already exists
					if (ed->next == 0) break;
					ed = ed->next;
				}
				ed->next = new PointEdge();
				ed->next->edge_other_end = pi;
				ed->next->edge_id = backward_id;
				ed->next->next = 0; pj->num_edges;
			}
			#pragma endregion
		};

		// false means not sure
		inline bool connected_points(Point* pi, Point* pj)
		{		
			#pragma region Check if two points are connected with an edge:
			if (pi->edges == 0 || pj->edges == 0) return false;
			
			PointEdge* ed = pi->edges;
			while (true)
			{
				if (ed->edge_other_end == pj) return true; // a constrained edge
				if (ed->next == 0 || ed->next == pi->edges) return false;
				ed = ed->next;
			}
			return false;
			#pragma endregion
		};

		inline bool oriented_edge_exists(Point* pi, Point* pj, bool &forward)
		{		
			#pragma region Check if two points are connected with an edge:
			if (pi->edges == 0 || pj->edges == 0) return false;
			
			PointEdge* ed = pi->edges;
			size_t edge_index;
			while (true)
			{
				if (ed->edge_other_end == pj) 
				{
					if (ed->edge_id == 0) return false;
					edge_index = (ed->edge_id - 1) >> 1;
					if (edge_index >= _num_oriented_edges) return false;
					forward = bool((ed->edge_id - 1) & 1);
					return true; // a constrained edge
				}

				if (ed->next == 0 || ed->next == pi->edges) return false;
				ed = ed->next;
			}
			return false;
			#pragma endregion
		};

		inline bool boundary_edge_exists(Point* pi, Point* pj)
		{		
			#pragma region Check if two points are connected with an edge:
			if (pi->edges == 0 || pj->edges == 0) return false;
			
			PointEdge* ed = pi->edges;
			while (true)
			{
				if (ed->edge_other_end == pj) 
				{
					if (ed->edge_id == 0) return false;
					return true; // a non-oriented edge
				}
				if (ed->next == 0 || ed->next == pi->edges) return false;
				ed = ed->next;
			}
			return false;
			#pragma endregion
		};


		// false means not sure
		inline bool notconnected_points(Point* pi, Point* pj)
		{		
			#pragma region Check if two points are connected with an edge:
			if (pi->edges == 0 || pj->edges == 0) return false;
				
			bool connected = connected_points(pi, pj);
			if ((pi->done || pj->done) && !connected) return true;
			return false;
			#pragma endregion
		};

		inline bool is_boundary_edge(Point* pi, Point* pj, size_t& edge_id)
		{		
			#pragma region Check if two points are connected with a boundary edge:
			if (pi->edges == 0 || pj->edges == 0) return false;
			
			PointEdge* ed = pi->edges;
			while (true)
			{
				if (ed->edge_other_end == pj && ed->edge_id != 0)
				{
					edge_id = ed->edge_id;
					return true; // a constrained edge
				}
				if (ed->next == 0 || ed->next == pi->edges) return false;
				ed = ed->next;
			}
			return false;
			#pragma endregion
		};

		inline void invalidate_extrenal_cell(size_t icell, bool N, bool W, bool S, bool E, bool firstcell)
		{
			#pragma region Invalidate Extrenal Cell via propagation:
			
			if (!firstcell && _bad_cells[icell]) return;

			if (!_bad_cells[icell])
			{
				_bad_cells[icell] = true; _num_bad_cells++;
			}
			
			size_t i, j;
			get_cell_indices(icell, i, j);
			if (firstcell && S && W)
			{				
				#pragma region South West Corner:
				size_t ied = _cell_first_edge_map[icell];
				size_t jed = _cell_second_edge_map[icell];
				double x1 = _edges_x1[ied]; double y1 = _edges_y1[ied];
				double x2 = _edges_x2[ied]; double y2 = _edges_y2[ied];
				double x3 = _edges_x2[jed]; double y3 = _edges_y2[jed];
				double L12 = sqrt(distance_squared(x2 - x1, y2 - y1));
				x1 = x2 + _dm * (x1 - x2) / L12; 
				y1 = y2 + _dm * (y1 - y2) / L12;
				double L23 = sqrt(distance_squared(x2-x3, y2 - y3));
				x3 = x2 + _dm * (x3 - x2) / L23; 
				y3 = y2 + _dm * (y3 - y2) / L23;
				double xmid = 0.5 * (x1 + x3); double ymid = 0.5 * (y1 + y3);
				if (xmid < x2 && ymid < y2) return;			
				#pragma endregion
			}
			else if (firstcell && S && E)
			{				
				#pragma region South East Corner:
				size_t ied = _cell_first_edge_map[icell];
				size_t jed = _cell_second_edge_map[icell];
				double x1 = _edges_x1[ied]; double y1 = _edges_y1[ied];
				double x2 = _edges_x2[ied]; double y2 = _edges_y2[ied];
				double x3 = _edges_x2[jed]; double y3 = _edges_y2[jed];
				double L12 = sqrt(distance_squared(x2 - x1, y2 - y1));
				x1 = x2 + _dm * (x1 - x2) / L12; 
				y1 = y2 + _dm * (y1 - y2) / L12;
				double L23 = sqrt(distance_squared(x2-x3, y2 - y3));
				x3 = x2 + _dm * (x3 - x2) / L23; 
				y3 = y2 + _dm * (y3 - y2) / L23;
				double xmid = 0.5 * (x1 + x3); double ymid = 0.5 * (y1 + y3);
				if (xmid > x2 && ymid < y2) return;			
				#pragma endregion
			}
			else if (firstcell && N && E)
			{				
				#pragma region North East Corner:
				size_t ied = _cell_first_edge_map[icell];
				size_t jed = _cell_second_edge_map[icell];
				double x1 = _edges_x1[ied]; double y1 = _edges_y1[ied];
				double x2 = _edges_x2[ied]; double y2 = _edges_y2[ied];
				double x3 = _edges_x2[jed]; double y3 = _edges_y2[jed];
				double L12 = sqrt(distance_squared(x2 - x1, y2 - y1));
				x1 = x2 + _dm * (x1 - x2) / L12; 
				y1 = y2 + _dm * (y1 - y2) / L12;
				double L23 = sqrt(distance_squared(x2-x3, y2 - y3));
				x3 = x2 + _dm * (x3 - x2) / L23; 
				y3 = y2 + _dm * (y3 - y2) / L23;
				double xmid = 0.5 * (x1 + x3); double ymid = 0.5 * (y1 + y3);
				if (xmid > x2 && ymid > y2) return;			
				#pragma endregion
			}
			else if (firstcell && N && W)
			{				
				#pragma region North West Corner:
				size_t ied = _cell_first_edge_map[icell];
				size_t jed = _cell_second_edge_map[icell];
				double x1 = _edges_x1[ied]; double y1 = _edges_y1[ied];
				double x2 = _edges_x2[ied]; double y2 = _edges_y2[ied];
				double x3 = _edges_x2[jed]; double y3 = _edges_y2[jed];
				double L12 = sqrt(distance_squared(x2 - x1, y2 - y1));
				x1 = x2 + _dm * (x1 - x2) / L12; 
				y1 = y2 + _dm * (y1 - y2) / L12;
				double L23 = sqrt(distance_squared(x2-x3, y2 - y3));
				x3 = x2 + _dm * (x3 - x2) / L23; 
				y3 = y2 + _dm * (y3 - y2) / L23;
				double xmid = 0.5 * (x1 + x3); double ymid = 0.5 * (y1 + y3);
				if (xmid < x2 && ymid > y2) return;			
				#pragma endregion
			}
			
			if (i == 0) W = false;
			if (i == _imax) E = false;
			if (j == 0) S = false;
			if (j == _jmax) N = false;

			if (N && W && S)
			{
				invalidate_extrenal_cell(icell + _neighbors[2], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[6], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[3], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[7], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[4], N, W, S, E, false);
			}
			else if (W && S && E)
			{
				invalidate_extrenal_cell(icell + _neighbors[3], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[7], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[4], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[8], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[1], N, W, S, E, false);
			}
			else if (S && E && N)
			{
				invalidate_extrenal_cell(icell + _neighbors[4], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[8], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[1], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[5], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[2], N, W, S, E, false);
			}
			else if (E && N && W)
			{
				invalidate_extrenal_cell(icell + _neighbors[1], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[5], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[2], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[6], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[3], N, W, S, E, false);
			}
			else if (E && N) 
			{
				invalidate_extrenal_cell(icell + _neighbors[1], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[5], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[2], N, W, S, E, false);
			}
			else if (N && W)
			{
				invalidate_extrenal_cell(icell + _neighbors[2], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[6], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[3], N, W, S, E, false);				
			}			
			else if (W && S)
			{
				invalidate_extrenal_cell(icell + _neighbors[3], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[7], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[4], N, W, S, E, false);
			}
			else if (S && E)
			{
				invalidate_extrenal_cell(icell + _neighbors[4], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[8], N, W, S, E, false);
				invalidate_extrenal_cell(icell + _neighbors[1], N, W, S, E, false);				
			}
			else if (N) invalidate_extrenal_cell(icell + _neighbors[2], N, W, S, E, false);
			else if (W) invalidate_extrenal_cell(icell + _neighbors[3], N, W, S, E, false);
			else if (S) invalidate_extrenal_cell(icell + _neighbors[4], N, W, S, E, false);
			else if (E) invalidate_extrenal_cell(icell + _neighbors[1], N, W, S, E, false);
			#pragma endregion
		};

		inline void invalidate_cell(size_t i, size_t j, bool N, bool W, bool S, bool E, size_t ied)
		{
			#pragma region Invalidate Cell Extended in a given direction:
			size_t icell; get_cell_handle(i, j, icell);
			
			_cell_edge_iter = _cell_first_edge_map.find(icell);
			if (_cell_edge_iter == _cell_first_edge_map.end())						
			{							 
				_cell_first_edge_map[icell] = ied;						 
			}					
			else if (_cell_edge_iter->second != ied && _cell_second_edge_map.find(icell) == _cell_second_edge_map.end())						 
			{							 
				_cell_second_edge_map[icell] = ied;
			}
			if (!_bad_cells[icell])
			{
				_bad_cells[icell] = true; _num_bad_cells++;
			}
			return;
			#pragma endregion
		};
	
        void invalidate_external_cells()
		{
			#pragma region invalidating External Cells:        
			size_t icell;
			for (size_t i = 0; i <= _imax; i++)
			{                
				for (size_t j = 0; j <= _jmax; j++)
				{
					get_cell_handle(i, j, icell);
					if (_bad_cells[icell]) break;
					_bad_cells[icell] = true; _num_bad_cells++;                    
				}
				for (size_t j = _jmax; j >= 0; j--)
				{
					get_cell_handle(i, j, icell);
					if (_bad_cells[icell]) break;
					_bad_cells[icell] = true; _num_bad_cells++;                    
					if (j == 0) break;
				}
			}

			for (size_t j = 0; j <= _jmax; j++)
			{                
				for (size_t i = 0; i <= _imax; i++)
				{
					get_cell_handle(i, j, icell);
					if (_bad_cells[icell]) break;
					_bad_cells[icell] = true; _num_bad_cells++;                    
				}
				for (size_t i = _imax; i >= 0; i--)
				{
					get_cell_handle(i, j, icell);
					if (_bad_cells[icell]) break;
					_bad_cells[icell] = true; _num_bad_cells++;                    
					if (i == 0) break;
				}
			}
			#pragma endregion
		};

		// the output is in _neighbor_points
		inline void get_neighbor_points(size_t icell)
		{
			#pragma region Retrieve Neighbors Points:						

			// reset _neighbor_points
			//for (size_t ii = 0; ii < _num_neighbor_points; ii++) _neighbor_points[ii] = 0;

			size_t jj(0); Point* pj; _num_neighbor_points_dynamic = 0;
			for (size_t ii = 0; ii < _num_neighbor_cells; ii++)
			{
				size_t jcell = icell + _neighbors[ii];

				//if (jcell >= _num_cells) continue;				

				//if (!_bad_cells[jcell]) continue; 

				if (_cell_points[jcell] == 0) continue; // _bad cells due to boundaries								

				pj = _cell_points[jcell];
				while (pj != 0)
				{
					_neighbor_points[jj] = pj; 
					_neighbor_points_cells[jj] = jcell;
					jj++; _num_neighbor_points_dynamic++;
					pj = pj->next;
				}
			}
			#pragma endregion
		};		

		inline bool is_valid_poly(PolyPoint* pp)
		{
			#pragma region Delete a poly:
			if (pp == 0) return true;
					
			double cx, cy;
			PolyPoint* pi = pp;
			PolyPoint* pj;
			bool valid;
			while (true)
			{
				valid = true;
				if (pi->center_L != 0)
				{
					valid = false;				
					cx = pi->center_L->x;					
					cy = pi->center_L->y;
					pj = pp;
					while (true)
					{
						// check if a pj lies outside of cirlce i					
						double dd = distance_squared(cx - pj->x, cy - pj->y);

						if (dd > _dm_squared + 1E-10) {valid = true; break;}
						if (pj->next == 0) break;
						pj = pj-> next;
						if (pj == pp) break;
					}
					if (!valid) 
					{
						//plot_cell("invalid_poly.ps", _icell, pp, 0); 
						return false;
					}
				}
				if (pi->center_R != 0)
				{
					valid = false;				
					cx = pi->center_R->x;					
					cy = pi->center_R->y;
					pj = pp;
					while (true)
					{
						// check if a pj lies outside of cirlce i					
						double dd = distance_squared(cx - pj->x, cy - pj->y);

						if (dd > _dm_squared + 1E-10) {valid = true; break;}
						if (pj->next == 0) break;
						pj = pj-> next;
						if (pj == pp) break;
					}
					if (!valid) 
					{
						//plot_cell("invalid_poly.ps", _icell, pp, 0); 
						return false;
					}
				}
				if (pi->next == 0) break;
				pi = pi-> next;
				if (pi == pp) break;
			}
			return true;
			#pragma endregion 
		};

		inline void delete_poly(PolyPoint* &pp)
		{
			#pragma region Delete a poly:
			if (pp == 0) return;
			
			PolyPoint* po = pp; PolyPoint* pi;
			while (true)
			{
				pi = pp->next; delete pp;
				if (pi == po || pi == 0) break;
				pp = pi;
			}
			pp = 0;
			#pragma endregion 
		};

		inline bool point_in_convex_hull(double xx, double yy, PolyPoint* pp)
		{			
			#pragma region Check if a given point lies in a convex hull:
			PolyPoint* po(pp);
			while (true)
			{
				if (pp->x * (pp->next->y - yy) + pp->next->x * (yy - pp->y) + xx * (pp->y - pp->next->y) < 0.0) 					
				{						
					return false;
				}
				pp = pp->next;
				if (pp == po) break;
			}
			return true;
			#pragma endregion
		};

		inline bool insert_point_in_poly(PolyPoint* pp, size_t icell)
		{
			#pragma region Insert a point in a poly:
			double xx, yy;
							
			pick_a_random_point_from_a_convex_hull(pp, xx, yy); // 1				
			if (is_valid_internal_point(icell, xx, yy))
			{
				Point* p = new Point();
				p->x = xx; p->y = yy;
				p->next = 0; p->edges = 0;
				p->done = false; p->on_boundary = false; p->num_edges = 0;
				Point* q = _cell_points[icell];
				if (q == 0) _cell_points[icell] = p;
				else
				{
					while (q->next != 0) q = q->next;
					q->next = p;
				}
				_num_sprinkled++;
				_bad_cells[icell] = true; _num_bad_cells++;											
				return true;				
			}			
			return false;
			#pragma endregion			
		};

		inline bool insert_point_in_polys(size_t icell)
		{
			#pragma region insert a point in a number of separated polys:
			double total_area(0.0);
			for (size_t ipoly = 0; ipoly < _num_polys; ipoly++)
			{
				if (_polys[ipoly] == 0) break; // the remaining polys are null

				_polys_area[ipoly] = 0.0;
				PolyPoint* po =  _polys[ipoly];
				PolyPoint* pp =  po;
				while (true)
				{
					_polys_area[ipoly] += pp->x * (pp->next->y - pp->prev->y);
					pp = pp->next;
					if (pp == po) break;
				}
				total_area += _polys_area[ipoly];
			}

			for (size_t ipoly = 0; ipoly < _num_polys; ipoly++)
			{
				if (_polys[ipoly] == 0) break; // the remaining polys are null

				_polys_area[ipoly] /= total_area;
			}

			double u = generate_a_random_number();

			total_area = 0.0;
			for (size_t ipoly = 0; ipoly < _num_polys; ipoly++)
			{
				if (_polys[ipoly] == 0) continue;

				total_area +=  _polys_area[ipoly];
				if (u < total_area)
				{
					return insert_point_in_poly(_polys[ipoly], icell);
				}
			}
			return false;
			#pragma endregion
		};

		inline bool isolated_cell(size_t icell)
		{
			#pragma region Check If a cell is isolated:
			for (size_t ii = 0; ii < _num_neighbor_cells; ii++)
			{
				if (!_bad_cells[icell + _neighbors[ii]])  return false;				
			}
			return true;
			#pragma endregion
		};

		inline void create_cell_polygon(size_t icell, PolyPoint* &pp, size_t &num_poly_pnts)
		{
			#pragma region create cell polygon:
			num_poly_pnts = 4;
			// constrcuct the approximating polygon
			PolyPoint* p1 = new PolyPoint();
			PolyPoint* p2 = new PolyPoint();
			PolyPoint* p3 = new PolyPoint();
			PolyPoint* p4 = new PolyPoint();

			size_t i,j;
			get_cell_indices(icell, i, j);
			double xj = _bg_xo + i * _bg_s;
			double yj = _bg_yo + j * _bg_s;

			p1->x = xj;         p1->y = yj;         p1->next = p2; p1->prev = p4; p1->center_L = 0; p1->center_R = 0;
			p2->x = xj + _bg_s; p2->y = yj;         p2->next = p3; p2->prev = p1; p2->center_L = 0; p2->center_R = 0;
			p3->x = xj + _bg_s; p3->y = yj + _bg_s; p3->next = p4; p3->prev = p2; p3->center_L = 0; p3->center_R = 0;
			p4->x = xj;         p4->y = yj + _bg_s; p4->next = p1; p4->prev = p3; p4->center_L = 0; p4->center_R = 0;
			
			pp = p1;
			#pragma endregion
		};			

		inline void generate_approximate_poly(size_t icell, PolyPoint* &pp, size_t &num_poly_pnts)
		{
			#pragma region Generate Initial Polyline
			_icell = icell;

			pp = 0; size_t i,j;
					
			num_poly_pnts = 4;

			// constrcuct the approximating polygon
			PolyPoint* p1 = new PolyPoint();
			PolyPoint* p2 = new PolyPoint();
			PolyPoint* p3 = new PolyPoint();
			PolyPoint* p4 = new PolyPoint();

			get_cell_indices(icell, i, j);
			double xj = _bg_xo + i * _bg_s;
			double yj = _bg_yo + j * _bg_s;

			p1->x = xj;         p1->y = yj;         p1->next = p2; p1->prev = p4; p1->center_L = 0; p1->center_R = 0;
			p2->x = xj + _bg_s; p2->y = yj;         p2->next = p3; p2->prev = p1; p2->center_L = 0; p2->center_R = 0;
			p3->x = xj + _bg_s; p3->y = yj + _bg_s; p3->next = p4; p3->prev = p2; p3->center_L = 0; p3->center_R = 0;
			p4->x = xj;         p4->y = yj + _bg_s; p4->next = p1; p4->prev = p3; p4->center_L = 0; p4->center_R = 0;

			pp = p1;				

			for (size_t ii = 0; ii < 120; ii++) _circles[ii] = 0;

			int jj = 0;
			for (size_t ii = 0; ii < _num_neighbor_cells; ii++)
			{				

				size_t jcell = icell + _neighbors[ii];

				if (_cell_points[jcell] == 0) continue; // _bad cells due to boundaries								

				_circles[jj] = _cell_points[jcell]; jj++;

				Point* q =  _cell_points[jcell]->next;
				while (q != 0)
				{
					_circles[jj] = q; jj++;
					q = q->next;
				}
			}

			if (_blue_cells[icell])
			{
				size_t ied;
				double ed_x1, ed_y1, ed_x2, ed_y2, ed_x3, ed_y3;
				_cell_edge_iter = _cell_second_edge_map.find(icell);
				if (_cell_edge_iter != _cell_second_edge_map.end())
				{
					delete_poly(pp);
					num_poly_pnts = 0;
					return;

					// cell is intersected using two edges
					ied = _cell_edge_iter->second;
					ed_x1 = _edges_x1[ied]; ed_y1 = _edges_y1[ied];
					ed_x2 = _edges_x2[ied]; ed_y2 = _edges_y2[ied];			

					_cell_edge_iter = _cell_first_edge_map.find(icell);
					ied = _cell_edge_iter->second;
					ed_x3 = _edges_x1[ied]; ed_y3 = _edges_y1[ied];				

					if (distance_squared(ed_x3 - ed_x2, ed_y3 - ed_y2) < 1E-10)
					{
						ed_x3 = _edges_x2[ied]; ed_y3 = _edges_y2[ied];				
					}
					else
					{
						double tmp = ed_x1; ed_x1 = ed_x3; ed_x3 = ed_x2; ed_x2 = tmp;
						tmp = ed_y1; ed_y1 = ed_y3; ed_y3 = ed_y2; ed_y2 = tmp;					
					}
				}
				else
				{
					_cell_edge_iter = _cell_first_edge_map.find(icell);
					ied = _cell_edge_iter->second;
					ed_x1 = _edges_x1[ied]; ed_y1 = _edges_y1[ied];
					ed_x2 = _edges_x2[ied]; ed_y2 = _edges_y2[ied];
				
					cut_polygon(pp, num_poly_pnts, ed_x1, ed_y1, ed_x2, ed_y2);

					if (num_poly_pnts == 0) return;					
				}

			}

			cut_polygon(pp, num_poly_pnts, _circles, jj);						
			#pragma endregion
		};

		inline void cut_polygon(PolyPoint* &pp, size_t &num_poly_pnts, double x1, double y1, double x2, double y2)
		{
			#pragma region Cut Polygon With an Edges:
			PolyPoint* po; PolyPoint* p1; PolyPoint* p2; PolyPoint* pbad; PolyPoint* ptmp;
			size_t num_in(0);
			for (size_t ipnt = 0; ipnt < num_poly_pnts; ipnt++)
			{
				if (invalid_triangle(pp->x, pp->y, x1, y1, x2, y2, 1E-10)) pp->u = -1.0;
				else 
				{
					pp->u = 1.0; num_in++;
				}
				pp = pp->next;
			}

			if (num_in == 0)
			{
				delete_poly(pp); num_poly_pnts = 0; return;
			}

			for (size_t ipnt = 0; ipnt < num_poly_pnts; ipnt++)
			{				
				if (pp->u > 0.5 && pp->next->u < 0.5)
				{
					// pp is in pp->next is out
					pbad = pp->next;
					po = pp->next->next;
					while (po->u < 0.5) po = po->next;

					double xst, yst;					
					p1 = new PolyPoint(); num_poly_pnts++;
					get_segments_intersection(x1, y1, x2, y2, pp->x, pp->y, pp->next->x, pp->next->y, xst, yst);
					p1->x = xst; p1->y = yst;
					get_segments_intersection(x1, y1, x2, y2, po->x, po->y, po->prev->x, po->prev->y, xst, yst);
					p2 = new PolyPoint(); num_poly_pnts++;
					p2->x = xst; p2->y = yst;

					p1->next = p2; p2->prev = p1;
					p1->prev = pp; pp->next = p1;
					p2->next = po; po->prev = p2;
					p1->center_L = 0; p1->center_R = 0; 
					p2->center_L = 0; p2->center_R = 0; 

					// delete all the points between pbad and po
					while (pbad != po)
					{
						ptmp = pbad;
						pbad = pbad->next;
						delete ptmp; num_poly_pnts--;
					}
					break;
				}
				pp = pp->next;
			}		
			#pragma endregion
		};
		 
		inline double distance_squared(double dx, double dy)
        {
            return dx * dx + dy * dy;        
        };

		inline bool invalid_triangle(double x1, double y1, double x2, double y2, double x3, double y3, double TOL)
		{
			if ((x1 - x3) * (y2 - y3) - (y1 - y3) * (x2 - x3) < -TOL) return true;
			else return false;
		};


		inline bool ZeroAngle(double x1, double y1, double x2, double y2, double x3, double y3, double TOL)
		{
			if ((x3 - x1) * (x2 - x1) < TOL) return false;
			if ((y3 - y1) * (y2 - y1) < TOL) return false;
			if (fabs((x3 - x1) * (y2 - y1) - (x2 - x1) * (y3 - y1)) < TOL) return true;
			else return false;
		};

		inline double area_triangle(double x1, double y1, double x2, double y2, double x3, double y3)
        {
            return 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
        };

		inline double area_parallelogram(double x1, double y1, double x2, double y2, double x3, double y3)
        {
            return (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
        };

		inline bool point_in_triangle(double xx, double yy, double x1, double y1, double x2, double y2, double x3, double y3)
		{
			#pragma region Point In Triangle test:
			double T((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3));
			double uT = (xx - x3) * (y2 - y3) + (x3 - x2) * (yy - y3);
			if ((uT < 0.0 && T > 0.0) || (uT > 0.0 && T < 0.0)) return false;
			if ((uT > T && T > 0.0) || (uT < T && T < 0.0)) return false;

			double vT = (xx - x3) * (y3 - y1) + (x1 - x3) * (yy - y3);
			if ((vT < 0.0 && T > 0.0) || (vT > 0.0 && T < 0.0)) return false;
			if ((vT > T && T > 0.0) || (vT < T && T < 0.0)) return false;

			if ((uT + vT > T && T > 0.0) || (uT + vT < T && T < 0.0)) return false;

			return true;
			#pragma endregion
		};

		inline bool crossing_segments(double x1, double y1, double x2, double y2,
                                      double x3, double y3, double x4, double y4)
        {
			#pragma region Crossing segments tests:
            double a1 = area_triangle(x1, y1, x3, y3, x4, y4);
			double a2 = area_triangle(x2, y2, x3, y3, x4, y4);
			if (a1 > -1E-10 && a2 > -1E-10) return false;
			if (a1 < 1E-10 && a2 < 1E-10) return false;			
			a1 = area_triangle(x3, y3, x1, y1, x2, y2);
			a2 = area_triangle(x4, y4, x1, y1, x2, y2);
			if (a1 > -1E-10 && a2 > -1E-10) return false;
			if (a1 < 1E-10 && a2 < 1E-10) return false;
			return true;
			#pragma endregion
        };

		inline bool crossing_segments(double x1, double y1, double x2, double y2,
                                      double x3, double y3, double x4, double y4, double &xx, double &yy)
        {
			#pragma region Crossing segments tests:
            double a1 = area_triangle(x1, y1, x3, y3, x4, y4);
			double a2 = area_triangle(x2, y2, x3, y3, x4, y4);
			if (a1 > -1E-10 && a2 > -1E-10) return false;
			if (a1 <  1E-10 && a2 <  1E-10) return false;
			a1 = area_triangle(x3, y3, x1, y1, x2, y2);
			a2 = area_triangle(x4, y4, x1, y1, x2, y2);
			if (a1 > -1E-10 && a2 > -1E-10) return false;
			if (a1 <  1E-10 && a2 <  1E-10) return false;
			double L12 = sqrt(distance_squared(x2 - x1, y2 - y1));
			double h3 = fabs(a1 / L12);
			double h4 = fabs(a2 / L12);
			double r = h3 / (h3 + h4);
			xx = x3 + r * (x4 - x3);
			yy = y3 + r * (y4 - y3);
			return true;
			#pragma endregion
        };

		inline bool intersecting_segments(double x1, double y1, double x2, double y2,
                                      double x3, double y3, double x4, double y4, double &xx, double &yy)
        {
			#pragma region Intersecting segments tests:
			// Intersection does not include end points --> used only to adjust boundary edges
			double Lsq = distance_squared(x1 - x3, y1 - y3);
			if (Lsq < 1E-10) return false;
			Lsq = distance_squared(x1 - x4, y1 - y4);
			if (Lsq < 1E-10) return false;
			Lsq = distance_squared(x2 - x3, y2 - y3);
			if (Lsq < 1E-10) return false;
			Lsq = distance_squared(x2 - x4, y2 - y4);
			if (Lsq < 1E-10) return false;

            double a1 = area_triangle(x1, y1, x3, y3, x4, y4);
			double a2 = area_triangle(x2, y2, x3, y3, x4, y4);
			if (a1 > 1E-10 && a2 > 1E-10) return false;
			if (a1 < -1E-10 && a2 < -1E-10) return false;
			if (fabs(a1) < 1E-10 && fabs(a2) < 1E-10) return false;
			a1 = area_triangle(x3, y3, x1, y1, x2, y2);
			a2 = area_triangle(x4, y4, x1, y1, x2, y2);
			if (a1 > 1E-10 && a2 > 1E-10) return false;
			if (a1 < -1E-10 && a2 < -1E-10) return false;
			if (fabs(a1) < 1E-10 && fabs(a2) < 1E-10) return false;
			double L12 = sqrt(distance_squared(x2 - x1, y2 - y1));
			double h3 = fabs(a1 / L12);
			double h4 = fabs(a2 / L12);
			double r = h3 / (h3 + h4);
			xx = x3 + r * (x4 - x3);
			yy = y3 + r * (y4 - y3);
			return true;
			#pragma endregion
        };

		inline bool get_segments_intersection(double x1, double y1, double x2, double y2,
                                  double x3, double y3, double x4, double y4,
                                  double &px, double &py)
        {
			#pragma region retrieve the intersection of two segments
			if (fabs(x3-x1) < 1E-10 && fabs(y3-y1) < 1E-10) {px = x1; py = y1; return true;}
			else if (fabs(x4-x1) < 1E-10 && fabs(y4-y1) < 1E-10) {px = x1; py = y1; return true;}
			else if (fabs(x3-x2) < 1E-10 && fabs(y3-y2) < 1E-10) {px = x2; py = y2; return true;}
			else if (fabs(x4-x2) < 1E-10 && fabs(y4-y2) < 1E-10) {px = x2; py = y2; return true;}

            double rx(x4 - x3), ry(y4 - y3);		
            double r = sqrt(rx * rx + ry * ry);			
            rx /= r; ry /= r;
            double nx(y2 - y1), ny(x1 - x2);			
            double dot = nx * rx + ny * ry;
            if (fabs(dot) < 1E-10) return false; // ray is parallel to line segment
			
            double dst = (nx * (x1 - x3) + ny * (y1 - y3)) / dot;            
            px = x3 + dst * rx; py = y3 + dst * ry;           
            if (px < std::min(x1, x2) - 1E-10) return false; // segments do not intersect
            if (px > std::max(x1, x2) + 1E-10) return false; // segments do not intersect
            if (py < std::min(y1, y2) - 1E-10) return false; // segments do not intersect
            if (py > std::max(y1, y2) + 1E-10) return false; // segments do not intersect           
            return true;
			#pragma endregion
        };

		inline bool incircumcircle(double ax, double ay, double bx, double by, double cx, double cy, double dx, double dy)
		{
			#pragma region Incircle test:
			double adx, ady, bdx, bdy, cdx, cdy;
			double abdet, bcdet, cadet;
			double alift, blift, clift;
	
			// 6 flops
			adx = ax - dx;			 
			ady = ay - dy;			  
			bdx = bx - dx;
			bdy = by - dy;
			cdx = cx - dx;
			cdy = cy - dy;

			// 18 flops
			abdet = adx * bdy - bdx * ady;
			bcdet = bdx * cdy - cdx * bdy;
			cadet = cdx * ady - adx * cdy;
			alift = adx * adx + ady * ady;	
			blift = bdx * bdx + bdy * bdy;
			clift = cdx * cdx + cdy * cdy;

			// 5 flops
			if (alift * bcdet + blift * cadet + clift * abdet > 1E-10) return false;
			
			return true;
			#pragma endregion
		};
		

		inline bool circumcircle(double x1, double y1, double x2, double y2, double x3, double y3, 
                          double &cx, double &cy, double &cr_sq)
		{			
			#pragma region find the circumcircle of three points:			
			// translate c.s. to (x1,y1) --> 4 flops
			x2 -= x1; y2 -= y1;
			x3 -= x1; y3 -= y1;
			
			// calculate (x0p, y0p)    --> 3 flops
			double det(x3 * y2 - x2 * y3);
    
			if (fabs(det) < 1E-10) {cr_sq = -1.0; return false;} // a degenerate circle			

			det = 1.0 / det; // 1 flop   

			double b1, b2; // 8 flops
			b1 = 0.5 * (x3 * x3 + y3 * y3);
			b2 = 0.5 * (x2 * x2 + y2 * y2);
    
			// 11 flops
			cx = (b1 * y2 - b2 * y3) * det;
			cy = (x3 * b2 - x2 * b1) * det;    
			cr_sq = cx * cx + cy * cy;
    
			// 2 flops
			// shift axes back    
			cx += x1; cy += y1; 
			return true;
			#pragma endregion
		};				

		bool get_edge(size_t icell, size_t jcell, size_t &ied)
		{
			#pragma region Retrieve an Edge passing through two cells:
			if (!_blue_cells[icell] || !_blue_cells[jcell]) return false;
			size_t ed_11 = _cell_first_edge_map[icell];
			size_t ed_21 = _cell_first_edge_map[jcell];
			if (ed_11 == ed_21)
			{
				ied = ed_11; return true;
			}
			size_t ed_12(ed_11);
			if (_cell_second_edge_map.find(icell) != _cell_second_edge_map.end()) ed_12 =  _cell_second_edge_map[icell];
			if (ed_12 == ed_21)
			{
				ied = ed_12; return true;
			}
			size_t ed_22(ed_21);
			if (_cell_second_edge_map.find(jcell) != _cell_second_edge_map.end()) ed_22 =  _cell_second_edge_map[jcell];
			if (ed_22 == ed_12 || ed_22 == ed_11)
			{
				ied = ed_22; return true;
			}
			return false;								
			#pragma endregion
		};

		size_t generate_edge_id(size_t edge_index, bool forward_direction)
		{
			size_t edge_id  = edge_index << 1;
			edge_id += forward_direction;
			edge_id++;
			return edge_id;
		};

		void get_edge_data(size_t edge_id, size_t &edge_index, bool &forward_direction)
		{			
			forward_direction = bool((edge_id - 1) & 1);
			edge_index = (edge_id - 1) >> 1;			
		};

		int save_point_cloud(std::string file_name)
		{
			#pragma region Save Point Cloud:
			std::fstream file(file_name.c_str(), std::ios::out);
			Point* q;
			for (size_t icell = 0; icell < _num_cells; icell++)
			{
				q =  _cell_points[icell];
				if (q == 0) continue;		
				file  << icell << " " << q->x << " " << q->y << std::endl;				
				
				while (q->next!=0)
				{
					q = q->next;
					file  << icell << " " << q->x << " " << q->y << std::endl;				
				}
			}
			return 0;		
			#pragma endregion
		};

		int save_boundary_points(std::string file_name, std::list<Point*> &boundary_points, std::list<size_t> &boundary_edges)
		{
			#pragma region Save Boundary Points:
			std::fstream file(file_name.c_str(), std::ios::out);
			std::list<Point*>::iterator iter;
			std::list<size_t>::iterator eiter = boundary_edges.begin();
			for (iter = boundary_points.begin(); iter != boundary_points.end(); iter++)
			{
				file << (*iter)->x << " " << (*iter)->y << " " << *eiter << std::endl;		
				eiter++;
			}
			return 0;		
			#pragma endregion
		};

		int load_boundary_points(std::string file_name, std::vector<double> &x, std::vector<double> &y, std::vector<size_t> &edges)
		{
			#pragma region Load Boundary Points:
			x.clear(); y.clear();
			std::string line; size_t num_points(0);
			std::ifstream myfile(file_name.c_str());
			if (myfile.is_open()) 
			{    
				while (!myfile.eof()) 
				{
					getline(myfile, line);
					if (line.size() == 0) break;
            
					std::istringstream iss(line);
					std::vector<std::string> tokens;
					copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
						std::back_inserter<std::vector<std::string> >(tokens));
					double xx = string_to_double(tokens[0]);
					double yy = string_to_double(tokens[1]);
					size_t iedge = (size_t) (string_to_double(tokens[2]));
					x.push_back(xx);
					y.push_back(yy);
					edges.push_back(iedge);
				}
			}
			return 0;
			#pragma endregion
		};
		
		// utilized in debugging the CDT/Voronoi algroithm
		int debug_point_cloud(std::string file_name)
		{
			#pragma region Check the tessellation code:		
			_num_neighbor_points = 100;
			_neighbor_points = new Point*[_num_neighbor_points];
			_neighbor_points_cells = new size_t[_num_neighbor_points];
			for (size_t i = 0; i < _num_neighbor_points; i++) _neighbor_points[i] = 0;					

			std::string line; size_t num_points(0);
			std::ifstream myfile(file_name.c_str());
			if (myfile.is_open()) 
			{    
				while (!myfile.eof()) 
				{
					getline(myfile, line);
					if (line.size() == 0) break;
            
					std::istringstream iss(line);
					std::vector<std::string> tokens;
					copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
						std::back_inserter<std::vector<std::string> >(tokens));
					size_t icell = size_t(string_to_double(tokens[0]));
					double x = string_to_double(tokens[1]);
					double y = string_to_double(tokens[2]);
								
					if (_cell_points[icell] == 0)
					{
						_cell_points[icell] = new Point();
						_cell_points[icell]->x = x;
						_cell_points[icell]->y = y;	
						_cell_points[icell]->edges = 0;
						_cell_points[icell]->next = 0;
						_cell_points[icell]->done = false;
						_cell_points[icell]->on_boundary = false;
						_bad_cells[icell] = true;
						num_points++;
					}
					else
					{
						Point* q = _cell_points[icell];
						while (q->next != 0) q = q->next;

						q->next = new Point();
						q->next->x = x;
						q->y = y;						
						q->edges = 0;
						q->next = 0;
						q->done = false;
						q->on_boundary = false;
						_bad_cells[icell] = true;
						num_points++;
					}
				}
			}
			std::vector<double> vx; std::vector<double> vy; std::vector<size_t> ve;
			load_boundary_points("boundary_points.dat", vx, vy, ve);
			size_t icell; double dst_sq;
			int num(vx.size());
			for (int i = 0; i < num; i++)
			{
				Point* p = closest_point(vx[i], vy[i], 1E-5, dst_sq, icell); i++;
				Point* q = closest_point(vx[i], vy[i], 1E-5, dst_sq, icell);
				add_constraint_edge(p, q, ve[i]);
			}
			plot_vertices("sprinkled_vertices_grid_circles_.ps", true, true);
			_circles = new Point*[120];
			for (int ii = 0; ii < 120; ii++) _circles[ii] = 0;
			eliminate_voids_via_polys();
			plot_vertices("sprinkled_vertices_grid_circles_after.ps", true, true);
			CDT(true);
			return 0;
			#pragma endregion
		};

		void create_element(std::list<double> &vxy, std::list<double>::iterator &iter,  
			                std::vector<bool> &v_on_boundary,
			                std::vector<Point*> &element, double TOL)
		{
			#pragma region create a voronoi cell:
			element.clear(); Point* q; Point* first_point(0); Point* last_point(0);
			iter = vxy.begin(); int i(0);
			while (iter != vxy.end())
			{
				double xx = *iter; iter++;
				double yy = *iter; iter++;
				get_closest_point(xx, yy, v_on_boundary[i], q, TOL);
				if (q != first_point && q != last_point) element.push_back(q); i++;
				if (first_point == 0) first_point = q;
				else last_point = q;
			}
			vxy.clear();
			#pragma endregion
		};

		Point* closest_point(double xx, double yy, double range_sq, double &dst_sq, size_t &icell)
		{
			#pragma region Retrieve Closest Point:
			Point* q = 0; size_t icell_;
			get_cell_handle(xx, yy, icell_);
			icell = icell_;
			double dmin = range_sq;
			bool first(true); Point* p;
			for (size_t i = 0; i < 9; i++)
			{
				size_t jcell = icell_ + _neighbors[i];
				if (_cell_points[jcell] == 0) continue;
				if (first) p = _cell_points[jcell];
				double dd = distance_squared(p->x - xx, p->y - yy);
				if (dd < dmin)
				{
					q = p;
					dmin = dd;
					icell = jcell;
				}
				if (p->next != 0) {first = false; p = p->next; i--;}
				else first = true;
			}
			dst_sq = dmin;
			return q;
			#pragma endregion
		};

		bool get_closest_point(double xx, double yy, bool on_boundary, Point* &q, double TOL)
		{
			#pragma region Retrieve Closest Point:
			q = 0; size_t icell;
			get_cell_handle(xx, yy, icell);
			bool first(true); Point* p;
			for (size_t i = 0; i < 9; i++)
			{
				size_t jcell = icell + _neighbors[i];
				if (_cell_nodes[jcell] == 0) continue;
				if (first) p = _cell_nodes[jcell];
				double dd = distance_squared(p->x - xx, p->y - yy);
				if (dd < TOL)
				{
					if (on_boundary && p->on_boundary) continue;
					q = p; 
					if (on_boundary)
					{
						q->x = xx; q->y = yy;
					}
					return true;
				}
				if (p->next != 0) {first = false; p = p->next; i--;}
				else first = true;
			}

			q = new Point();
			q->x = xx; q->y = yy; q->next = 0; q->edges = 0;
			q->done = false; q->on_boundary = false;
			
			if (_cell_nodes[icell] == 0) {_cell_nodes[icell] = q; return false;}
			
			p = _cell_nodes[icell];
			while (p->next!= 0) 
			{
				p = p->next;
			}
			p->next = q;
			return false; // a new node was created
			#pragma endregion
		};

    private:       

		bool _failed; 
		bool _debug_code;
		bool _save_point_cloud;

		// input:
        double _dm, _dm_squared; // uniform distance function
		
		std::vector<double>				 _ExternalBoundaries;
		std::vector< std::vector<double> > _Holes;
		std::vector< std::vector<double> > _Cracks;
		int _output; // 0 : point cloud - 1: Voronoi Cells - 2: CDT - 3: All-quad
		bool _plot_results;
		double _xmin, _ymin, _xmax, _ymax;
		
		size_t _icell;

        // Data for boundary edges
        size_t                 _num_boundary_edges; // external boundaries
        size_t                 _num_holes_edges;
		size_t                 _num_oriented_edges;
        std::vector<double> _edges_x1; std::vector<double> _edges_y1;
        std::vector<double> _edges_x2; std::vector<double> _edges_y2;

		std::map<size_t, size_t> _cell_first_edge_map;
		std::map<size_t, size_t> _cell_second_edge_map;
		std::map<size_t, size_t>::iterator _cell_edge_iter;
		
        // background grid
        double              _bg_xo, _bg_yo, _bg_s;        
        size_t              _nr, _nc, _imax, _jmax;
		size_t              _num_cells, _num_bad_cells, _num_sprinkled;
        bool*               _bad_cells;
		bool*               _blue_cells;

		int*                _neighbors;
		Point**             _cell_points; // inserted points
		
		Point**				_cell_nodes; // for Voronoi mesh		

		size_t				_num_neighbor_cells;
		size_t				_num_neighbor_points;
		size_t				_num_neighbor_points_dynamic;
		Point**             _neighbor_points;
		size_t*             _neighbor_points_cells;

		Point**             _circles; // Points in the neighbors cells

		// polygonal voids
		PolyPoint**         _polys; // assume an upper bound of _num_polys
		size_t              _num_polys;
		double*             _polys_area;    

		std::vector< std::vector<size_t> > _side_sets;

		size_t           _max_num_poly_pnts;
        size_t           _MY_RAND_MAX;
		size_t           _fixed_seed;
};
                                
#endif	



