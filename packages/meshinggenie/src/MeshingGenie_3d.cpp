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


#include "MeshingGenie_3d.h"


using namespace std;

MeshingGenie_3d::MeshingGenie_3d(std::vector<double> &xb, std::vector<double> &yb, std::vector<double> &zb, 
				 std::vector<std::vector<size_t> > &bfaces, double r, double tol, 
				 std::vector< std::vector<double> > &xf, std::vector< std::vector<double> > &yf, std::vector< std::vector<double> > &zf)
{
	// Constants

	scale_input(xb, yb, zb, r, tol, xf, yf, zf);

	_xf = xf;
	_yf = yf;
	_zf = zf;
	_nfx.resize(_xf.size());
	_nfy.resize(_xf.size());
	_nfz.resize(_xf.size());

	_MY_RAND_MAX = RAND_MAX; // 15 bits
    _MY_RAND_MAX = _MY_RAND_MAX << 15;
    _MY_RAND_MAX += RAND_MAX; // it's now 30 bits

	_SQRT_3 = 1.73205080756;
	_SQRT_3_INV = 0.577350269189;

	_RF = 0.7;

	_tol = tol;

	_dm = r; _xb = xb; _yb = yb; _zb = zb;
	_bfaces = bfaces;
	
	// computing normals
	size_t num_faces(bfaces.size());
	_bfaces_nx.clear(); _bfaces_nx.resize(num_faces);
	_bfaces_ny.clear(); _bfaces_ny.resize(num_faces);
	_bfaces_nz.clear(); _bfaces_nz.resize(num_faces);
	_nonconvex_bfaces.clear(); _nonconvex_bfaces.resize(num_faces);

	_b_node_faces.resize(xb.size());

	for (size_t iface = 0; iface < num_faces; iface++) 
	{
		_nonconvex_bfaces[iface] = false;
		size_t num_face_nodes(_bfaces[iface].size());
		for (size_t inode = 0; inode < num_face_nodes; inode++)
		{
			size_t i = _bfaces[iface][inode];
			_b_node_faces[i].push_back(iface);
		}
	}

	for (size_t iface = 0; iface < num_faces; iface++)
	{
		#pragma region Calculate Face normals:
		double x1 = _xb[bfaces[iface][0]]; double y1 = _yb[bfaces[iface][0]]; double z1 = _zb[bfaces[iface][0]];
		double x2 = _xb[bfaces[iface][1]]; double y2 = _yb[bfaces[iface][1]]; double z2 = _zb[bfaces[iface][1]];
		
		size_t iipp(2); size_t num_iface_corners(_bfaces[iface].size());
		double x3, y3, z3; bool invalid_face(true);
		for (size_t inode = 2; inode < num_iface_corners; inode++)
		{							
			#pragma region pick a point that is not colinear with jj and jjp:
			x3 = _xb[_bfaces[iface][inode]];
			y3 = _yb[_bfaces[iface][inode]];
			z3 = _zb[_bfaces[iface][inode]];
							
			double ax(x3 - x1), ay(y3 - y1), az(z3 - z1);
			double bx(x3 - x2), by(y3 - y2), bz(z3 - z2);

			double cx(ay * bz - az * by);
			double cy(az * bx - ax * bz);
			double cz(ax * by - ay * bx);

			if (fabs(cx) > 1E-10 || fabs(cy)> 1E-10 || fabs(cz) > 1E-10)
			{
				invalid_face = false;
				double ax(x2-x1), ay(y2-y1), az(z2-z1);
				double bx(x3-x1), by(y3-y1), bz(z3-z1);
				double cx = ay * bz - az * by;
				double cy = az * bx - ax * bz;
				double cz = ax * by - ay * bx;
				double c_inv = 1.0 / sqrt(cx * cx + cy * cy + cz * cz);
				_bfaces_nx[iface] = cx * c_inv;
				_bfaces_ny[iface] = cy * c_inv;
				_bfaces_nz[iface] = cz * c_inv;
				break;
			}
			#pragma endregion
		}
		#pragma endregion
	}

	size_t num_i_faces(_xf.size());
	for (size_t iface = 0; iface < num_i_faces; iface++)
	{
		#pragma region calculate internal face normal:
		double nx, ny, nz;
		size_t num_face_corners(_xf[iface].size());

		double x1 = _xf[iface][0]; double y1 = _yf[iface][0]; double z1 = _zf[iface][0];
		double x2 = _xf[iface][1]; double y2 = _yf[iface][1]; double z2 = _zf[iface][1];
		double x3, y3, z3; bool invalid_face(true);
	
		for (size_t inode = 2; inode < num_face_corners; inode++)
		{							
			#pragma region pick a point that is not colinear with jj and jjp:
			x3 = _xf[iface][inode];
			y3 = _yf[iface][inode];
			z3 = _zf[iface][inode];
							
			double ax(x3 - x1), ay(y3 - y1), az(z3 - z1);
			double bx(x3 - x2), by(y3 - y2), bz(z3 - z2);

			double cx(ay * bz - az * by);
			double cy(az * bx - ax * bz);
			double cz(ax * by - ay * bx);

			if (fabs(cx) > 1E-10 || fabs(cy)> 1E-10 || fabs(cz) > 1E-10)
			{
				invalid_face = false;
				double ax(x2-x1), ay(y2-y1), az(z2-z1);
				double bx(x3-x1), by(y3-y1), bz(z3-z1);
				double cx = ay * bz - az * by;
				double cy = az * bx - ax * bz;
				double cz = ax * by - ay * bx;
				double c_inv = 1.0 / sqrt(cx * cx + cy * cy + cz * cz);
				nx = cx * c_inv;
				ny = cy * c_inv;
				nz = cz * c_inv;
				break;
			}
			#pragma endregion
		}

		if (invalid_face) return;

		_nfx[iface] = nx;
		_nfy[iface] = ny;
		_nfz[iface] = nz;
		#pragma endregion
	}

	std::vector<size_t> edge_faces;
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		#pragma region retrieve non_convex_faces and sharp edges:
		bool non_convex(false);		
		double x1 = _xb[bfaces[iface][0]]; double y1 = _yb[bfaces[iface][0]]; double z1 = _zb[bfaces[iface][0]];
		size_t num_iface_corners(_bfaces[iface].size());
		for (size_t ii = 0; ii < num_iface_corners; ii++)
		{
			#pragma region retrieve non convex faces:
			size_t iip(ii + 1);
			if (iip == num_iface_corners) iip = 0;

			size_t node_i(_bfaces[iface][ii]), node_j(_bfaces[iface][iip]);
			if (node_i > node_j) continue; 
			
			get_b_edge_faces(node_i, node_j, edge_faces);
			
			if (edge_faces.size() == 1 || edge_faces.size() > 2) 
			{
				size_t num_edge_faces(edge_faces.size());
				for (size_t jj = 0; jj < num_edge_faces; jj++)
				{
					_nonconvex_bfaces[edge_faces[jj]] = true;
				
				}
				_sharp_edges.push_back(_bfaces[iface][ii]);
				_sharp_edges.push_back(_bfaces[iface][iip]);
				continue;
			}

			// edge has two faces
			size_t jface(edge_faces[0]);
			if (jface == iface) jface = edge_faces[1];

			// find common edge:
			size_t num_jface_corners(_bfaces[jface].size());			
			for (size_t jj = 0; jj < num_jface_corners; jj++)
			{
				size_t jjp(jj + 1);
				if (jjp == num_jface_corners) jjp = 0;
										
				if (_bfaces[iface][ii] == _bfaces[jface][jjp] && _bfaces[iface][iip] == _bfaces[jface][jj])
				{
					size_t jjpp(jjp + 1);
					double x1_(_xb[_bfaces[jface][jj]]), y1_(_yb[_bfaces[jface][jj]]), z1_(_zb[_bfaces[jface][jj]]);
					double x2_(_xb[_bfaces[jface][jjp]]), y2_(_yb[_bfaces[jface][jjp]]), z2_(_zb[_bfaces[jface][jjp]]);
					double xx_, yy_, zz_;
					while (true)
					{							
						#pragma region pick a point that is not colinear with jj and jjp:
						if (jjpp == num_jface_corners) jjpp = 0; 
						xx_ = _xb[_bfaces[jface][jjpp]];
						yy_ = _yb[_bfaces[jface][jjpp]];
						zz_ = _zb[_bfaces[jface][jjpp]];
							
						double ax(xx_ - x1_), ay(yy_ - y1_), az(zz_ - z1_);
						double bx(xx_ - x2_), by(yy_ - y2_), bz(zz_ - z2_);

						double cx(ay * bz - az * by);
						double cy(az * bx - ax * bz);
						double cz(ax * by - ay * bx);

						if (fabs(cx) > 1E-10 || fabs(cy)> 1E-10 || fabs(cz) > 1E-10)
						{
							break;
						}

						jjpp = jjpp + 1;
						#pragma endregion
					}

					// A common edge form a tetrahedron						
					double dot = (xx_ - x1) * _bfaces_nx[iface] + (yy_ - y1) * _bfaces_ny[iface] + (zz_ - z1) * _bfaces_nz[iface];
					if (dot > -1E-10)
					{
						// incluse palanar neighbors
						_nonconvex_bfaces[iface] = true;
						_nonconvex_bfaces[jface] = true;

						dot = _bfaces_nx[iface] * _bfaces_nx[jface] + 
								_bfaces_ny[iface] * _bfaces_ny[jface] + 
								_bfaces_nz[iface] * _bfaces_nz[jface];

						if (fabs(dot) < 0.5 && _bfaces[iface][ii] < _bfaces[iface][iip])
						{
							_sharp_edges.push_back(_bfaces[iface][ii]);
							_sharp_edges.push_back(_bfaces[iface][iip]);
						}

						non_convex = true;
					}
					break; // skip jface next corner
				}
			}
			#pragma endregion		
		}
		#pragma endregion
	}

}

int MeshingGenie_3d::execute(size_t fixed_seed)
{
	#pragma region Execute Function:

	if (fixed_seed == 0) fixed_seed = size_t(time(NULL));
	
	_fixed_seed = fixed_seed;

	srand(fixed_seed);

	std::cout<< "\n*** Random Seed = " << fixed_seed << "\n" << std::endl;


	clock_t start_time, end_time; double cpu_time;	 
	start_time = clock();

	/*
	save_vtk_polygons(_xb, _yb, _zb, _bfaces);
	std::vector< std::vector<size_t> > faces(1);
	faces[0].push_back(0);
	faces[0].push_back(1);
	faces[0].push_back(2);
	save_vtk_polygons(_xf[0], _yf[0], _zf[0], faces);
	save_vtk_polygons(_xf[1], _yf[1], _zf[1], faces);
	*/

	generate_background_grid();

	identify_external_cell();

	_num_inserted_points = 0;
	_cell_points = new Point*[_n];
	for (size_t icell = 0; icell < _n; icell++)
	{
		_cell_points[icell] = 0;
	}

	for (size_t icell = 0; icell < _n; icell++)
	{
		if (_invalid_cells[icell]) continue;

		if (_boundary_cells[icell] || _external_cells[icell]) 
		{
			_invalid_cells[icell] = true; _num_valid--;
		}
	}
				
	size_t jcell_o;
	get_cell_index(2, 2, 2, jcell_o);
	size_t ii, jj, kk, jcell;
	_neighbors = new int[124];
	for (int index = 0; index < 124; index++)
	{
		get_neighbor_indices(index, 2, 2, 2, ii, jj, kk);	
		get_cell_index(ii, jj, kk, jcell);
		_neighbors[index] = int(jcell) - int(jcell_o);
	}

	//plot_boundary_cells();
	//std::cout << " boundary cells plot done" << std::endl;
	//return 0;

	sprinkle_points_along_sharp_edges(0.5 * _dm);
	
	detect_internal_cracks_on_internal_faces();

	size_t num_internal_faces(_xf.size());
	for (size_t iface = 0; iface < num_internal_faces; iface++)
	{
		sample_face(iface);
	}

	size_t ref_levels(30);
	for (size_t iref = 0; iref <= ref_levels; iref++) 
	{
		use_grid_level(iref);
		if (_num_valid == 0) break;
	}

	end_time = clock();
	cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

	std::cout<< "==================================================================" << std::endl;
	std::cout<< "Number of inserted points = " << _num_inserted_points << std::endl;
	std::cout<< "Execution Time = " << cpu_time << " seconds." << std::endl;
	std::cout<< "Remaining Cells Ratio = " << _num_valid * 1.0 / _no << std::endl;
	std::cout<< "Point Insertion Rate = " << _num_inserted_points * 1.0 / cpu_time << " points / sec" << std::endl;
	std::cout<< "Point Insertion Rate = " << _num_inserted_points * 1.0 / _no << " points / cell" << std::endl;
	std::cout<< "==================================================================" << std::endl;

   
	generate_CVM();	

	std::cout<< "\n*** Mission Accomplished ***." << std::endl;

	return 0;
	#pragma endregion
}

void MeshingGenie_3d::generate_background_grid()
{
	#pragma region generate background grid:

	double xmin(_xb[0]), xmax(_xb[0]);
	double ymin(_yb[0]), ymax(_yb[0]);
	double zmin(_zb[0]), zmax(_zb[0]);

	size_t num_b_points(_xb.size());
	for (size_t i = 1; i < num_b_points; i++)
	{
		if (_xb[i] < xmin) xmin = _xb[i];
		if (_xb[i] > xmax) xmax = _xb[i];
		if (_yb[i] < ymin) ymin = _yb[i];
		if (_yb[i] > ymax) ymax = _yb[i];
		if (_zb[i] < zmin) zmin = _zb[i];
		if (_zb[i] > zmax) zmax = _zb[i];
	}
	_box_xmin = xmin;
	_box_xmax = xmax;
	_box_ymin = ymin;
	_box_ymax = ymax;
	_box_zmin = zmin;
	_box_zmax = zmax;

	double Lx = xmax - xmin;
	double Ly = ymax - ymin;
	double Lz = zmax - zmin;

	_s = _SQRT_3_INV * _dm;

	// 3 cells are added to each side so that each active cell have the same number of neighors 
	_ni = size_t(ceil(Lx / _s)) + 10; 
	_nj = size_t(ceil(Ly / _s)) + 10;
	_nk = size_t(ceil(Lz / _s)) + 10;

	_njk = _nj * _nk;
	_n = _ni * _njk;

	_dm_squared = _dm * _dm;
	
	_xo = xmin - 5 * _s;
	_yo = ymin - 5 * _s;
	_zo = zmin - 5 * _s;

	std::cout<< "*** Number of cells = " << _n << std::endl;

	_boundary_cells = new bool[_n];
	_external_cells = new bool[_n];
	_invalid_cells = new bool[_n];

	for (size_t i = 0; i < _n; i++)
	{
		_boundary_cells[i] = false;
		_external_cells[i] = false;
		_invalid_cells[i] = true;
	}
	_num_valid = 0;
	for (size_t i = 3; i < _ni - 2; i++)
	{
		for (size_t j = 3; j < _nj - 3; j++)
		{
			for (size_t k = 3; k < _nk - 3; k++)
			{
				size_t icell;
				get_cell_index(i, j, k, icell);
				_invalid_cells[icell] = false;
				_num_valid++;
			}
		}
	}

	#pragma endregion
}

int MeshingGenie_3d::sprinkle_points_along_sharp_edges(double d)
{
	#pragma region Sprinkle Points Along Boundaries:
	double max_dist_sq = 4 * d * d;

	// form loop
	size_t num_sharp_edges(_sharp_edges.size() / 2);
	std::vector<bool> taken(num_sharp_edges);
	for (size_t ied = 0; ied < num_sharp_edges; ied++) taken[ied] = false;
	std::list<size_t> loop;
	std::list<size_t>::iterator iter;
	std::vector<size_t> loop_v;
	std::vector<double> loop_d;
	
	while (true)
	{
		loop.clear(); size_t j1(0), j2(0);
		for (size_t ied = 0; ied < num_sharp_edges; ied++)
		{        
			if (taken[ied]) continue;

			size_t i1(_sharp_edges[2 * ied]), i2(_sharp_edges[2 * ied + 1]);
			taken[ied] = true;
			j1 = i1; j2 = i2;
			loop.push_back(j1); loop.push_back(j2);
			break;			
		}

		if (j1 == j2) break;

		for (size_t ied = 0; ied < num_sharp_edges; ied++)
		{
			if (taken[ied]) continue;

			size_t i1(_sharp_edges[2 * ied]), i2(_sharp_edges[2 * ied + 1]);
			if (i1 == j2 && i2 != j1)
			{
				loop.push_back(i2); j2 = i2; taken[ied] = true;			
			}
			else if (i2 == j2 && i1 != j1)
			{
				loop.push_back(i1); j2 = i2; taken[ied] = true;
			}
			else if (i1 == j1 && i2 != j2)
			{
				loop.push_front(i2); j1 = i2; taken[ied] = true;
			}
			else if (i2 == j1 && i1 != j2)
			{
				loop.push_front(i1); j1 = i1; taken[ied] = true;
			}
		}

		loop_v.clear(); loop_v.reserve(loop.size());
		for (iter = loop.begin(); iter!=loop.end(); iter++)
		{
			loop_v.push_back(*iter);
		}
		loop_d.resize(loop.size());
	
		size_t num_nodes(loop_v.size() - 1); double L(0.0);
		size_t inode_st(0); loop_d[0] = 0.0;
		for (size_t inode = 0; inode < num_nodes; inode++)
		{
			size_t i1 = loop_v[inode]; size_t i2 = loop_v[inode + 1];
			double xo(_xb[i1]), yo(_yb[i1]), zo(_zb[i1]);
			double xn(_xb[i2]), yn(_yb[i2]), zn(_zb[i2]);
			double dx(xn - xo), dy(yn - yo), dz(zn - zo);
			L += sqrt(dx * dx + dy * dy + dz * dz);
			loop_d[inode + 1] = L;
		}

		double px(_xb[loop_v[0]]), py(_yb[loop_v[0]]), pz(_zb[loop_v[0]]);
		double qx(_xb[loop_v[num_nodes]]), qy(_yb[loop_v[num_nodes]]), qz(_zb[loop_v[num_nodes]]);

		EdgePoint* po = new EdgePoint(); EdgePoint* pn = new EdgePoint();
		po->x = px; po->y = py; po->z = pz; po->inode = 0; po->h = L; po->v = 0.0;
		pn->x = qx; pn->y = qy; pn->z = qz; pn->inode = num_nodes - 1; pn->h = 0.0; pn->v = L - loop_d[num_nodes - 1];  
		po->next = pn; pn->next = 0;

		EdgePoint* pp(po);
		size_t num_points(2);
		while (pp != 0)
		{
			if (pp->h <= 2 * d)
			{
				pp = pp->next;
				continue;				
			}

			double ds = pp->h - 2 * d;
			double u = generate_a_random_number();
			double tmpdst = d + u * ds;
			
			double dst(- pp->v);
			for (size_t inode = pp->inode; inode < num_nodes; inode++)
			{
				size_t i1 = loop_v[inode]; size_t i2 = loop_v[inode + 1];
				dst += loop_d[inode + 1] - loop_d[inode];
				if (dst > tmpdst)
				{
					double xo(_xb[i1]), yo(_yb[i1]), zo(_zb[i1]);
					double xn(_xb[i2]), yn(_yb[i2]), zn(_zb[i2]);
					double dx(xn - xo), dy(yn - yo), dz(zn - zo);

					double r = (dst - tmpdst) / (loop_d[inode + 1] - loop_d[inode]);
			
					double xx = xo + r * (xn - xo);
					double yy = yo + r * (yn - yo);
					double zz = zo + r * (zn - zo);
					
					pn = new EdgePoint();
					pn->inode = inode;
					pn->h = pp->h - tmpdst; pp->h = tmpdst;
					pn->v = (loop_d[inode + 1] - loop_d[inode]) - (dst - tmpdst);
					pn->x = xx;
					pn->y = yy;
					pn->z = zz;
					pn->next = pp->next;
					num_points++;

					pp->next = pn; 

					if (tmpdst <= 2 * d) pp = pn; 

					break;
				}
			}
		}

		// sprinke points along discrete random edge
		pp = po; Point* newpoint;
		while (pp != 0)
		{
			size_t i, j, k, icell;
			i = size_t((pp->x - _xo) / _s);
			j = size_t((pp->y - _yo) / _s);
			k = size_t((pp->z - _zo) / _s);
			get_cell_index(i, j, k, icell);
			
			newpoint = closest_point(pp->x, pp->y, pp->z, 1E-10);
			if (newpoint == 0)
			{
				newpoint = new Point();
				newpoint->x = pp->x;
				newpoint->y = pp->y;
				newpoint->z = pp->z;
				newpoint->_on_reflex_corner = true;
				newpoint->_on_internal_boundaries = false;
				newpoint->next = 0;
				newpoint->poly = 0;

				if (_cell_points[icell] == 0) _cell_points[icell] = newpoint;
				else
				{
					Point* q = _cell_points[icell];
					while (q->next != 0) q = q->next;
					q->next = newpoint;
				}
					
				_num_inserted_points++;

				if (!_invalid_cells[icell]) 
				{
					_invalid_cells[icell] = true; _num_valid--;
				}						
			}
			po = pp;
			pp = pp->next;
			delete po;
		}
	}
	return 0;
	#pragma endregion
}

int MeshingGenie_3d::sample_edge(double xo, double yo, double zo, double xn, double yn, double zn, double d,
	                             std::vector<double> &ex, std::vector<double> &ey, std::vector<double> &ez)
{
	#pragma region Sprinkle Points Along Boundaries:
	double max_dist_sq = 4 * d * d;
	
	double L = sqrt((xn - xo) * (xn - xo) + (yn - yo) * (yn - yo) + (zn - zo) * (zn - zo));
	EdgePoint* po = new EdgePoint(); EdgePoint* pn = new EdgePoint();
	po->x = xo; po->y = yo; po->z = zo; po->h = L;
	pn->x = xn; pn->y = yn; pn->z = zn; pn->h = 0.0;
	po->next = pn; pn->next = 0;

	EdgePoint* pp(po);
	size_t num_points(2);
	while (pp->next != 0)
	{
		if (pp->h <= 2 * d)
		{
			pp = pp->next;
			continue;				
		}

		double ds = pp->h - 2 * d;
		double u = generate_a_random_number();
		double tmpdst = d + u * ds;
			
		double r = tmpdst / pp->h;
		double xx = pp->x + r * (pp->next->x - pp->x);
		double yy = pp->y + r * (pp->next->y - pp->y);
		double zz = pp->z + r * (pp->next->z - pp->z);

		EdgePoint* newpoint = new EdgePoint();
		newpoint->x = xx; newpoint->y = yy; newpoint->z = zz;
		newpoint->h = pp->h - tmpdst; pp->h = tmpdst;
		newpoint->next = pp->next;
		num_points++;

		pp->next = newpoint; pp->h = tmpdst;
	}

	ex.clear(); ey.clear(); ez.clear();
	ex.reserve(num_points);
	ey.reserve(num_points);
	ez.reserve(num_points);
	pp = po;
	while (pp != 0)
	{
		ex.push_back(pp->x);
		ey.push_back(pp->y);
		ez.push_back(pp->z);
		pp = pp->next;
	}
	return 0;
	#pragma endregion
}

void MeshingGenie_3d::identify_external_cell()
{
	#pragma region Identify External Cells:
	// loop over faces and sample points using 0.01 r spacing
	size_t num_faces = _bfaces.size();
	for (size_t iface = 0; iface < num_faces; iface++)
	{
		#pragma region Sample faces:
		size_t num_face_nodes(_bfaces[iface].size());
		size_t io(_bfaces[iface][0]), i1(_bfaces[iface][1]), i2(_bfaces[iface][2]);
		double xo(_xb[io]), yo(_yb[io]), zo(_zb[io]);		
		double x1(_xb[i1]), y1(_yb[i1]), z1(_zb[i1]);

		for (size_t inode = 2; inode < num_face_nodes; inode++)
		{
			i2 = _bfaces[iface][inode];
			double x2(_xb[i2]), y2(_yb[i2]), z2(_zb[i2]);

			// sample face
			double dx(x2 - x1), dy(y2 - y1), dz(z2 - z1);
			double L(sqrt(dx * dx + dy * dy + dz * dz));
			size_t m(size_t(L / (0.05 *_dm))); double s = L / m;

			double xx(x1), yy(y1), zz(z1);
			dx /= m; dy /= m; dz /= m;
			for (size_t ipnt = 0; ipnt <= m; ipnt++)
			{
				double dx_(xx - xo), dy_(yy-yo), dz_(zz-zo);
				double L_(sqrt(dx_ * dx_ + dy_ * dy_ + dz_ * dz_));
				size_t m_(size_t(L_/(0.05 *_dm))); double s_ = L_ / m_;
				double xx_(xo), yy_(yo), zz_(zo);
				dx_ /= m_; dy_ /= m_; dz_ /= m_;
				for (size_t jpnt = 0; jpnt <= m_; jpnt++)
				{
					// sample point
					size_t i, j, k, icell;
					
					i = size_t((xx_ - _xo) / _s);
					j = size_t((yy_ - _yo) / _s);
					k = size_t((zz_ - _zo) / _s);
					get_cell_index(i, j, k, icell);
					_boundary_cells[icell] = true;

					_cell_faces_iter = _cell_faces.find(icell);
					if (_cell_faces_iter == _cell_faces.end())
					{
						std::set<size_t> newset;
						newset.insert(iface);
						_cell_faces[icell] = newset;
					}
					else _cell_faces_iter->second.insert(iface);

					if (fabs(xx_ -_xo - i * _s) < 1E-10) 
					{
						get_cell_index(i - 1, j, k, icell);
						_boundary_cells[icell] = true;
						_cell_faces_iter = _cell_faces.find(icell);
						if (_cell_faces_iter == _cell_faces.end())
						{
							std::set<size_t> newset;
							newset.insert(iface);
							_cell_faces[icell] = newset;
						}
						else _cell_faces_iter->second.insert(iface);
					}
					if (fabs(xx_ -_xo - (i + 1) * _s) < 1E-10) 
					{
						get_cell_index(i + 1, j, k, icell);
						_boundary_cells[icell] = true;

						_cell_faces_iter = _cell_faces.find(icell);
						if (_cell_faces_iter == _cell_faces.end())
						{
							std::set<size_t> newset;
							newset.insert(iface);
							_cell_faces[icell] = newset;
						}
						else _cell_faces_iter->second.insert(iface);
					}

					if (fabs(yy_ -_yo - j * _s) < 1E-10) 
					{
						get_cell_index(i, j - 1, k, icell);
						_boundary_cells[icell] = true;

						_cell_faces_iter = _cell_faces.find(icell);
						if (_cell_faces_iter == _cell_faces.end())
						{
							std::set<size_t> newset;
							newset.insert(iface);
							_cell_faces[icell] = newset;
						}
						else _cell_faces_iter->second.insert(iface);
					}
					if (fabs(yy_ -_yo - (j + 1) * _s) < 1E-10) 
					{
						get_cell_index(i, j + 1, k, icell);
						_boundary_cells[icell] = true;

						_cell_faces_iter = _cell_faces.find(icell);
						if (_cell_faces_iter == _cell_faces.end())
						{
							std::set<size_t> newset;
							newset.insert(iface);
							_cell_faces[icell] = newset;
						}
						else _cell_faces_iter->second.insert(iface);
					}

					if (fabs(zz_ -_zo - k * _s) < 1E-10) 
					{
						get_cell_index(i, j, k - 1, icell);
						_boundary_cells[icell] = true;

						_cell_faces_iter = _cell_faces.find(icell);
						if (_cell_faces_iter == _cell_faces.end())
						{
							std::set<size_t> newset;
							newset.insert(iface);
							_cell_faces[icell] = newset;
						}
						else _cell_faces_iter->second.insert(iface);
					}

					if (fabs(zz_ -_zo - (k + 1) * _s) < 1E-10) 
					{
						get_cell_index(i, j, k + 1, icell);
						_boundary_cells[icell] = true;

						_cell_faces_iter = _cell_faces.find(icell);
						if (_cell_faces_iter == _cell_faces.end())
						{
							std::set<size_t> newset;
							newset.insert(iface);
							_cell_faces[icell] = newset;
						}
						else _cell_faces_iter->second.insert(iface);
					}

					xx_+= dx_; yy_+= dy_; zz_+= dz_;
				}
				xx+= dx; yy+= dy; zz+= dz;
			}
			// move to the next triangle
			x1 = x2; y1 = y2; z1 = z2;			
		}
		#pragma endregion
	}

	size_t i, j, k, jcell;

	std::set<size_t>* faces;
	std::set<size_t>::iterator face_iter;
	for (size_t icell = 0; icell < _n; icell++)
	{
		#pragma region Identify Some External cells:

		if (!_boundary_cells[icell]) continue;

		faces = (&_cell_faces[icell]);

		bool sharp_feature(false);
		for (face_iter = faces->begin(); face_iter != faces->end(); face_iter++)
		{
			size_t fb(*face_iter);
			if (_nonconvex_bfaces[fb]) {sharp_feature = true; break;}
		}

		if (sharp_feature) continue;
		
		size_t fa(*faces->begin());
		double ax(_bfaces_nx[fa]);
		double ay(_bfaces_ny[fa]);
		double az(_bfaces_nz[fa]);

		if (fabs(ax) > fabs(ay) &&  fabs(ax) > fabs(az) && ax > 0.0) 
		{
			get_cell_indices(icell, i, j, k);
			get_cell_index(i + 1, j, k, jcell);
			if (!_boundary_cells[jcell]) _external_cells[jcell] = true;
		}
			
		if (fabs(ax) > fabs(ay) &&  fabs(ax) > fabs(az) && ax < 0.0) 
		{
			get_cell_indices(icell, i, j, k);
			get_cell_index(i - 1, j, k, jcell);
			if (!_boundary_cells[jcell]) _external_cells[jcell] = true;
		}

		if (fabs(ay) > fabs(az) &&  fabs(ay) > fabs(ax) && ay > 0.0) 
		{
			get_cell_indices(icell, i, j, k);
			get_cell_index(i, j + 1, k, jcell);
			if (!_boundary_cells[jcell]) _external_cells[jcell] = true;
		}
			
		if (fabs(ay) > fabs(az) &&  fabs(ay) > fabs(ax) && ay < 0.0) 
		{
			get_cell_indices(icell, i, j, k);
			get_cell_index(i, j - 1, k, jcell);
			if (!_boundary_cells[jcell]) _external_cells[jcell] = true;
		}
		
		if (fabs(az) > fabs(ax) &&  fabs(az) > fabs(ay) && az > 0.0) 
		{
			get_cell_indices(icell, i, j, k);
			get_cell_index(i, j, k + 1, jcell);
			if (!_boundary_cells[jcell]) _external_cells[jcell] = true;
		}
			
		if (fabs(az) > fabs(ax) &&  fabs(az) > fabs(ay) && az < 0.0) 
		{
			get_cell_indices(icell, i, j, k);
			get_cell_index(i, j, k - 1, jcell);
			if (!_boundary_cells[jcell]) _external_cells[jcell] = true;
		}

		#pragma endregion
	}

	while (true)
	{
		#pragma region Apply flooding:
		bool nothing_new(true);
		for (size_t icell = 0; icell < _n; icell++) 
		{
			if (!_external_cells[icell]) continue;

			if (_invalid_cells[icell]) continue;

			get_cell_indices(icell, i, j, k);

			get_cell_index(i + 1, j, k, jcell);
			if (!_boundary_cells[jcell] && !_external_cells[jcell]) {_external_cells[jcell] = true; nothing_new = false;}

			get_cell_index(i - 1, j, k, jcell);
			if (!_boundary_cells[jcell] && !_external_cells[jcell]) {_external_cells[jcell] = true; nothing_new = false;}
		
			get_cell_index(i, j + 1, k, jcell);
			if (!_boundary_cells[jcell] && !_external_cells[jcell]) {_external_cells[jcell] = true; nothing_new = false;}

			get_cell_index(i, j - 1, k, jcell);
			if (!_boundary_cells[jcell] && !_external_cells[jcell]) {_external_cells[jcell] = true; nothing_new = false;}

			get_cell_index(i, j, k + 1, jcell);
			if (!_boundary_cells[jcell] && !_external_cells[jcell]) {_external_cells[jcell] = true; nothing_new = false;}

			get_cell_index(i, j, k - 1, jcell);
			if (!_boundary_cells[jcell] && !_external_cells[jcell]) {_external_cells[jcell] = true; nothing_new = false;}
		}
		if (nothing_new) break;
		#pragma endregion
	}
	#pragma endregion
}

void MeshingGenie_3d::use_grid_level(size_t refLevel)
{
	#pragma region Grid Level refLevel:
	
	clock_t start_time, end_time; double cpu_time;	 
	start_time = clock();	

	size_t icell, jcell;
	size_t no(0); size_t num_inserted(_num_inserted_points);
	if (refLevel == 0)
	{
		// create an active pool
		_pool_size = _num_valid * 3;
		_active_pool_i = new size_t[_pool_size];
		_active_pool_j = new size_t[_pool_size];
		_active_pool_k = new size_t[_pool_size];
		_active_pool_icell = new size_t[_pool_size];
		size_t io, jo, ko;
		for (size_t icell = 0; icell < _n; icell++) 
		{
			if (!_invalid_cells[icell])			
			{
				get_cell_indices(icell, io, jo, ko);
				_active_pool_i[no] = io;
				_active_pool_j[no] = jo; 
				_active_pool_k[no] = ko;
				_active_pool_icell[no] = icell;
				no++;
			}			
		}
		_no = no;
		_TWO_POW_REF_LEVEL = 1;
		_ss = _s;
	}
	else
	{
		size_t new_two_pow_ref_lef  = _TWO_POW_REF_LEVEL * 2;
		
		bool c1(false), c2(false), c3(false), c4(false);
		bool c5(false), c6(false), c7(false), c8(false);
		
		size_t i, j, k;
		size_t ii(_pool_size - 1), kk(_num_valid - 1);
		size_t prevLevel(refLevel - 1);
		while (true)
		{
			i = _active_pool_i[kk]; j = _active_pool_j[kk]; k = _active_pool_k[kk]; icell = _active_pool_icell[kk];
			if (!_invalid_cells[icell])
			{
				get_active_children(i, j, k, icell, c1, c2, c3, c4, c5, c6, c7, c8);
				i *= 2; j *= 2; k *= 2;
				if (!c1)
				{					
					_active_pool_i[ii] = i; 
					_active_pool_j[ii] = j;
					_active_pool_k[ii] = k;
					_active_pool_icell[ii] = icell;
					ii--;
				}
				if (!c2)
				{
					_active_pool_i[ii] = i; 
					_active_pool_j[ii] = j;
					_active_pool_k[ii] = k + 1;
					_active_pool_icell[ii] = icell;
					ii--;
				}
				if (!c3) 
				{
					_active_pool_i[ii] = i; 
					_active_pool_j[ii] = j + 1;
					_active_pool_k[ii] = k;
					_active_pool_icell[ii] = icell;
					ii--;
				}
				if (!c4) 
				{
					_active_pool_i[ii] = i; 
					_active_pool_j[ii] = j + 1;
					_active_pool_k[ii] = k + 1;
					_active_pool_icell[ii] = icell;
					ii--;
				}
				if (!c5)
				{					
					_active_pool_i[ii] = i + 1; 
					_active_pool_j[ii] = j;
					_active_pool_k[ii] = k;
					_active_pool_icell[ii] = icell;
					ii--;
				}
				if (!c6)
				{
					_active_pool_i[ii] = i + 1; 
					_active_pool_j[ii] = j;
					_active_pool_k[ii] = k + 1;
					_active_pool_icell[ii] = icell;
					ii--;
				}
				if (!c7) 
				{
					_active_pool_i[ii] = i + 1; 
					_active_pool_j[ii] = j + 1;
					_active_pool_k[ii] = k;
					_active_pool_icell[ii] = icell;
					ii--;
				}
				if (!c8) 
				{
					_active_pool_i[ii] = i + 1; 
					_active_pool_j[ii] = j + 1;
					_active_pool_k[ii] = k + 1;
					_active_pool_icell[ii] = icell;
					ii--;
				}
			}
			if (kk == 0) break;
			
			kk--;

			if (ii < kk)
			{
				std::cout << "ERROR in filling active pool" << std::endl;
			}
		}
		ii++;
		_num_valid = _pool_size - ii; no = _num_valid;
		for (size_t i = 0; i < _num_valid; i++)
		{
			_active_pool_i[i] = _active_pool_i[i + ii];
			_active_pool_j[i] = _active_pool_j[i + ii];
			_active_pool_k[i] = _active_pool_k[i + ii];
			_active_pool_icell[i] = _active_pool_icell[i + ii];
		}
		_TWO_POW_REF_LEVEL *= 2;
		_ss *= 0.5;
	}

	if (_num_valid == 0) 
	{
		end_time = clock();
		cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
	
		std::cout<< "Refinement Level : (" << refLevel << ")" << std::endl;
		std::cout<< "spacing : (" << _ss << ")" << std::endl;
		std::cout<< "Pool is empty .. maximal condition is achieved!" << std::endl;
		std::cout<< "Execution Time = " << cpu_time << " seconds." << std::endl;
		std::cout<< std::endl;
		return;
	}
	
	Point* pp;
	size_t i, j, k;
	size_t irand;	
	double u, v, w, xx, yy, zz;
	size_t num_invalidated(0), num_invalidated_invalid(0);

	size_t nit(0);
	size_t nit_max(size_t(_RF * no));
	while (nit < nit_max)	
	{
		// pick a random cell from active pool:
		u = generate_a_random_number();
		irand = size_t(u * (_num_valid - 1));

		i = _active_pool_i[irand];
		j = _active_pool_j[irand];
		k = _active_pool_k[irand];
		icell = _active_pool_icell[irand];
		
		if (_invalid_cells[icell])
		{
			_num_valid--;			
			size_t tmp =  _active_pool_i[irand];
			_active_pool_i[irand] = _active_pool_i[_num_valid];
			_active_pool_i[_num_valid] = tmp;
			tmp =  _active_pool_j[irand];
			_active_pool_j[irand] = _active_pool_j[_num_valid];
			_active_pool_j[_num_valid] = tmp;
			tmp =  _active_pool_k[irand];
			_active_pool_k[irand] = _active_pool_k[_num_valid];
			_active_pool_k[_num_valid] = tmp;
			tmp =  _active_pool_icell[irand];
			_active_pool_icell[irand] = _active_pool_icell[_num_valid];
			_active_pool_icell[_num_valid] = tmp;
			num_invalidated_invalid++;
			nit++;
			if (_num_valid == 0) break;
			continue;
		}

		// pick a random point
		u = generate_a_random_number();
		v = generate_a_random_number();
		w = generate_a_random_number();
		
		xx = _xo + (i + u) * _ss;
		yy = _yo + (j + v) * _ss;
		zz = _zo + (k + w) * _ss;

		// check for conflicts
		bool valid_point(true);
		for (size_t index = 0; index < 124; index++) 
		{
			jcell = icell + _neighbors[index];
			pp = _cell_points[jcell];
			if (pp == 0) continue;

			double dd = distance_squared(pp->x - xx, pp->y - yy, pp->z - zz);
			if (dd < _dm_squared)
			{
				valid_point = false;
				break;
			}
		}

		if (valid_point)
		{
			// Create a new points and assign it to that cell
			_invalid_cells[icell] = true;

			_cell_points[icell] = new Point();
			_cell_points[icell]->x = xx;
			_cell_points[icell]->y = yy;
			_cell_points[icell]->z = zz;
			_cell_points[icell]->_on_reflex_corner = false;
			_cell_points[icell]->_on_internal_boundaries = false;
			_cell_points[icell]->next = 0;
			_cell_points[icell]->poly = 0;
			_num_inserted_points++;

			// deactivae that cell
			_num_valid--;
			size_t tmp =  _active_pool_i[irand];
			_active_pool_i[irand] = _active_pool_i[_num_valid];
			_active_pool_i[_num_valid] = tmp;
			tmp =  _active_pool_j[irand];
			_active_pool_j[irand] = _active_pool_j[_num_valid];
			_active_pool_j[_num_valid] = tmp;
			tmp =  _active_pool_k[irand];
			_active_pool_k[irand] = _active_pool_k[_num_valid];
			_active_pool_k[_num_valid] = tmp;
			tmp =  _active_pool_icell[irand];
			_active_pool_icell[irand] = _active_pool_icell[_num_valid];
			_active_pool_icell[_num_valid] = tmp;
						
			if (_num_valid == 0) break;
		}
		nit++;
	}

	end_time = clock();
	cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
	
	std::cout<< "Refinement Level : (" << refLevel << ")" << std::endl;
	std::cout<< "Active cells ratio = " << no * 1.0 / _pool_size << std::endl;
	std::cout<< "Number of inserted points = " << (_num_inserted_points - num_inserted) << std::endl;
	std::cout<< "Execution Time = " << cpu_time << " seconds." << std::endl;
	std::cout<< std::endl;
	std::cout<< "Iterations ratio = " << nit * 1.0 / no << std::endl;
	std::cout<< "Remaining Cells Ratio = " << _num_valid * 1.0 / no << std::endl;
	std::cout<< "Invalidated Ratio due to point insertion = " << num_invalidated * 1.0 / no << std::endl;
	std::cout<< "Invalidated Ratio due to invalid cells = " << num_invalidated_invalid * 1.0 / no << std::endl;
	std::cout<< "Insertion Rate = " << (_num_inserted_points - num_inserted) * 1.0 / cpu_time << " points/sec." << std::endl;
	std::cout<< "--------------------------------------------------------------------" << std::endl;
	std::cout<< std::endl;
	#pragma endregion
}

inline double MeshingGenie_3d::distance_squared(double dx, double dy, double dz)
{
    return dx * dx + dy * dy + dz * dz;        
}


inline void MeshingGenie_3d::get_neighbor_indices(size_t index, size_t io, size_t jo, size_t ko, size_t &ii, size_t &jj, size_t &kk)
{
	#pragma region get neighbor indices:

	ii = 0; jj = 0; kk = 0;
	// sharing a face
	if (index ==   0) {ii = io + 1; jj = jo    ; kk = ko    ; return;}
	if (index ==   1) {ii = io - 1; jj = jo    ; kk = ko    ; return;}
	if (index ==   2) {ii = io    ; jj = jo + 1; kk = ko    ; return;}
	if (index ==   3) {ii = io    ; jj = jo - 1; kk = ko    ; return;}	
	if (index ==   4) {ii = io    ; jj = jo    ; kk = ko + 1; return;}
	if (index ==   5) {ii = io    ; jj = jo    ; kk = ko - 1; return;}	

	// sharing an edge 
	if (index ==   6) {ii = io    ; jj = jo + 1; kk = ko + 1; return;}
	if (index ==   7) {ii = io    ; jj = jo + 1; kk = ko - 1; return;}
	if (index ==   8) {ii = io    ; jj = jo - 1; kk = ko + 1; return;}
	if (index ==   9) {ii = io    ; jj = jo - 1; kk = ko - 1; return;}

	if (index ==  10) {ii = io + 1; jj = jo    ; kk = ko + 1; return;}
	if (index ==  11) {ii = io + 1; jj = jo    ; kk = ko - 1; return;}
	if (index ==  12) {ii = io - 1; jj = jo    ; kk = ko + 1; return;}
	if (index ==  13) {ii = io - 1; jj = jo    ; kk = ko - 1; return;}

	if (index ==  14) {ii = io + 1; jj = jo + 1; kk = ko    ; return;}
	if (index ==  15) {ii = io + 1; jj = jo - 1; kk = ko    ; return;}
	if (index ==  16) {ii = io - 1; jj = jo + 1; kk = ko    ; return;}
	if (index ==  17) {ii = io - 1; jj = jo - 1; kk = ko    ; return;}

	// sharing a vertex
	if (index ==  18) {ii = io + 1; jj = jo + 1; kk = ko + 1; return;}
	if (index ==  19) {ii = io + 1; jj = jo + 1; kk = ko - 1; return;}
	if (index ==  20) {ii = io + 1; jj = jo - 1; kk = ko + 1; return;}
	if (index ==  21) {ii = io + 1; jj = jo - 1; kk = ko - 1; return;}

	if (index ==  22) {ii = io - 1; jj = jo + 1; kk = ko + 1; return;}
	if (index ==  23) {ii = io - 1; jj = jo + 1; kk = ko - 1; return;}
	if (index ==  24) {ii = io - 1; jj = jo - 1; kk = ko + 1; return;}
	if (index ==  25) {ii = io - 1; jj = jo - 1; kk = ko - 1; return;}

	// Second Layer :

	if (index ==  26) {ii = io    ; jj = jo    ; kk = ko + 2; return;} // Up
	if (index ==  27) {ii = io    ; jj = jo + 1; kk = ko + 2; return;}
	if (index ==  28) {ii = io    ; jj = jo + 2; kk = ko + 2; return;} // Left
	if (index ==  29) {ii = io    ; jj = jo + 2; kk = ko + 1; return;}
	if (index ==  30) {ii = io    ; jj = jo + 2; kk = ko    ; return;}
	if (index ==  31) {ii = io    ; jj = jo + 2; kk = ko - 1; return;}
	if (index ==  32) {ii = io    ; jj = jo + 2; kk = ko - 2; return;} // down
	if (index ==  33) {ii = io    ; jj = jo + 1; kk = ko - 2; return;}
	if (index ==  34) {ii = io    ; jj = jo    ; kk = ko - 2; return;}
	if (index ==  35) {ii = io    ; jj = jo - 1; kk = ko - 2; return;}
	if (index ==  36) {ii = io    ; jj = jo - 2; kk = ko - 2; return;} // right
	if (index ==  37) {ii = io    ; jj = jo - 2; kk = ko - 1; return;}
	if (index ==  38) {ii = io    ; jj = jo - 2; kk = ko    ; return;}
	if (index ==  39) {ii = io    ; jj = jo - 2; kk = ko + 1; return;}
	if (index ==  40) {ii = io    ; jj = jo - 2; kk = ko + 2; return;} // up
	if (index ==  41) {ii = io    ; jj = jo - 1; kk = ko + 2; return;}

	if (index ==  42) {ii = io + 1; jj = jo    ; kk = ko + 2; return;} // Up
	if (index ==  43) {ii = io + 1; jj = jo + 1; kk = ko + 2; return;}
	if (index ==  44) {ii = io + 1; jj = jo + 2; kk = ko + 2; return;} // Left
	if (index ==  45) {ii = io + 1; jj = jo + 2; kk = ko + 1; return;}
	if (index ==  46) {ii = io + 1; jj = jo + 2; kk = ko    ; return;}
	if (index ==  47) {ii = io + 1; jj = jo + 2; kk = ko - 1; return;}
	if (index ==  48) {ii = io + 1; jj = jo + 2; kk = ko - 2; return;} // down
	if (index ==  49) {ii = io + 1; jj = jo + 1; kk = ko - 2; return;}
	if (index ==  50) {ii = io + 1; jj = jo    ; kk = ko - 2; return;}
	if (index ==  51) {ii = io + 1; jj = jo - 1; kk = ko - 2; return;}
	if (index ==  52) {ii = io + 1; jj = jo - 2; kk = ko - 2; return;} // right
	if (index ==  53) {ii = io + 1; jj = jo - 2; kk = ko - 1; return;}
	if (index ==  54) {ii = io + 1; jj = jo - 2; kk = ko    ; return;}
	if (index ==  55) {ii = io + 1; jj = jo - 2; kk = ko + 1; return;}
	if (index ==  56) {ii = io + 1; jj = jo - 2; kk = ko + 2; return;} // up
	if (index ==  57) {ii = io + 1; jj = jo - 1; kk = ko + 2; return;}

	if (index ==  58) {ii = io - 1; jj = jo    ; kk = ko + 2; return;} // Up
	if (index ==  59) {ii = io - 1; jj = jo + 1; kk = ko + 2; return;}
	if (index ==  60) {ii = io - 1; jj = jo + 2; kk = ko + 2; return;} // Left
	if (index ==  61) {ii = io - 1; jj = jo + 2; kk = ko + 1; return;}
	if (index ==  62) {ii = io - 1; jj = jo + 2; kk = ko    ; return;}
	if (index ==  63) {ii = io - 1; jj = jo + 2; kk = ko - 1; return;}
	if (index ==  64) {ii = io - 1; jj = jo + 2; kk = ko - 2; return;} // down
	if (index ==  65) {ii = io - 1; jj = jo + 1; kk = ko - 2; return;}
	if (index ==  66) {ii = io - 1; jj = jo    ; kk = ko - 2; return;}
	if (index ==  67) {ii = io - 1; jj = jo - 1; kk = ko - 2; return;}
	if (index ==  68) {ii = io - 1; jj = jo - 2; kk = ko - 2; return;} // right
	if (index ==  69) {ii = io - 1; jj = jo - 2; kk = ko - 1; return;}
	if (index ==  70) {ii = io - 1; jj = jo - 2; kk = ko    ; return;}
	if (index ==  71) {ii = io - 1; jj = jo - 2; kk = ko + 1; return;}
	if (index ==  72) {ii = io - 1; jj = jo - 2; kk = ko + 2; return;} // up
	if (index ==  73) {ii = io - 1; jj = jo - 1; kk = ko + 2; return;}

	if (index ==  74) {ii = io + 2; jj = jo    ; kk = ko    ; return;} // right
	if (index ==  75) {ii = io + 2; jj = jo    ; kk = ko + 1; return;} // up
	if (index ==  76) {ii = io + 2; jj = jo + 1; kk = ko + 1; return;} // Left
	if (index ==  77) {ii = io + 2; jj = jo + 1; kk = ko    ; return;}
	if (index ==  78) {ii = io + 2; jj = jo + 1; kk = ko - 1; return;} // down
	if (index ==  79) {ii = io + 2; jj = jo    ; kk = ko - 1; return;}
	if (index ==  80) {ii = io + 2; jj = jo - 1; kk = ko - 1; return;} // right
	if (index ==  81) {ii = io + 2; jj = jo - 1; kk = ko    ; return;}
	if (index ==  82) {ii = io + 2; jj = jo - 1; kk = ko + 1; return;}

	if (index ==  83) {ii = io - 2; jj = jo    ; kk = ko    ; return;} // right
	if (index ==  84) {ii = io - 2; jj = jo    ; kk = ko + 1; return;} // up
	if (index ==  85) {ii = io - 2; jj = jo + 1; kk = ko + 1; return;} // Left
	if (index ==  86) {ii = io - 2; jj = jo + 1; kk = ko    ; return;}
	if (index ==  87) {ii = io - 2; jj = jo + 1; kk = ko - 1; return;} // down
	if (index ==  88) {ii = io - 2; jj = jo    ; kk = ko - 1; return;}
	if (index ==  89) {ii = io - 2; jj = jo - 1; kk = ko - 1; return;} // right
	if (index ==  90) {ii = io - 2; jj = jo - 1; kk = ko    ; return;}
	if (index ==  91) {ii = io - 2; jj = jo - 1; kk = ko + 1; return;}

	if (index ==  92) {ii = io + 2; jj = jo    ; kk = ko + 2; return;} // Up
	if (index ==  93) {ii = io + 2; jj = jo + 1; kk = ko + 2; return;}
	if (index ==  94) {ii = io + 2; jj = jo + 2; kk = ko + 2; return;} // Left
	if (index ==  95) {ii = io + 2; jj = jo + 2; kk = ko + 1; return;}
	if (index ==  96) {ii = io + 2; jj = jo + 2; kk = ko    ; return;}
	if (index ==  97) {ii = io + 2; jj = jo + 2; kk = ko - 1; return;}
	if (index ==  98) {ii = io + 2; jj = jo + 2; kk = ko - 2; return;} // down
	if (index ==  99) {ii = io + 2; jj = jo + 1; kk = ko - 2; return;}
	if (index == 100) {ii = io + 2; jj = jo    ; kk = ko - 2; return;}
	if (index == 101) {ii = io + 2; jj = jo - 1; kk = ko - 2; return;}
	if (index == 102) {ii = io + 2; jj = jo - 2; kk = ko - 2; return;} // right
	if (index == 103) {ii = io + 2; jj = jo - 2; kk = ko - 1; return;}
	if (index == 104) {ii = io + 2; jj = jo - 2; kk = ko    ; return;}
	if (index == 105) {ii = io + 2; jj = jo - 2; kk = ko + 1; return;}
	if (index == 106) {ii = io + 2; jj = jo - 2; kk = ko + 2; return;} // up
	if (index == 107) {ii = io + 2; jj = jo - 1; kk = ko + 2; return;}

	if (index == 108) {ii = io - 2; jj = jo    ; kk = ko + 2; return;} // Up
	if (index == 109) {ii = io - 2; jj = jo + 1; kk = ko + 2; return;}
	if (index == 110) {ii = io - 2; jj = jo + 2; kk = ko + 2; return;} // Left
	if (index == 111) {ii = io - 2; jj = jo + 2; kk = ko + 1; return;}
	if (index == 112) {ii = io - 2; jj = jo + 2; kk = ko    ; return;}
	if (index == 113) {ii = io - 2; jj = jo + 2; kk = ko - 1; return;}
	if (index == 114) {ii = io - 2; jj = jo + 2; kk = ko - 2; return;} // down
	if (index == 115) {ii = io - 2; jj = jo + 1; kk = ko - 2; return;}
	if (index == 116) {ii = io - 2; jj = jo    ; kk = ko - 2; return;}
	if (index == 117) {ii = io - 2; jj = jo - 1; kk = ko - 2; return;}
	if (index == 118) {ii = io - 2; jj = jo - 2; kk = ko - 2; return;} // right
	if (index == 119) {ii = io - 2; jj = jo - 2; kk = ko - 1; return;}
	if (index == 120) {ii = io - 2; jj = jo - 2; kk = ko    ; return;}
	if (index == 121) {ii = io - 2; jj = jo - 2; kk = ko + 1; return;}
	if (index == 122) {ii = io - 2; jj = jo - 2; kk = ko + 2; return;} // up
	if (index == 123) {ii = io - 2; jj = jo - 1; kk = ko + 2; return;}

	// Third Layer
	if (index == 124) {ii = io    ; jj = jo    ; kk = ko + 3; return;} // Up
	if (index == 125) {ii = io    ; jj = jo + 1; kk = ko + 3; return;} //
	if (index == 126) {ii = io    ; jj = jo + 2; kk = ko + 3; return;} //
	if (index == 127) {ii = io    ; jj = jo + 3; kk = ko + 3; return;} // Left
	if (index == 128) {ii = io    ; jj = jo + 3; kk = ko + 2; return;} // 
	if (index == 129) {ii = io    ; jj = jo + 3; kk = ko + 1; return;} //
	if (index == 130) {ii = io    ; jj = jo + 3; kk = ko    ; return;} //
	if (index == 131) {ii = io    ; jj = jo + 3; kk = ko - 1; return;} //
	if (index == 132) {ii = io    ; jj = jo + 3; kk = ko - 2; return;} //
	if (index == 133) {ii = io    ; jj = jo + 3; kk = ko - 3; return;} // down
	if (index == 134) {ii = io    ; jj = jo + 2; kk = ko - 3; return;} //
	if (index == 135) {ii = io    ; jj = jo + 1; kk = ko - 3; return;} //
	if (index == 136) {ii = io    ; jj = jo    ; kk = ko - 3; return;} //
	if (index == 137) {ii = io    ; jj = jo - 1; kk = ko - 3; return;} //
	if (index == 138) {ii = io    ; jj = jo - 2; kk = ko - 3; return;} //
	if (index == 139) {ii = io    ; jj = jo - 3; kk = ko - 3; return;} // Right
	if (index == 140) {ii = io    ; jj = jo - 3; kk = ko - 2; return;} //
	if (index == 141) {ii = io    ; jj = jo - 3; kk = ko - 1; return;} //
	if (index == 142) {ii = io    ; jj = jo - 3; kk = ko    ; return;} //
	if (index == 143) {ii = io    ; jj = jo - 3; kk = ko + 1; return;} //
	if (index == 144) {ii = io    ; jj = jo - 3; kk = ko + 2; return;} //
	if (index == 145) {ii = io    ; jj = jo - 3; kk = ko + 3; return;} // Up
	if (index == 146) {ii = io    ; jj = jo - 2; kk = ko + 3; return;} //
	if (index == 147) {ii = io    ; jj = jo - 1; kk = ko + 3; return;} //

	if (index == 148) {ii = io + 1; jj = jo    ; kk = ko + 3; return;} // Up
	if (index == 149) {ii = io + 1; jj = jo + 1; kk = ko + 3; return;} //
	if (index == 150) {ii = io + 1; jj = jo + 2; kk = ko + 3; return;} //
	if (index == 151) {ii = io + 1; jj = jo + 3; kk = ko + 3; return;} // Left
	if (index == 152) {ii = io + 1; jj = jo + 3; kk = ko + 2; return;} // 
	if (index == 153) {ii = io + 1; jj = jo + 3; kk = ko + 1; return;} //
	if (index == 154) {ii = io + 1; jj = jo + 3; kk = ko    ; return;} //
	if (index == 155) {ii = io + 1; jj = jo + 3; kk = ko - 1; return;} //
	if (index == 156) {ii = io + 1; jj = jo + 3; kk = ko - 2; return;} //
	if (index == 157) {ii = io + 1; jj = jo + 3; kk = ko - 3; return;} // down
	if (index == 158) {ii = io + 1; jj = jo + 2; kk = ko - 3; return;} //
	if (index == 159) {ii = io + 1; jj = jo + 1; kk = ko - 3; return;} //
	if (index == 160) {ii = io + 1; jj = jo    ; kk = ko - 3; return;} //
	if (index == 161) {ii = io + 1; jj = jo - 1; kk = ko - 3; return;} //
	if (index == 162) {ii = io + 1; jj = jo - 2; kk = ko - 3; return;} //
	if (index == 163) {ii = io + 1; jj = jo - 3; kk = ko - 3; return;} // Right
	if (index == 164) {ii = io + 1; jj = jo - 3; kk = ko - 2; return;} //
	if (index == 165) {ii = io + 1; jj = jo - 3; kk = ko - 1; return;} //
	if (index == 166) {ii = io + 1; jj = jo - 3; kk = ko    ; return;} //
	if (index == 167) {ii = io + 1; jj = jo - 3; kk = ko + 1; return;} //
	if (index == 168) {ii = io + 1; jj = jo - 3; kk = ko + 2; return;} //
	if (index == 169) {ii = io + 1; jj = jo - 3; kk = ko + 3; return;} // Up
	if (index == 170) {ii = io + 1; jj = jo - 2; kk = ko + 3; return;} //
	if (index == 171) {ii = io + 1; jj = jo - 1; kk = ko + 3; return;} //

	if (index == 172) {ii = io - 1; jj = jo    ; kk = ko + 3; return;} // Up
	if (index == 173) {ii = io - 1; jj = jo + 1; kk = ko + 3; return;} //
	if (index == 174) {ii = io - 1; jj = jo + 2; kk = ko + 3; return;} //
	if (index == 175) {ii = io - 1; jj = jo + 3; kk = ko + 3; return;} // Left
	if (index == 176) {ii = io - 1; jj = jo + 3; kk = ko + 2; return;} // 
	if (index == 177) {ii = io - 1; jj = jo + 3; kk = ko + 1; return;} //
	if (index == 178) {ii = io - 1; jj = jo + 3; kk = ko    ; return;} //
	if (index == 179) {ii = io - 1; jj = jo + 3; kk = ko - 1; return;} //
	if (index == 180) {ii = io - 1; jj = jo + 3; kk = ko - 2; return;} //
	if (index == 181) {ii = io - 1; jj = jo + 3; kk = ko - 3; return;} // down
	if (index == 182) {ii = io - 1; jj = jo + 2; kk = ko - 3; return;} //
	if (index == 183) {ii = io - 1; jj = jo + 1; kk = ko - 3; return;} //
	if (index == 184) {ii = io - 1; jj = jo    ; kk = ko - 3; return;} //
	if (index == 185) {ii = io - 1; jj = jo - 1; kk = ko - 3; return;} //
	if (index == 186) {ii = io - 1; jj = jo - 2; kk = ko - 3; return;} //
	if (index == 187) {ii = io - 1; jj = jo - 3; kk = ko - 3; return;} // Right
	if (index == 188) {ii = io - 1; jj = jo - 3; kk = ko - 2; return;} //
	if (index == 189) {ii = io - 1; jj = jo - 3; kk = ko - 1; return;} //
	if (index == 190) {ii = io - 1; jj = jo - 3; kk = ko    ; return;} //
	if (index == 191) {ii = io - 1; jj = jo - 3; kk = ko + 1; return;} //
	if (index == 192) {ii = io - 1; jj = jo - 3; kk = ko + 2; return;} //
	if (index == 193) {ii = io - 1; jj = jo - 3; kk = ko + 3; return;} // Up
	if (index == 194) {ii = io - 1; jj = jo - 2; kk = ko + 3; return;} //
	if (index == 195) {ii = io - 1; jj = jo - 1; kk = ko + 3; return;} //

	if (index == 196) {ii = io + 2; jj = jo    ; kk = ko + 3; return;} // Up
	if (index == 197) {ii = io + 2; jj = jo + 1; kk = ko + 3; return;} //
	if (index == 198) {ii = io + 2; jj = jo + 2; kk = ko + 3; return;} //
	if (index == 199) {ii = io + 2; jj = jo + 3; kk = ko + 3; return;} // Left
	if (index == 200) {ii = io + 2; jj = jo + 3; kk = ko + 2; return;} // 
	if (index == 201) {ii = io + 2; jj = jo + 3; kk = ko + 1; return;} //
	if (index == 202) {ii = io + 2; jj = jo + 3; kk = ko    ; return;} //
	if (index == 203) {ii = io + 2; jj = jo + 3; kk = ko - 1; return;} //
	if (index == 204) {ii = io + 2; jj = jo + 3; kk = ko - 2; return;} //
	if (index == 205) {ii = io + 2; jj = jo + 3; kk = ko - 3; return;} // down
	if (index == 206) {ii = io + 2; jj = jo + 2; kk = ko - 3; return;} //
	if (index == 207) {ii = io + 2; jj = jo + 1; kk = ko - 3; return;} //
	if (index == 208) {ii = io + 2; jj = jo    ; kk = ko - 3; return;} //
	if (index == 209) {ii = io + 2; jj = jo - 1; kk = ko - 3; return;} //
	if (index == 210) {ii = io + 2; jj = jo - 2; kk = ko - 3; return;} //
	if (index == 211) {ii = io + 2; jj = jo - 3; kk = ko - 3; return;} // Right
	if (index == 212) {ii = io + 2; jj = jo - 3; kk = ko - 2; return;} //
	if (index == 213) {ii = io + 2; jj = jo - 3; kk = ko - 1; return;} //
	if (index == 214) {ii = io + 2; jj = jo - 3; kk = ko    ; return;} //
	if (index == 215) {ii = io + 2; jj = jo - 3; kk = ko + 1; return;} //
	if (index == 216) {ii = io + 2; jj = jo - 3; kk = ko + 2; return;} //
	if (index == 217) {ii = io + 2; jj = jo - 3; kk = ko + 3; return;} // Up
	if (index == 218) {ii = io + 2; jj = jo - 2; kk = ko + 3; return;} //
	if (index == 219) {ii = io + 2; jj = jo - 1; kk = ko + 3; return;} //

	if (index == 220) {ii = io - 2; jj = jo    ; kk = ko + 3; return;} // Up
	if (index == 221) {ii = io - 2; jj = jo + 1; kk = ko + 3; return;} //
	if (index == 222) {ii = io - 2; jj = jo + 2; kk = ko + 3; return;} //
	if (index == 223) {ii = io - 2; jj = jo + 3; kk = ko + 3; return;} // Left
	if (index == 224) {ii = io - 2; jj = jo + 3; kk = ko + 2; return;} // 
	if (index == 225) {ii = io - 2; jj = jo + 3; kk = ko + 1; return;} //
	if (index == 226) {ii = io - 2; jj = jo + 3; kk = ko    ; return;} //
	if (index == 227) {ii = io - 2; jj = jo + 3; kk = ko - 1; return;} //
	if (index == 228) {ii = io - 2; jj = jo + 3; kk = ko - 2; return;} //
	if (index == 229) {ii = io - 2; jj = jo + 3; kk = ko - 3; return;} // down
	if (index == 230) {ii = io - 2; jj = jo + 2; kk = ko - 3; return;} //
	if (index == 231) {ii = io - 2; jj = jo + 1; kk = ko - 3; return;} //
	if (index == 232) {ii = io - 2; jj = jo    ; kk = ko - 3; return;} //
	if (index == 233) {ii = io - 2; jj = jo - 1; kk = ko - 3; return;} //
	if (index == 234) {ii = io - 2; jj = jo - 2; kk = ko - 3; return;} //
	if (index == 235) {ii = io - 2; jj = jo - 3; kk = ko - 3; return;} // Right
	if (index == 236) {ii = io - 2; jj = jo - 3; kk = ko - 2; return;} //
	if (index == 237) {ii = io - 2; jj = jo - 3; kk = ko - 1; return;} //
	if (index == 238) {ii = io - 2; jj = jo - 3; kk = ko    ; return;} //
	if (index == 239) {ii = io - 2; jj = jo - 3; kk = ko + 1; return;} //
	if (index == 240) {ii = io - 2; jj = jo - 3; kk = ko + 2; return;} //
	if (index == 241) {ii = io - 2; jj = jo - 3; kk = ko + 3; return;} // Up
	if (index == 242) {ii = io - 2; jj = jo - 2; kk = ko + 3; return;} //
	if (index == 243) {ii = io - 2; jj = jo - 1; kk = ko + 3; return;} //


	if (index == 244) {ii = io + 3; jj = jo    ; kk = ko    ; return;} // right
	if (index == 245) {ii = io + 3; jj = jo    ; kk = ko + 1; return;} // up
	if (index == 246) {ii = io + 3; jj = jo + 1; kk = ko + 1; return;} // Left
	if (index == 247) {ii = io + 3; jj = jo + 1; kk = ko    ; return;}
	if (index == 248) {ii = io + 3; jj = jo + 1; kk = ko - 1; return;} // down
	if (index == 249) {ii = io + 3; jj = jo    ; kk = ko - 1; return;}
	if (index == 250) {ii = io + 3; jj = jo - 1; kk = ko - 1; return;} // right
	if (index == 251) {ii = io + 3; jj = jo - 1; kk = ko    ; return;}
	if (index == 252) {ii = io + 3; jj = jo - 1; kk = ko + 1; return;}

	if (index == 253) {ii = io - 3; jj = jo    ; kk = ko    ; return;} // right
	if (index == 254) {ii = io - 3; jj = jo    ; kk = ko + 1; return;} // up
	if (index == 255) {ii = io - 3; jj = jo + 1; kk = ko + 1; return;} // Left
	if (index == 256) {ii = io - 3; jj = jo + 1; kk = ko    ; return;}
	if (index == 257) {ii = io - 3; jj = jo + 1; kk = ko - 1; return;} // down
	if (index == 258) {ii = io - 3; jj = jo    ; kk = ko - 1; return;}
	if (index == 259) {ii = io - 3; jj = jo - 1; kk = ko - 1; return;} // right
	if (index == 260) {ii = io - 3; jj = jo - 1; kk = ko    ; return;}
	if (index == 261) {ii = io - 3; jj = jo - 1; kk = ko + 1; return;}

	if (index == 262) {ii = io + 3; jj = jo    ; kk = ko + 2; return;} // Up
	if (index == 263) {ii = io + 3; jj = jo + 1; kk = ko + 2; return;}
	if (index == 264) {ii = io + 3; jj = jo + 2; kk = ko + 2; return;} // Left
	if (index == 265) {ii = io + 3; jj = jo + 2; kk = ko + 1; return;}
	if (index == 266) {ii = io + 3; jj = jo + 2; kk = ko    ; return;}
	if (index == 267) {ii = io + 3; jj = jo + 2; kk = ko - 1; return;}
	if (index == 268) {ii = io + 3; jj = jo + 2; kk = ko - 2; return;} // down
	if (index == 269) {ii = io + 3; jj = jo + 1; kk = ko - 2; return;}
	if (index == 270) {ii = io + 3; jj = jo    ; kk = ko - 2; return;}
	if (index == 271) {ii = io + 3; jj = jo - 1; kk = ko - 2; return;}
	if (index == 272) {ii = io + 3; jj = jo - 2; kk = ko - 2; return;} // right
	if (index == 273) {ii = io + 3; jj = jo - 2; kk = ko - 1; return;}
	if (index == 274) {ii = io + 3; jj = jo - 2; kk = ko    ; return;}
	if (index == 275) {ii = io + 3; jj = jo - 2; kk = ko + 1; return;}
	if (index == 276) {ii = io + 3; jj = jo - 2; kk = ko + 2; return;} // up
	if (index == 277) {ii = io + 3; jj = jo - 1; kk = ko + 2; return;}

	if (index == 278) {ii = io - 3; jj = jo    ; kk = ko + 2; return;} // Up
	if (index == 279) {ii = io - 3; jj = jo + 1; kk = ko + 2; return;}
	if (index == 280) {ii = io - 3; jj = jo + 2; kk = ko + 2; return;} // Left
	if (index == 281) {ii = io - 3; jj = jo + 2; kk = ko + 1; return;}
	if (index == 282) {ii = io - 3; jj = jo + 2; kk = ko    ; return;}
	if (index == 283) {ii = io - 3; jj = jo + 2; kk = ko - 1; return;}
	if (index == 284) {ii = io - 3; jj = jo + 2; kk = ko - 2; return;} // down
	if (index == 285) {ii = io - 3; jj = jo + 1; kk = ko - 2; return;}
	if (index == 286) {ii = io - 3; jj = jo    ; kk = ko - 2; return;}
	if (index == 287) {ii = io - 3; jj = jo - 1; kk = ko - 2; return;}
	if (index == 288) {ii = io - 3; jj = jo - 2; kk = ko - 2; return;} // right
	if (index == 289) {ii = io - 3; jj = jo - 2; kk = ko - 1; return;}
	if (index == 290) {ii = io - 3; jj = jo - 2; kk = ko    ; return;}
	if (index == 291) {ii = io - 3; jj = jo - 2; kk = ko + 1; return;}
	if (index == 292) {ii = io - 3; jj = jo - 2; kk = ko + 2; return;} // up
	if (index == 293) {ii = io - 3; jj = jo - 1; kk = ko + 2; return;}

	if (index == 294) {ii = io + 3; jj = jo    ; kk = ko + 3; return;} // Up
	if (index == 295) {ii = io + 3; jj = jo + 1; kk = ko + 3; return;} //
	if (index == 296) {ii = io + 3; jj = jo + 2; kk = ko + 3; return;} //
	if (index == 297) {ii = io + 3; jj = jo + 3; kk = ko + 3; return;} // Left
	if (index == 298) {ii = io + 3; jj = jo + 3; kk = ko + 2; return;} // 
	if (index == 299) {ii = io + 3; jj = jo + 3; kk = ko + 1; return;} //
	if (index == 300) {ii = io + 3; jj = jo + 3; kk = ko    ; return;} //
	if (index == 301) {ii = io + 3; jj = jo + 3; kk = ko - 1; return;} //
	if (index == 302) {ii = io + 3; jj = jo + 3; kk = ko - 2; return;} //
	if (index == 303) {ii = io + 3; jj = jo + 3; kk = ko - 3; return;} // down
	if (index == 304) {ii = io + 3; jj = jo + 2; kk = ko - 3; return;} //
	if (index == 305) {ii = io + 3; jj = jo + 1; kk = ko - 3; return;} //
	if (index == 306) {ii = io + 3; jj = jo    ; kk = ko - 3; return;} //
	if (index == 307) {ii = io + 3; jj = jo - 1; kk = ko - 3; return;} //
	if (index == 308) {ii = io + 3; jj = jo - 2; kk = ko - 3; return;} //
	if (index == 309) {ii = io + 3; jj = jo - 3; kk = ko - 3; return;} // Right
	if (index == 310) {ii = io + 3; jj = jo - 3; kk = ko - 2; return;} //
	if (index == 311) {ii = io + 3; jj = jo - 3; kk = ko - 1; return;} //
	if (index == 312) {ii = io + 3; jj = jo - 3; kk = ko    ; return;} //
	if (index == 313) {ii = io + 3; jj = jo - 3; kk = ko + 1; return;} //
	if (index == 314) {ii = io + 3; jj = jo - 3; kk = ko + 2; return;} //
	if (index == 315) {ii = io + 3; jj = jo - 3; kk = ko + 3; return;} // Up
	if (index == 316) {ii = io + 3; jj = jo - 2; kk = ko + 3; return;} //
	if (index == 317) {ii = io + 3; jj = jo - 1; kk = ko + 3; return;} //

	if (index == 318) {ii = io - 3; jj = jo    ; kk = ko + 3; return;} // Up
	if (index == 319) {ii = io - 3; jj = jo + 1; kk = ko + 3; return;} //
	if (index == 320) {ii = io - 3; jj = jo + 2; kk = ko + 3; return;} //
	if (index == 321) {ii = io - 3; jj = jo + 3; kk = ko + 3; return;} // Left
	if (index == 322) {ii = io - 3; jj = jo + 3; kk = ko + 2; return;} // 
	if (index == 323) {ii = io - 3; jj = jo + 3; kk = ko + 1; return;} //
	if (index == 324) {ii = io - 3; jj = jo + 3; kk = ko    ; return;} //
	if (index == 325) {ii = io - 3; jj = jo + 3; kk = ko - 1; return;} //
	if (index == 326) {ii = io - 3; jj = jo + 3; kk = ko - 2; return;} //
	if (index == 327) {ii = io - 3; jj = jo + 3; kk = ko - 3; return;} // down
	if (index == 328) {ii = io - 3; jj = jo + 2; kk = ko - 3; return;} //
	if (index == 329) {ii = io - 3; jj = jo + 1; kk = ko - 3; return;} //
	if (index == 330) {ii = io - 3; jj = jo    ; kk = ko - 3; return;} //
	if (index == 331) {ii = io - 3; jj = jo - 1; kk = ko - 3; return;} //
	if (index == 332) {ii = io - 3; jj = jo - 2; kk = ko - 3; return;} //
	if (index == 333) {ii = io - 3; jj = jo - 3; kk = ko - 3; return;} // Right
	if (index == 334) {ii = io - 3; jj = jo - 3; kk = ko - 2; return;} //
	if (index == 335) {ii = io - 3; jj = jo - 3; kk = ko - 1; return;} //
	if (index == 336) {ii = io - 3; jj = jo - 3; kk = ko    ; return;} //
	if (index == 337) {ii = io - 3; jj = jo - 3; kk = ko + 1; return;} //
	if (index == 338) {ii = io - 3; jj = jo - 3; kk = ko + 2; return;} //
	if (index == 339) {ii = io - 3; jj = jo - 3; kk = ko + 3; return;} // Up
	if (index == 340) {ii = io - 3; jj = jo - 2; kk = ko + 3; return;} //
	if (index == 341) {ii = io - 3; jj = jo - 1; kk = ko + 3; return;} //

	#pragma endregion
}

inline double MeshingGenie_3d::generate_a_random_number()
{   
	#pragma region Generate a Random Number:
    size_t num(rand());
    num = num << 15;
    num += rand(); // another random number
    return static_cast<double>(num * 1.0 / _MY_RAND_MAX);
	#pragma endregion
}

inline void MeshingGenie_3d::get_active_children(size_t i, size_t j, size_t k, size_t icell,
	                            bool &c1, bool &c2, bool &c3, bool &c4, bool &c5, bool &c6, bool &c7, bool &c8)
{
	#pragma region Retrieve Active Children:

	double xx, yy, zz;
	xx = _xo + i * _ss;
	yy = _yo + j * _ss;
	zz = _zo + k * _ss;

	Point* pp;

	double xmin(xx), ymin(yy), zmin(zz);
	c1 = false; c2 = false; c3 = false; c4 = false;
	c5 = false; c6 = false; c7 = false; c8 = false;
	double ss = 0.5 * _ss;

	size_t jcell;
	for (size_t index = 0; index < 124; index++) 
	{
		jcell = icell + _neighbors[index];
		pp = _cell_points[jcell];
		if (pp == 0) continue;
		
		for (size_t child = 0; child < 8; child++)
		{
			if (child == 0 && c1) continue;
			if (child == 1 && c2) continue;
			if (child == 2 && c3) continue;
			if (child == 3 && c4) continue;
			if (child == 4 && c5) continue;
			if (child == 5 && c6) continue;
			if (child == 6 && c7) continue;
			if (child == 7 && c8) continue;

			if (child == 0) {xmin = xx; ymin = yy; zmin = zz;}
			else if (child == 1) {xmin = xx; ymin = yy;      zmin = zz + ss;}
			else if (child == 2) {xmin = xx; ymin = yy + ss; zmin = zz;}
			else if (child == 3) {xmin = xx; ymin = yy + ss; zmin = zz + ss;}
			else if (child == 4) {xmin = xx + ss; ymin = yy; zmin = zz;}
			else if (child == 5) {xmin = xx + ss; ymin = yy;      zmin = zz + ss;}
			else if (child == 6) {xmin = xx + ss; ymin = yy + ss; zmin = zz;}
			else if (child == 7) {xmin = xx + ss; ymin = yy + ss; zmin = zz + ss;}

			double dx(xmin - pp->x), dxx(fabs(dx + ss)); dx = fabs(dx);
			if (dxx > dx) dx = dxx;

			double dy(ymin - pp->y), dyy(fabs(dy + ss)); dy = fabs(dy);
			if (dyy > dy) dy = dyy;

			double dz(zmin - pp->z), dzz(fabs(dz + ss)); dz = fabs(dz);
			if (dzz > dz) dz = dzz;

			if (dx * dx + dy * dy + dz * dz < _dm_squared)
			{
				if (child == 0) c1 = true;
				else if (child == 1) c2 = true;
				else if (child == 2) c3 = true;
				else if (child == 3) c4 = true;
				else if (child == 4) c5 = true;
				else if (child == 5) c6 = true;
				else if (child == 6) c7 = true;
				else if (child == 7) c8 = true;

				if (c1 && c2 && c3 && c4 && c5 && c6 && c7 && c8) return;
			}
		}
	}
	#pragma endregion
}


void MeshingGenie_3d::test_output()
{
	#pragma region Test Output:
	double dmin = 100 * _dm_squared;
	Point* p; Point* q; size_t icell, jcell;
	for (size_t i = 0; i < _ni; i++)
	{
		for (size_t j = 0; j < _nj; j++)
		{
			for (size_t k = 0; k < _nk; k++)
			{
				get_cell_index(i, j, k, icell);
				p = _cell_points[icell];
				if (p == 0) continue;
				// inner loop
				for (size_t ii = 0; ii < _ni; ii++)
				{
					for (size_t jj = 0; jj < _nj; jj++)
					{
						for (size_t kk = 0; kk < _nk; kk++)
						{
							get_cell_index(ii, jj, kk, jcell);
							q = _cell_points[jcell];
							if (q == 0) continue;

							if (p == q) continue;

							// check the empty sphere principal
							double dstsq = distance_squared(p->x - q->x, p->y - q->y, p->z - q->z);
							if (dstsq < dmin) dmin = dstsq;
						}
					}		
				}
			}
		}		
	}
	std::cout<<"Minimum distance between any two points is " << sqrt(dmin) << std::endl;

	// checking that all cells are covered
	_ss = _s / pow(2.0, 14);
	size_t ni(_ni - 2), nj(_nj - 2), nk(_nk - 2);
	for (size_t i = 2; i < ni; i++)
	{
		for (size_t j = 2; j < nj; j++)
		{
			for (size_t k = 2; k < nk; k++)
			{
				// check if these grid points are covered

				for (size_t nit = 0; nit < 100; nit++)
				{
					double u = generate_a_random_number();
					double v = generate_a_random_number();
					double w = generate_a_random_number();
					double xx = _xo + (i + u) * _s;
					double yy = _yo + (j + v) * _s;
					double zz = _zo + (k + w) * _s;

					// inner loop
					double dmin = 100 * _dm_squared;
					for (size_t ii = 0; ii < _ni; ii++)
					{
						for (size_t jj = 0; jj < _nj; jj++)
						{
							for (size_t kk = 0; kk < _nk; kk++)
							{
								get_cell_index(ii, jj, kk, jcell);
								q = _cell_points[jcell];
								if (q == 0) continue;

								// check the empty sphere proncipal
								double dstsq = distance_squared(xx - q->x, yy - q->y, zz - q->z);
								if (dstsq < dmin) 
								{
									dmin = dstsq;
								}
							}
						}
					}
				}
				if (dmin > _dm_squared + 1E-7)
				{
					std::cout<< "Cell (" << i << " , " << j << " , " << k << " ) is not covered!!" << std::endl;					
				}				
			}
		}
	}
	#pragma endregion
}

void MeshingGenie_3d::generate_CVM()
{
	#pragma region Generating Voronoi Cells:
	std::cout << "\n*** Generating Voronoi Cells ...";

	clock_t start_time, end_time; double cpu_time;	 
	start_time = clock();

	_v_cells = new Polyhedron*[_num_inserted_points];

	_cell_vpoints = new VoronoiPoint*[_n];
	for (size_t icell = 0; icell< _n; icell++) _cell_vpoints[icell] = 0;

	delete [] _neighbors;

	size_t jcell, jcell_o;
	get_cell_index(5, 5, 5, jcell_o);

	std::vector<int> neighbors;
	for (size_t i = 0; i <= 10; i++)
	{
		for (size_t j = 0; j <= 10; j++)
		{
			for (size_t k = 0; k <= 10; k++)
			{
				get_cell_index(i, j, k, jcell);
				neighbors.push_back(jcell - jcell_o);
			}
		}
	}

	size_t num_neighbors(0); int num_layers(5);
	if (true)
	{
		#pragma region order layers of neighbors:
		// order neighbors based on distance from jcell center (5 layers)
		double xo(_xo + (num_layers + 0.5) * _s), yo(_yo + (num_layers + 0.5) * _s), zo(_zo + (num_layers + 0.5) * _s);
		std::vector<double> dd_fc; // furthest corner
		std::vector<double> dd_cc; // closest corner
		
		for (int i = 0; i <= 2 * num_layers; i++)
		{
			double xmin(_xo + i * _s), xmax(_xo + (i + 1) * _s);
			
			double xf(xmin);
			if (fabs(xmax - xo) > fabs(xmin - xo)) xf = xmax;
			double dxf(xf - xo), dxc((abs(int(i) - num_layers) - 1) * _s);
			
			for (int j = 0; j <= 2 * num_layers; j++)
			{
				double ymin(_yo + j * _s), ymax(_yo + (j + 1) * _s);
				
				double yf(ymin);
				if (fabs(ymax - yo) > fabs(ymin - yo)) yf = ymax;
				double dyf(yf - yo), dyc((abs(int(j) - num_layers) - 1) * _s);
				
				for (int k = 0; k <= 2 * num_layers; k++)
				{
					double zmin(_zo + k * _s), zmax(_zo + (k + 1) * _s);

					double zf(zmin);
					if (fabs(zmax - zo) > fabs(zmin - zo)) zf = zmax;
					double dzf(zf - zo), dzc((abs(int(k) - num_layers) - 1) * _s);

					dd_fc.push_back(dxf * dxf + dyf * dyf + dzf * dzf);
					dd_cc.push_back(dxc * dxc + dyc * dyc + dzc * dzc);
				}
			}
		}
		size_t num(neighbors.size());
		for (size_t i = 0; i < num; i++)
		{
			for (size_t j = i + 1; j < num; j++)
			{
				if (dd_fc[j] < dd_fc[i])
				{
					int di = neighbors[i]; neighbors[i] = neighbors[j]; neighbors[j] = di;
					double df = dd_fc[i]; dd_fc[i] = dd_fc[j]; dd_fc[j] = df;
					double dc = dd_cc[i]; dd_cc[i] = dd_cc[j]; dd_cc[j] = dc;
				}			
			}
		}
		num_neighbors = num;
		_neighbors = new int[num_neighbors];
		_num_neighbors = num;;
		for (size_t i = 0; i < num; i++)
		{
			_neighbors[i] = neighbors[i];
		}
		#pragma endregion
	}	
	
	PolyFace* f;
	VoronoiPoint* vpnt;
	FacePoint* pnt;
	FacePoint* fpnt;

	size_t iv(0); bool failed(false); double ratio(0.1);

	// working on internal boundary nodes first
	for (size_t icell = 0; icell < _n; icell++)
	{	
		if (_cell_points[icell] == 0) continue;

		construct_internal_boundary_voronoi_cell(icell, iv);
		
		if (iv > ratio * _num_inserted_points)
		{
			std::cout << "\n " << ratio * 100 << "% ...";
			ratio += 0.1;
		}
	}

	// working on reflex edges second
	for (size_t icell = 0; icell < _n; icell++)
	{	
		if (_cell_points[icell] == 0) continue;

		construct_reflex_voronoi_cell(icell, iv);
		
		if (iv > ratio * _num_inserted_points)
		{
			std::cout << "\n " << ratio * 100 << "% ...";
			ratio += 0.1;
		}
	}

	// work on regular voronoi cells finally
	for (size_t icell = 0; icell < _n; icell++)
	{	
		if (_cell_points[icell] == 0) continue;
		
		construct_regular_voronoi_cell(icell, iv);

		if (iv > ratio * _num_inserted_points)
		{
			std::cout << "\n " << ratio * 100 << "% ...";
			ratio += 0.1;
		}
	}

	end_time = clock();
	cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

	 
	std::cout << " done! Execution Time = " << cpu_time << " seconds." << std::endl;

	std::cout<< "\n*** Extracting Output ...";

	start_time = clock();

	// output:
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	std::vector<std::vector<size_t> > faces;
	std::vector<double> faces_normals;
	std::vector<size_t> elements;
	std::vector<std::vector<size_t> > bfaces;

	std::map<VoronoiPoint* , size_t> vpoints_map;
	
	size_t ipoint(0);
	for (size_t icell = 0; icell < _n; icell++) // loop over grid cells to collect Voronoi points
	{
		vpnt =  _cell_vpoints[icell];
		while (vpnt != 0)
		{
			vpoints_map[vpnt] = ipoint; ipoint++; 
			vpnt = vpnt->nextVpoint;
		}
	}
	x.resize(ipoint); y.resize(ipoint); z.resize(ipoint);

	ipoint = 0;
	for (size_t icell = 0; icell < _n; icell++) // loop over grid cells to collect Voronoi points
	{
		vpnt =  _cell_vpoints[icell];
		while (vpnt != 0)
		{
			x[ipoint] = vpnt->x; y[ipoint] = vpnt->y; z[ipoint] = vpnt->z; ipoint++; 
			vpnt = vpnt->nextVpoint;
		}
	}

	size_t num_vcells(_num_inserted_points);

	// count number of faces
	size_t num_faces(0), num_b_faces(0);
	for (size_t icell = 0; icell < num_vcells; icell++) // loop over voronoi cells
	{
		// loop over element faces
		Polyhedron* poly = _v_cells[icell];
		while (poly != 0)
		{
			f = poly->face; 
			while (f != 0)
			{			
				if (f->on_boundary) num_b_faces++;
				f = f->nextface; num_faces++;
			}
			poly = poly->nextpoly;
		}
	}

	size_t iface(0);
	faces.resize(num_faces); bfaces.resize(num_b_faces);
	faces_normals.reserve(num_faces * 3);
	elements.reserve(num_vcells + 1);
	for (size_t icell = 0; icell < num_vcells; icell++) // loop over voronoi cells
	{
		// loop over element faces
		Polyhedron* poly = _v_cells[icell];
		while (poly != 0)
		{
			elements.push_back(iface);
			f = poly->face; size_t jface(0);
			while (f != 0)
			{
				pnt  = f->corner;
				fpnt = pnt;
				// loop over face corners
				while (true)
				{
					if (f->on_boundary) bfaces[jface].push_back(vpoints_map[pnt->Vpoint]);
					faces[iface].push_back(vpoints_map[pnt->Vpoint]);
					pnt = pnt->nextcorner;
					if (pnt == fpnt)break;
				}
				faces_normals.push_back(f->nx);
				faces_normals.push_back(f->ny);
				faces_normals.push_back(f->nz);
				if (f->on_boundary) jface++;
				f = f->nextface; iface++;
			}
			poly = poly->nextpoly;
		}
	}
	elements.push_back(iface);

	//save_face_normals(faces_normals);

	end_time = clock();
	cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

	std::cout<< " done! Execution Time = " << cpu_time << " seconds." << std::endl;

	std::cout << "\n*** Number of Voronoi Cells  = " << num_vcells << std::endl;

	save_vcells(x, y, z, faces, elements);

	save_vtk_polygons(x, y, z, faces);

	if (true)
	{
		// Viewer clipping
		/*
		std::vector< std::vector<size_t> > clippedfaces;
		double xc(0.5 * (_box_xmin + _box_xmax));
		double yc(0.5 * (_box_ymin + _box_ymax));
		double zc(0.5 * (_box_zmin + _box_zmax));
		size_t num_elements(elements.size() - 1);
		for (size_t ie = 0; ie < num_elements; ie++)
		{
			bool valid_element(true);
			for (size_t i = elements[ie]; i < elements[ie + 1]; i++)
			{
				size_t num_corners(faces[i].size());
				for (size_t j = 0; j < num_corners; j++)
				{
					double xx = x[faces[i][j]];
					double yy = y[faces[i][j]];
					double zz = z[faces[i][j]];
					if (xx < xc && yy < yc && zz < zc)
					{
						valid_element = false; break;
					}
				}
				if (!valid_element) break;
			}
			if (valid_element)
			{
				for (size_t i = elements[ie]; i < elements[ie + 1]; i++)			
				{				
					clippedfaces.push_back(faces[i]);				
				}
			}
		}
		save_vtk_polygons(x, y, z, clippedfaces);
		*/
	}

	//save_ply_mesh("final_mesh.ply", x, y, z, faces);

	//test_voronoi_cells(x, y, z, faces, faces_normals);

	return;
	#pragma endregion
}

inline bool MeshingGenie_3d::construct_internal_boundary_voronoi_cell(size_t icell, size_t &v_cell_index)
{
	#pragma region Construct a Voronoi Cell around an internal boundary point:

	size_t i, j, k;
	get_cell_indices(icell, i, j, k);

	// cell center:
	double xc(_xo + (i + 0.5) * _s);
	double yc(_yo + (j + 0.5) * _s);
	double zc(_zo + (k + 0.5) * _s);

	size_t jcell;

	Point* p; Point* q;
	double xmid, ymid, zmid;
	double nx, ny, nz;
	PolyFace* f;
	PolyFace* ff;
	PolyFace* fn;
	PolyFace* face;
	FacePoint* pnt;
	Polyhedron* poly;

	p = _cell_points[icell];

	while (p != 0)
	{
		if (!p->_on_internal_boundaries)
		{
			p = p->next;
			continue;
		}

		//std::cout<< "icell = " << icell << std::endl;

		poly = new Polyhedron();
		poly->face = 0; poly->nextpoly =  0;

		std::vector<size_t> faces_v;
		std::set<size_t> faces;
		std::vector<size_t> int_faces_v;
		std::set<size_t> int_faces;

		std::set<size_t>* cell_faces;
		std::set<size_t>::iterator faces_iter;
		std::set<size_t>::iterator faces_iter_j;

		for (size_t index = 0; index < _num_neighbors; index++)
		{
			#pragma region Retrieve Neighbor boundary facess:			
			jcell = icell + _neighbors[index];

			if (!_boundary_cells[jcell]) continue;

			cell_faces = (&_cell_faces[jcell]);

			for (faces_iter = cell_faces->begin(); faces_iter != cell_faces->end(); faces_iter++)
			{
				size_t iface(*faces_iter);
				if (faces.find(iface) == faces.end())
				{
					faces.insert(iface); 
					faces_v.push_back(iface);
				}			
			}
			#pragma endregion
		}

		for (size_t index = 0; index < _num_neighbors; index++)
		{
			#pragma region Retrieve Neighbor boundary facess:			
			jcell = icell + _neighbors[index];

			if (_cell_internal_faces.find(jcell) == _cell_internal_faces.end()) continue;

			cell_faces = (&_cell_internal_faces[jcell]);

			for (faces_iter = cell_faces->begin(); faces_iter != cell_faces->end(); faces_iter++)
			{
				size_t iface(*faces_iter);
				if (int_faces.find(iface) == int_faces.end())
				{
					int_faces.insert(iface); 
					int_faces_v.push_back(iface);
				}
			}
			#pragma endregion
		}
	
		size_t num_i_faces(int_faces.size());
		std::vector<PolyFace*> internal_faces(num_i_faces);
		for (size_t i = 0; i < num_i_faces; i++)
		{
			#pragma region Creating boundary faces:
			
			size_t iface(int_faces_v[i]);
			
			bool duplicated(false);

			for (size_t j = 0; j < i; j++)
			{
				face = internal_faces[j];
				if (face == 0) continue;

				double dot(_nfx[iface] * face->nx + _nfy[iface] * face->ny + _nfz[iface] * face->nz);
				if (dot <= 0.9) continue;

				double ax = _xf[iface][0] - face->corner->Vpoint->x;
				double ay = _yf[iface][0] - face->corner->Vpoint->y;
				double az = _zf[iface][0] - face->corner->Vpoint->z;
				double a = sqrt(ax * ax + ay * ay + az * az);
				if (fabs(ax * face->nx + ay * face->ny + az * face->nz) > 1E-8 * a) continue; 
				
				duplicated = true;
				break;				
			}

			if (duplicated)
			{
				internal_faces[i] = 0;
			}
			else
			{
				internal_faces[i] = construct_super_triangle(_xf[iface][0], _yf[iface][0], _zf[iface][0], 
															 _nfx[iface], _nfy[iface], _nfz[iface],
															 xc, yc, zc, 0);

				if (internal_faces[i]->num_corners < 3)
				{
					delete_face_corners(internal_faces[i]);
					delete internal_faces[i];
					internal_faces[i] = 0;
				}
			}
			#pragma endregion
		}

		for (size_t index = 0; index < _num_neighbors; index++)
		{
			#pragma region Adding Neighbor points:			
			jcell = icell + _neighbors[index];

			q = _cell_points[jcell];
		
			if (q == 0) continue;

			while (q != 0)
			{
				if (q == p)
				{
					q = q->next;
					continue;
				}

				Polyhedron* q_poly(q->poly);
				while (q_poly != 0)
				{
					#pragma region Retrieve common faces if any:
					f = poly->face; fn = 0;
					while (f != 0)  			
					{
						fn = f;
						f = f->nextface;
					}

					bool first(true);
					f = q_poly->face;
					while (f != 0)
					{
						if (f->q == p)	
						{
							pnt = f->corner; size_t num_face_nodes(f->num_corners);
							FacePoint** fpnts = new FacePoint*[num_face_nodes];
							for (size_t inode = 0; inode < num_face_nodes; inode++)
							{
								fpnts[inode] = new FacePoint();
								fpnts[inode]->Vpoint = new VoronoiPoint();
								fpnts[inode]->Vpoint->on_boundary = pnt->Vpoint->on_boundary;
								fpnts[inode]->Vpoint->x = pnt->Vpoint->x;
								fpnts[inode]->Vpoint->y = pnt->Vpoint->y;
								fpnts[inode]->Vpoint->z = pnt->Vpoint->z;
								fpnts[inode]->Vpoint->polys = pnt->Vpoint->polys;
								fpnts[inode]->Vpoint->nextVpoint = pnt->Vpoint->nextVpoint;
								pnt = pnt->prevcorner;
							}

							for (size_t inode = 0; inode < num_face_nodes; inode++)
							{
								size_t inodep = inode + 1;
								if (inodep == num_face_nodes) inodep = 0;
								fpnts[inode]->nextcorner = fpnts[inodep];
								fpnts[inodep]->prevcorner = fpnts[inode];
							}

							face = copy_polyface(f, true);
			
							if (first)
							{
							
								// trim old faces using cf
								ff = poly->face;
								while (ff != 0)  
								{
									if (ff->num_corners > 0)
									{
										double xx = face->corner->Vpoint->x;
										double yy = face->corner->Vpoint->y;
										double zz = face->corner->Vpoint->z;
										double px, py, pz, qx, qy, qz;
										trim_face(ff, xx, yy, zz, face->nx, face->ny, face->nz, px, py, pz, qx, qy, qz);
									}
									ff = ff->nextface;
								}
							
								first = false;
							}

							// add face to p poly
							if (fn != 0) fn->nextface = face;
							else poly->face = face;

							fn = face;

							delete[] fpnts;
						}

						f = f->nextface; // next face in q polyhedron
					}
					q_poly = q_poly->nextpoly;
					#pragma endregion
				}

				if (q->poly != 0)
				{
					q = q->next;
					continue;				
				}

				// construct a super triangle representing the mid plane
				xmid = 0.5 * (q->x + p->x); ymid = 0.5 * (q->y + p->y); zmid = 0.5 * (q->z + p->z);
				nx = (q->x - p->x); ny = (q->y - p->y); nz = (q->z - p->z);
				
				face = construct_super_triangle(xmid, ymid, zmid, nx, ny, nz, xc, yc, zc, q);

				if (face->num_corners < 3)
				{
					q = q->next;
					continue;				
				}

				// trim new face using previous internal faces:
				f = poly->face;
				while (f != 0)  
				{
					if (face->num_corners > 0 && f->num_corners > 0)
					{
						double xx = f->corner->Vpoint->x;
						double yy = f->corner->Vpoint->y;
						double zz = f->corner->Vpoint->z;
						double px, py, pz, qx, qy, qz;
						trim_face(face, xx, yy, zz, f->nx, f->ny, f->nz, px, py, pz, qx, qy, qz);
					}
					if (face->num_corners < 3) break;
					f = f->nextface;
				}

				if (face->num_corners < 3) 
				{
					q  = q->next;
					continue;
				}
								
				// trim old faces using face
				f = poly->face;
				while (f != 0)  
				{
					if (!f->preserved && f->num_corners > 0)
					{
						double xx = face->corner->Vpoint->x;
						double yy = face->corner->Vpoint->y;
						double zz = face->corner->Vpoint->z;
						double px, py, pz, qx, qy, qz;
						trim_face(f, xx, yy, zz, face->nx, face->ny, face->nz, px, py, pz, qx, qy, qz);
					}
					f = f->nextface;
				}				

				f = poly->face; fn = 0;
				while (f != 0)  			
				{
					fn = f;
					f = f->nextface;
				}
			
				if (fn == 0) poly->face = face;
				else fn->nextface = face;

				q  = q->next;
			}
			#pragma endregion
		}

		size_t num_b_faces(faces.size());
		std::vector<PolyFace*> boundary_faces(num_b_faces);
		for (size_t i = 0; i < num_b_faces; i++)
		{
			#pragma region Creating boundary faces:
			
			// create a face and trim it using existing faces
			size_t iface(faces_v[i]);

			face = create_boundary_face(iface);			
			PolyFace* f = poly->face;
			while (f != 0)
			{
				// trim face with internal faces only listed first
				double rx, ry, rz, sx, sy, sz;
				if (f->num_corners >= 3) 
				{
					trim_face(face, f->corner->Vpoint->x, f->corner->Vpoint->y, f->corner->Vpoint->z, f->nx, f->ny, f->nz,
								rx, ry, rz, sx, sy, sz);
					
					if (face->num_corners < 3) break;
				}
				f = f->nextface;
			}

			if (face->num_corners < 3) 
			{
				delete_face_corners(face);
				delete face;
				continue;
			}


			size_t ii = _bfaces[iface][0];
			bool duplicated(false);

			for (size_t j = 0; j < i; j++)
			{
				face = boundary_faces[j];
				if (face == 0) continue;
				double dot(_bfaces_nx[iface] * face->nx + _bfaces_ny[iface] * face->ny + _bfaces_nz[iface] * face->nz);
				if (dot > 0.9) 
				{
					duplicated = true;
					break;
				}
			}

			if (duplicated)
			{
				boundary_faces[i] = 0;
			}
			else
			{
				boundary_faces[i] = construct_super_triangle(_xb[ii], _yb[ii], _zb[ii], 
															 _bfaces_nx[iface], _bfaces_ny[iface], _bfaces_nz[iface],
															 xc, yc, zc, 0);

				if (boundary_faces[i]->num_corners < 3)
				{
					delete_face_corners(boundary_faces[i]);
					delete boundary_faces[i];
					boundary_faces[i] = 0;
				}
			}
			#pragma endregion
		}		

		for (size_t iface = 0; iface < num_b_faces; iface++)
		{	
			#pragma region Adding boundary faces:
			face = boundary_faces[iface];

			if (face == 0) continue;

			face->on_boundary = true;
			if (!add_face_to_poly(poly, face))
			{
				delete_face_corners(face);
				delete face;
			}
			#pragma endregion
		}
		

		for (size_t iface = 0; iface < num_i_faces; iface++)
		{	
			#pragma region Adding boundary faces:
			face = internal_faces[iface];

			if (face == 0) continue;

			face->on_boundary = true;

			Polyhedron* ppoly(poly);
			Polyhedron* firstnewpoly(0);
			Polyhedron* lastnewpoly(0);
			Polyhedron* lastpoly(0);
			while (ppoly != 0)
			{
				Polyhedron* new_poly = copy_polyhedron(ppoly);				
				if (!add_face_to_poly(new_poly, copy_polyface(face, true)))
				{
					// ppoly is not intersected by face 
					delete_poly_faces(new_poly);
					delete new_poly;
					lastpoly = ppoly;
					ppoly = ppoly->nextpoly;
					continue;
				}
				
				if (firstnewpoly == 0)
				{
					firstnewpoly = new_poly;
					lastnewpoly = new_poly;
				}
				else
				{
					lastnewpoly->nextpoly = new_poly;
					lastnewpoly = new_poly;
				}
				
				add_face_to_poly(ppoly, copy_polyface(face, false));
				lastpoly = ppoly;
				ppoly = ppoly->nextpoly;
			}
			// Delete original copies of face and face inverted
			delete_face_corners(face); delete face;
			lastpoly->nextpoly = firstnewpoly;
			#pragma endregion
		}
		
		// weld corners of the new polys to old Voronoi points
		weld_poly_corners(poly);

		// remove degenerate faces from created polys
		remove_degenerate_faces(poly);

		/*
		if (true)
		{
			Polyhedron* ppoly(poly);
			while (ppoly!=0)
			{
				int bug(0);
				plot_polyhedron(ppoly);
				ppoly = ppoly->nextpoly;
				bug++;
			}
		}
		*/

		_v_cells[v_cell_index] = poly; v_cell_index++;
		p->poly = poly;

		p = p->next;
	}
	return true;
	#pragma endregion
}

inline bool MeshingGenie_3d::construct_reflex_voronoi_cell(size_t icell, size_t &v_cell_index)
{
	#pragma region Construct a Voronoi Cell around a reflex edge point:

	size_t i, j, k;
	get_cell_indices(icell, i, j, k);

	// cell center:
	double xc(_xo + (i + 0.5) * _s);
	double yc(_yo + (j + 0.5) * _s);
	double zc(_zo + (k + 0.5) * _s);

	size_t jcell;

	Point* p; Point* q;
	double xmid, ymid, zmid;
	double nx, ny, nz;
	PolyFace* f;
	PolyFace* ff;
	PolyFace* fn;
	PolyFace* face;
	FacePoint* pnt;
	Polyhedron* poly;

	p = _cell_points[icell];

	while (p != 0)
	{
		if (!p->_on_reflex_corner)
		{
			p = p->next;
			continue;
		}

		//std::cout<< "icell = " << icell << std::endl;

		poly = new Polyhedron();
		poly->face = 0; poly->nextpoly =  0;

		std::vector<size_t> faces_v;
		std::set<size_t> faces;
		std::set<size_t>* cell_faces;
		std::set<size_t>::iterator faces_iter;
		std::set<size_t>::iterator faces_iter_j;

		for (size_t index = 0; index < _num_neighbors; index++)
		{
			#pragma region Retrieve Neighbor boundary facess:			
			jcell = icell + _neighbors[index];

			if (!_boundary_cells[jcell]) continue;

			cell_faces = (&_cell_faces[jcell]);

			for (faces_iter = cell_faces->begin(); faces_iter != cell_faces->end(); faces_iter++)
			{
				size_t iface(*faces_iter);
				if (faces.find(iface) == faces.end())
				{
					faces.insert(iface); 
					faces_v.push_back(iface);
				}			
			}
			#pragma endregion
		}		
	
		size_t num_b_faces(faces.size());
		std::vector<PolyFace*> boundary_faces(num_b_faces);
		for (size_t i = 0; i < num_b_faces; i++)
		{
			#pragma region Creating boundary faces:
			
			size_t iface(faces_v[i]);
			size_t ii = _bfaces[iface][0];
			bool duplicated(false);

			for (size_t j = 0; j < i; j++)
			{
				face = boundary_faces[j];
				if (face == 0) continue;

				double dot(_bfaces_nx[iface] * face->nx + _bfaces_ny[iface] * face->ny + _bfaces_nz[iface] * face->nz);
				if (dot <= 0.9) continue;

				double ax = _xb[ii] - face->corner->Vpoint->x;
				double ay = _yb[ii] - face->corner->Vpoint->y;
				double az = _zb[ii] - face->corner->Vpoint->z;
				double a = sqrt(ax * ax + ay * ay + az * az);
				if (fabs(ax * face->nx + ay * face->ny + az * face->nz) > 1E-8 * a) continue; 
				
				duplicated = true;
				break;				
			}

			if (duplicated)
			{
				boundary_faces[i] = 0;
			}
			else
			{
				boundary_faces[i] = construct_super_triangle(_xb[ii], _yb[ii], _zb[ii], 
															 _bfaces_nx[iface], _bfaces_ny[iface], _bfaces_nz[iface],
															 xc, yc, zc, 0);

				if (boundary_faces[i]->num_corners < 3)
				{
					delete_face_corners(boundary_faces[i]);
					delete boundary_faces[i];
					boundary_faces[i] = 0;
				}
			}
			#pragma endregion
		}

		for (size_t index = 0; index < _num_neighbors; index++)
		{
			#pragma region Adding Neighbor points:			
			jcell = icell + _neighbors[index];

			q = _cell_points[jcell];
		
			if (q == 0) continue;

			while (q != 0)
			{
				if (q == p)
				{
					q = q->next;
					continue;
				}

				Polyhedron* q_poly(q->poly);
				while (q_poly != 0)
				{
					#pragma region Retrieve common faces if any:
					f = poly->face; fn = 0;
					while (f != 0)  			
					{
						fn = f;
						f = f->nextface;
					}

					bool first(true);
					f = q_poly->face;
					while (f != 0)
					{
						if (f->q == p)	
						{
							pnt = f->corner; size_t num_face_nodes(f->num_corners);
							FacePoint** fpnts = new FacePoint*[num_face_nodes];
							for (size_t inode = 0; inode < num_face_nodes; inode++)
							{
								fpnts[inode] = new FacePoint();
								fpnts[inode]->Vpoint = new VoronoiPoint();
								fpnts[inode]->Vpoint->on_boundary = pnt->Vpoint->on_boundary;
								fpnts[inode]->Vpoint->x = pnt->Vpoint->x;
								fpnts[inode]->Vpoint->y = pnt->Vpoint->y;
								fpnts[inode]->Vpoint->z = pnt->Vpoint->z;
								fpnts[inode]->Vpoint->polys = pnt->Vpoint->polys;
								fpnts[inode]->Vpoint->nextVpoint = pnt->Vpoint->nextVpoint;
								pnt = pnt->prevcorner;
							}

							for (size_t inode = 0; inode < num_face_nodes; inode++)
							{
								size_t inodep = inode + 1;
								if (inodep == num_face_nodes) inodep = 0;
								fpnts[inode]->nextcorner = fpnts[inodep];
								fpnts[inodep]->prevcorner = fpnts[inode];
							}

							face = copy_polyface(f, true);
			
							if (first)
							{
							
								// trim old faces using cf
								ff = poly->face;
								while (ff != 0)  
								{
									if (ff->num_corners > 0)
									{
										double xx = face->corner->Vpoint->x;
										double yy = face->corner->Vpoint->y;
										double zz = face->corner->Vpoint->z;
										double px, py, pz, qx, qy, qz;
										trim_face(ff, xx, yy, zz, face->nx, face->ny, face->nz, px, py, pz, qx, qy, qz);
									}
									ff = ff->nextface;
								}
							
								first = false;
							}

							// add face to p poly
							if (fn != 0) fn->nextface = face;
							else poly->face = face;

							fn = face;

							delete[] fpnts;
						}

						f = f->nextface; // next face in q polyhedron
					}
					q_poly = q_poly->nextpoly;
					#pragma endregion
				}

				if (q->poly != 0)
				{
					q = q->next;
					continue;				
				}

				// construct a super triangle representing the mid plane
				xmid = 0.5 * (q->x + p->x); ymid = 0.5 * (q->y + p->y); zmid = 0.5 * (q->z + p->z);
				nx = (q->x - p->x); ny = (q->y - p->y); nz = (q->z - p->z);
				
				face = construct_super_triangle(xmid, ymid, zmid, nx, ny, nz, xc, yc, zc, q);

				if (face->num_corners < 3)
				{
					q = q->next;
					continue;				
				}

				// trim new face using previous internal faces:
				f = poly->face;
				while (f != 0)  
				{
					if (face->num_corners > 0 && f->num_corners > 0)
					{
						double xx = f->corner->Vpoint->x;
						double yy = f->corner->Vpoint->y;
						double zz = f->corner->Vpoint->z;
						double px, py, pz, qx, qy, qz;
						trim_face(face, xx, yy, zz, f->nx, f->ny, f->nz, px, py, pz, qx, qy, qz);
					}
					if (face->num_corners < 3) break;
					f = f->nextface;
				}

				if (face->num_corners < 3) 
				{
					q  = q->next;
					continue;
				}
								
				// trim old faces using face
				f = poly->face;
				while (f != 0)  
				{
					if (!f->preserved && f->num_corners > 0)
					{
						double xx = face->corner->Vpoint->x;
						double yy = face->corner->Vpoint->y;
						double zz = face->corner->Vpoint->z;
						double px, py, pz, qx, qy, qz;
						trim_face(f, xx, yy, zz, face->nx, face->ny, face->nz, px, py, pz, qx, qy, qz);
					}
					f = f->nextface;
				}				

				f = poly->face; fn = 0;
				while (f != 0)  			
				{
					fn = f;
					f = f->nextface;
				}
			
				if (fn == 0) poly->face = face;
				else fn->nextface = face;

				q  = q->next;
			}
			#pragma endregion
		}

		for (size_t iface = 0; iface < num_b_faces; iface++)
		{	
			#pragma region Adding boundary faces:
			face = boundary_faces[iface];

			if (face == 0) continue;

			face->on_boundary = true;

			Polyhedron* ppoly(poly);
			Polyhedron* firstnewpoly(0);
			Polyhedron* lastnewpoly(0);
			Polyhedron* lastpoly(0);
			while (ppoly != 0)
			{
				Polyhedron* new_poly = copy_polyhedron(ppoly);				
				if (!add_face_to_poly(new_poly, copy_polyface(face, true)))
				{
					// ppoly is not intersected by face 
					delete_poly_faces(new_poly);
					delete new_poly;
					lastpoly = ppoly;
					ppoly = ppoly->nextpoly;
					continue;
				}
				
				if (firstnewpoly == 0)
				{
					firstnewpoly = new_poly;
					lastnewpoly = new_poly;
				}
				else
				{
					lastnewpoly->nextpoly = new_poly;
					lastnewpoly = new_poly;
				}
				
				add_face_to_poly(ppoly, copy_polyface(face, false));
				lastpoly = ppoly;
				ppoly = ppoly->nextpoly;
			}
			// Delete original copies of face and face inverted
			delete_face_corners(face); delete face;
			lastpoly->nextpoly = firstnewpoly;
			#pragma endregion
		}
		
		poly = remove_external_polys(poly);

		/*
		if (true)
		{
			Polyhedron* ppoly(poly);
			while (ppoly!=0)
			{
				int bug(0);
				plot_polyhedron(ppoly);
				ppoly = ppoly->nextpoly;
				bug++;
			}
		}
		*/
		
		// weld corners of the new polys to old Voronoi points
		weld_poly_corners(poly);

		// remove degenerate faces from created polys
		remove_degenerate_faces(poly);

		_v_cells[v_cell_index] = poly; v_cell_index++;
		p->poly = poly;

		p = p->next;
	}
	return true;
	#pragma endregion
}

inline bool MeshingGenie_3d::construct_regular_voronoi_cell(size_t icell, size_t &v_cell_index)
{
	#pragma region Construct a Voronoi Cell around a regular edge point:

	size_t i, j, k;
	get_cell_indices(icell, i, j, k);

	// cell center:
	double xc(_xo + (i + 0.5) * _s);
	double yc(_yo + (j + 0.5) * _s);
	double zc(_zo + (k + 0.5) * _s);

	size_t jcell;

	Point* p; Point* q;
	double xmid, ymid, zmid;
	double nx, ny, nz;
	PolyFace* f;
	PolyFace* ff;
	PolyFace* fn;
	PolyFace* face;
	FacePoint* pnt;
	Polyhedron* poly;

	p = _cell_points[icell];

	while (p != 0)
	{
		if (p->poly != 0)
		{
			p = p->next;
			continue;
		}

		poly = new Polyhedron();
		poly->face = 0; poly->nextpoly =  0;

		std::vector<size_t> faces_v;
		std::set<size_t> faces;
		std::set<size_t>* cell_faces;
		std::set<size_t>::iterator faces_iter;
		std::set<size_t>::iterator faces_iter_j;

		for (size_t index = 0; index < _num_neighbors; index++)
		{
			#pragma region Retrieve Neighbor boundary facess:			
			jcell = icell + _neighbors[index];

			if (!_boundary_cells[jcell]) continue;

			cell_faces = (&_cell_faces[jcell]);

			for (faces_iter = cell_faces->begin(); faces_iter != cell_faces->end(); faces_iter++)
			{
				size_t iface(*faces_iter);
				if (faces.find(iface) == faces.end())
				{
					faces.insert(iface); 
					faces_v.push_back(iface);
				}			
			}
			#pragma endregion
		}		

		for (size_t index = 0; index < _num_neighbors; index++)
		{
			#pragma region Adding Neighbor points:			
			jcell = icell + _neighbors[index];

			q = _cell_points[jcell];
		
			if (q == 0) continue;

			while (q != 0)
			{
				if (q == p)
				{
					q = q->next;
					continue;
				}

				Polyhedron* q_poly(q->poly);
				while (q_poly != 0)
				{
					#pragma region Retrieve common faces if any:
					f = poly->face; fn = 0;
					while (f != 0)  			
					{
						fn = f;
						f = f->nextface;
					}

					bool first(true);
					f = q_poly->face;
					while (f != 0)
					{
						if (f->q == p)	
						{
							pnt = f->corner; size_t num_face_nodes(f->num_corners);
							FacePoint** fpnts = new FacePoint*[num_face_nodes];
							for (size_t inode = 0; inode < num_face_nodes; inode++)
							{
								fpnts[inode] = new FacePoint();
								fpnts[inode]->Vpoint = new VoronoiPoint();
								fpnts[inode]->Vpoint->on_boundary = pnt->Vpoint->on_boundary;
								fpnts[inode]->Vpoint->x = pnt->Vpoint->x;
								fpnts[inode]->Vpoint->y = pnt->Vpoint->y;
								fpnts[inode]->Vpoint->z = pnt->Vpoint->z;
								fpnts[inode]->Vpoint->polys = pnt->Vpoint->polys;
								fpnts[inode]->Vpoint->nextVpoint = pnt->Vpoint->nextVpoint;
								pnt = pnt->prevcorner;
							}

							for (size_t inode = 0; inode < num_face_nodes; inode++)
							{
								size_t inodep = inode + 1;
								if (inodep == num_face_nodes) inodep = 0;
								fpnts[inode]->nextcorner = fpnts[inodep];
								fpnts[inodep]->prevcorner = fpnts[inode];
							}

							face = copy_polyface(f, true);
			
							if (first)
							{
							
								// trim old faces using cf
								ff = poly->face;
								while (ff != 0)  
								{
									if (ff->num_corners > 0)
									{
										double xx = face->corner->Vpoint->x;
										double yy = face->corner->Vpoint->y;
										double zz = face->corner->Vpoint->z;
										double px, py, pz, qx, qy, qz;
										trim_face(ff, xx, yy, zz, face->nx, face->ny, face->nz, px, py, pz, qx, qy, qz);
									}
									ff = ff->nextface;
								}
							
								first = false;
							}

							// add face to p poly
							if (fn != 0) fn->nextface = face;
							else poly->face = face;

							fn = face;

							delete[] fpnts;
						}

						f = f->nextface; // next face in q polyhedron
					}
					q_poly = q_poly->nextpoly;
					#pragma endregion
				}

				if (q->poly != 0)
				{
					q = q->next;
					continue;					
				}

				// construct a super triangle representing the mid plane
				xmid = 0.5 * (q->x + p->x); ymid = 0.5 * (q->y + p->y); zmid = 0.5 * (q->z + p->z);
				nx = (q->x - p->x); ny = (q->y - p->y); nz = (q->z - p->z);
				
				face = construct_super_triangle(xmid, ymid, zmid, nx, ny, nz, xc, yc, zc, q);

				if (face->num_corners < 3)
				{
					q = q->next;
					continue;				
				}

				// trim new face using previous internal faces:
				f = poly->face;
				while (f != 0)  
				{
					if (face->num_corners > 0 && f->num_corners > 0)
					{
						double xx = f->corner->Vpoint->x;
						double yy = f->corner->Vpoint->y;
						double zz = f->corner->Vpoint->z;
						double px, py, pz, qx, qy, qz;
						trim_face(face, xx, yy, zz, f->nx, f->ny, f->nz, px, py, pz, qx, qy, qz);
					}
					if (face->num_corners < 3) break;
					f = f->nextface;
				}

				if (face->num_corners < 3) 
				{
					q  = q->next;
					continue;
				}
								
				// trim old faces using face
				f = poly->face;
				while (f != 0)  
				{
					if (!f->preserved && f->num_corners > 0)
					{
						double xx = face->corner->Vpoint->x;
						double yy = face->corner->Vpoint->y;
						double zz = face->corner->Vpoint->z;
						double px, py, pz, qx, qy, qz;
						trim_face(f, xx, yy, zz, face->nx, face->ny, face->nz, px, py, pz, qx, qy, qz);
					}
					f = f->nextface;
				}				

				f = poly->face; fn = 0;
				while (f != 0)  			
				{
					fn = f;
					f = f->nextface;
				}

				if (fn == 0) poly->face = face;
				else fn->nextface = face;

				q  = q->next;
			}
			#pragma endregion
		}

		size_t num_b_faces(faces.size());
		std::vector<PolyFace*> boundary_faces(num_b_faces);
		for (size_t i = 0; i < num_b_faces; i++)
		{
			#pragma region Creating boundary faces:
			
			// create a face and trim it using existing faces
			size_t iface(faces_v[i]);

			face = create_boundary_face(iface);			
			PolyFace* f = poly->face;
			while (f != 0)
			{
				// trim face with internal faces only listed first
				double rx, ry, rz, sx, sy, sz;
				if (f->num_corners >= 3) 
				{
					trim_face(face, f->corner->Vpoint->x, f->corner->Vpoint->y, f->corner->Vpoint->z, f->nx, f->ny, f->nz,
								rx, ry, rz, sx, sy, sz);
					
					if (face->num_corners < 3) break;
				}
				f = f->nextface;
			}

			if (face->num_corners < 3) 
			{
				delete_face_corners(face);
				delete face;
				continue;
			}


			size_t ii = _bfaces[iface][0];
			bool duplicated(false);

			for (size_t j = 0; j < i; j++)
			{
				face = boundary_faces[j];
				if (face == 0) continue;
				double dot(_bfaces_nx[iface] * face->nx + _bfaces_ny[iface] * face->ny + _bfaces_nz[iface] * face->nz);
				if (dot > 0.9) 
				{
					duplicated = true;
					break;
				}
			}

			if (duplicated)
			{
				boundary_faces[i] = 0;
			}
			else
			{
				boundary_faces[i] = construct_super_triangle(_xb[ii], _yb[ii], _zb[ii], 
															 _bfaces_nx[iface], _bfaces_ny[iface], _bfaces_nz[iface],
															 xc, yc, zc, 0);

				if (boundary_faces[i]->num_corners < 3)
				{
					delete_face_corners(boundary_faces[i]);
					delete boundary_faces[i];
					boundary_faces[i] = 0;
				}
			}
			#pragma endregion
		}
		

		for (size_t iface = 0; iface < num_b_faces; iface++)
		{	
			#pragma region Adding boundary faces:
			face = boundary_faces[iface];

			if (face == 0) continue;

			face->on_boundary = true;
			if (!add_face_to_poly(poly, face))
			{
				delete_face_corners(face);
				delete face;
			}
			#pragma endregion
		}
		
		
		// weld corners of the new polys to old Voronoi points
		weld_poly_corners(poly);

		// remove degenerate faces from created polys
		remove_degenerate_faces(poly);

		//plot_polyhedron(poly);

		_v_cells[v_cell_index] = poly; v_cell_index++;
		p->poly = poly;

		p = p->next;
	}
	return true;
	#pragma endregion
}


inline bool MeshingGenie_3d::trim_face(PolyFace* f,
	                  double x1, double y1, double z1,
					  double nx, double ny, double nz,
					  double &px, double &py, double &pz,
					  double &qx, double &qy, double &qz)
{
	#pragma region Trim Face using a plane:
	// mark face corners
	bool first(true);
	FacePoint* pnt; FacePoint* fpnt; FacePoint* tmppnt; VoronoiPoint* vpnt;
	FacePoint* pntp; FacePoint* pntm; VoronoiPoint* vpntp;
	pnt = f->corner; fpnt = pnt;
	size_t num_valid_corner(0);
	while (true)
	{
		// check if point is on the right side of plane
		vpnt = pnt->Vpoint;
		double dot = nx * (vpnt->x - x1) + ny * (vpnt->y - y1) + nz * (vpnt->z - z1);

		if (dot > 1E-10) pnt->invalid = true;
		else 
		{
			pnt->invalid = false;
			num_valid_corner++;
		}
		pnt = pnt->nextcorner;
		if (pnt == fpnt) break;
	}

	if (num_valid_corner == f->num_corners) return false;

	if (num_valid_corner == 0)
	{
		delete_face_corners(f);
		return false;
	}

	// insert intersection corners
	double xx, yy, zz, xp, yp, zp;
	while (true)
	{
		pntp = pnt->nextcorner;
		if (pnt->invalid && !pntp->invalid || !pnt->invalid && pntp->invalid)
		{
			// insert intersection point
			vpnt = pnt->Vpoint; vpntp = pntp->Vpoint;
			xx = vpnt->x; yy = vpnt->y; zz = vpnt->z;
			xp = vpntp->x; yp = vpntp->y; zp = vpntp->z;

			double dot_1 = nx * (x1 - xx) + ny * (y1 - yy) + nz * (z1 - zz);
			double dot_2 = nx * (xp - xx) + ny * (yp - yy) + nz * (zp - zz);
			
			double u = dot_1 / dot_2;
			if (u > 1E-10 && 1.0 - u > 1E-10)
			{
				xx = xx + u * (xp - xx);
				yy = yy + u * (yp - yy);
				zz = zz + u * (zp - zz);

				if (pnt->invalid)
				{
					qx = xx; qy = yy; qz = zz;
					if (first) first = false;
				}
				else
				{
					px = xx; py = yy; pz = zz; 
					if (first) first = false;
				}

				vpnt = new VoronoiPoint(); 
				if (f->on_boundary)	vpnt->on_boundary = true;
				else vpnt->on_boundary = false;
				
				vpnt->x = xx;
				vpnt->y = yy;
				vpnt->z = zz;

				tmppnt = new FacePoint();
				tmppnt->invalid = false;
				tmppnt->nextcorner = pntp; pntp->prevcorner = tmppnt;
				tmppnt->prevcorner = pnt; pnt->nextcorner = tmppnt;
				tmppnt->Vpoint = vpnt;
				pnt = tmppnt;
				f->num_corners++;
			}
		}
		pnt = pnt->nextcorner;
		if (pnt == fpnt) break;
	}

	// remove invalid points	
	while (true)
	{
		pntp = pnt->nextcorner;
		if (pnt->invalid)
		{
			pntp = pnt->nextcorner;
			pntm = pnt->prevcorner;
			pntp->prevcorner = pntm; pntm->nextcorner = pntp;
			if (f->corner == pnt) 
			{
				f->corner = pntm;
				fpnt = pntm;
			}

			delete pnt->Vpoint;
			delete pnt;
			f->num_corners--;
			pnt = pntm;
		}

		pnt = pnt->nextcorner;

		if (pnt == fpnt && !pnt->invalid) break;
	}
	return true;
	#pragma endregion
}


inline bool MeshingGenie_3d::get_trimming_points(PolyFace* f,
	                  double x1, double y1, double z1,
					  double nx, double ny, double nz,
					  double &px, double &py, double &pz,
					  double &qx, double &qy, double &qz)
{
	#pragma region Trim Face using a plane:
	// mark face corners
	bool first(true);
	FacePoint* pnt; FacePoint* fpnt; VoronoiPoint* vpnt;
	FacePoint* pntp; VoronoiPoint* vpntp;
	pnt = f->corner; fpnt = pnt;
	size_t num_valid_corner(0);
	while (true)
	{
		// check if point is on the right side of plane
		vpnt = pnt->Vpoint;
		double dot = nx * (vpnt->x - x1) + ny * (vpnt->y - y1) + nz * (vpnt->z - z1);

		if (dot > 1E-10) pnt->invalid = true;
		else 
		{
			pnt->invalid = false;
			num_valid_corner++;
		}
		pnt = pnt->nextcorner;
		if (pnt == fpnt) break;
	}

	if (num_valid_corner == f->num_corners) return false;

	if (num_valid_corner == 0)
	{
		delete_face_corners(f);
		return false;
	}

	// insert intersection corners
	double xx, yy, zz, xp, yp, zp;
	while (true)
	{
		pntp = pnt->nextcorner;
		if (pnt->invalid && !pntp->invalid || !pnt->invalid && pntp->invalid)
		{
			// insert intersection point
			vpnt = pnt->Vpoint; vpntp = pntp->Vpoint;
			xx = vpnt->x; yy = vpnt->y; zz = vpnt->z;
			xp = vpntp->x; yp = vpntp->y; zp = vpntp->z;

			double dot_1 = nx * (x1 - xx) + ny * (y1 - yy) + nz * (z1 - zz);
			double dot_2 = nx * (xp - xx) + ny * (yp - yy) + nz * (zp - zz);
			
			double u = dot_1 / dot_2;
			if (u > 1E-10 && 1.0 - u > 1E-10)
			{
				xx = xx + u * (xp - xx);
				yy = yy + u * (yp - yy);
				zz = zz + u * (zp - zz);

				if (pnt->invalid)
				{
					qx = xx; qy = yy; qz = zz;
					if (first) first = false;
					else return true;
				}
				else
				{
					px = xx; py = yy; pz = zz; 
					if (first) first = false;
					else return true;
				}
				
			}
		}
		pnt = pnt->nextcorner;
		if (pnt == fpnt) break;
	}
	return false;
	#pragma endregion
}



inline void MeshingGenie_3d::weld_Voronoi_Point(VoronoiPoint* &vpoint, Polyhedron* poly, double tol)
{
	#pragma region Weld Voronoi Point:
	double xx = vpoint->x;
	double yy = vpoint->y;
	double zz = vpoint->z;

	size_t i, j, k, icell, jcell;
	i = size_t((xx - _xo) / _s);
	j = size_t((yy - _yo) / _s);
	k = size_t((zz - _zo) / _s);
	get_cell_index(i, j, k, icell);
	VoronoiPoint* vpnt(_cell_vpoints[icell]);
	// check point in that cell
	while (vpnt != 0)
	{
		if (fabs(vpnt->x - xx) < tol && fabs(vpnt->y - yy) < tol && fabs(vpnt->z - zz) < tol)	
		{
			bool weld_points(true);

			if (vpoint->on_boundary && !vpnt->on_boundary) weld_points = false;
				
			if (!vpoint->on_boundary && vpnt->on_boundary) weld_points = false;

			if (weld_points && vpoint->on_boundary) // both points on boundary
			{
				// vpoint must not hold any sharp features information (not a sharo corner nor on a sharp edge)				
			}

			if (fabs(vpnt->x - xx) < 1E-10 && fabs(vpnt->y - yy) < 1E-10 || fabs(vpnt->z - zz) < 1E-10) weld_points = true;

			if (weld_points)
			{
				delete vpoint;
				vpoint = vpnt;
				vpnt->polys.insert(poly);						
				return;
			}
		}
		vpnt = vpnt->nextVpoint;		
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

	VoronoiPoint* vvpnt;
	if (!ip && !jp && !kp && !im && !jm && !km) 
	{
		// no need to search neighborhoods, create a new Voronoi corner
		vpoint->nextVpoint = 0;
		vpoint->polys.insert(poly);
		vvpnt = _cell_vpoints[icell];
		if (vvpnt == 0) {_cell_vpoints[icell] = vpoint; return;}
		while (vvpnt->nextVpoint != 0)
		{			
			vvpnt = vvpnt->nextVpoint;
		}
		vvpnt->nextVpoint = vpoint;
		return;
	}

	// check points in the neighbor 27 cells
	for (size_t index = 0; index < 27; index++)
	{
		jcell = icell + _neighbors[index];
		// check all point in that cell
		vpnt = _cell_vpoints[jcell];
		while (vpnt != 0)
		{
			if (fabs(vpnt->x - xx) < tol && fabs(vpnt->y - yy) < tol && fabs(vpnt->z - zz) < tol) 
			{
				bool weld_points(true);

				if (vpoint->on_boundary && !vpnt->on_boundary) weld_points = false;
				
				if (!vpoint->on_boundary && vpnt->on_boundary) weld_points = false;

				if (weld_points && vpoint->on_boundary) // both points on boundary
				{
					// vpoint must not hold any sharp features information (not a sharo corner nor on a sharp edge)					
					
				}

				if (fabs(vpnt->x - xx) < 1E-10 && fabs(vpnt->y - yy) < 1E-10 || fabs(vpnt->z - zz) < 1E-10) weld_points = true;

				if (weld_points)
				{
					delete vpoint;
					vpoint = vpnt;
					vpnt->polys.insert(poly);						
					return;
				}
			}
			vpnt = vpnt->nextVpoint;		
		}
	}
	// Vornoi Point does not exist in the background grid, add it
	vpoint->nextVpoint = 0;
	vpoint->polys.insert(poly);
	vvpnt = _cell_vpoints[icell];
	if (vvpnt == 0) {_cell_vpoints[icell] = vpoint; return;}
	while (vvpnt->nextVpoint != 0)
	{			
		vvpnt = vvpnt->nextVpoint;
	}
	vvpnt->nextVpoint = vpoint;
	return;
	#pragma endregion
}


inline void MeshingGenie_3d::delete_poly_faces(Polyhedron* poly)
{
	#pragma region Delete Polyfaces:
	PolyFace* f = poly->face;
	while (f != 0)
	{
		PolyFace* nextf = f->nextface;
		delete_face_corners(f);
		delete f;
		f = nextf;					
	}
	#pragma endregion
}

inline MeshingGenie_3d::Polyhedron* MeshingGenie_3d::remove_external_polys(Polyhedron* poly)
{
	#pragma region Remove External polys:
	Polyhedron* ppoly(poly);
	Polyhedron* nextpoly(poly->nextpoly);
	Polyhedron* lastpoly(0);
	while (ppoly != 0)
	{
		nextpoly = ppoly->nextpoly;
		if (is_external_poly(ppoly))
		{
			if (poly == ppoly) poly = nextpoly;
			delete_poly_faces(ppoly);
			delete ppoly;
			if (lastpoly!=0) lastpoly->nextpoly = nextpoly;
		}
		else 
		{
			lastpoly = ppoly;
		}

		ppoly = nextpoly;
	}
	return poly;
	#pragma endregion
}

inline bool MeshingGenie_3d::is_external_poly(Polyhedron* poly)
{
	#pragma region Detect External poly:
	// check boundary faces
	PolyFace* f(poly->face);
	double xc(0.0), yc(0.0), zc(0.0);
	size_t num(0);
	while (f != 0)
	{
		size_t num_face_corners(f->num_corners);
		if (num_face_corners < 3)
		{
			f = f->nextface;
			continue;
		}

		FacePoint* pnt(f->corner);
		for (size_t ipnt = 0; ipnt < num_face_corners; ipnt++)
		{
			xc+= pnt->Vpoint->x;
			yc+= pnt->Vpoint->y;
			zc+= pnt->Vpoint->z; num++;
			pnt = pnt->nextcorner;
		}
		f = f->nextface;
	}
	xc /= num; yc /= num; zc /= num;

	// bounding box test
	if (xc < _box_xmin || xc > _box_xmax) return true;
	if (yc < _box_ymin || yc > _box_ymax) return true;
	if (zc < _box_zmin || zc > _box_zmax) return true;


	// Retrieve the cell of that point:
	size_t i, j, k, icell;
	i = size_t((xc - _xo) / _s);
	j = size_t((yc - _yo) / _s);
	k = size_t((zc - _zo) / _s);
	get_cell_index(i, j, k, icell);

	if (_external_cells[icell]) return true;
	if (!_boundary_cells[icell]) return false;

	size_t num_in(0), num_out(0); num = 0;
	std::set<size_t>* cell_faces(&_cell_faces[icell]);
	std::set<size_t>::iterator faces_iter;
	for (faces_iter = cell_faces->begin(); faces_iter != cell_faces->end(); faces_iter++)
	{
		size_t iface(*faces_iter);
		size_t ii = _bfaces[iface][0];
		double ax(xc - _xb[ii]), ay(yc - _yb[ii]), az(zc - _zb[ii]);
		double dot(_bfaces_nx[iface]*ax + _bfaces_ny[iface] * ay + _bfaces_nz[iface] * az);
		if (dot > 0) num_out++;
		else num_in++;
		num++;
	}

	if (num_out == num) return true;
	if (num_in == num) return false;

	Polyhedron* cell_poly = create_a_cell(icell);
	
	//if (fabs(xc + 1.338224) < 0.001 && fabs(yc + 2.189202) < 0.001 && fabs(zc + 4.76132309)< 0.0001)
	//{
	//	plot_polyhedron(poly);
	//	plot_polyhedron(cell_poly);
	//}

	for (faces_iter = cell_faces->begin(); faces_iter != cell_faces->end(); faces_iter++)
	{
		size_t iface(*faces_iter);
		size_t ii = _bfaces[iface][0];
		f = construct_super_triangle(_xb[ii], _yb[ii], _zb[ii], 
												  _bfaces_nx[iface], _bfaces_ny[iface], _bfaces_nz[iface],
												  xc, yc, zc, 0);
		if (!add_face_to_poly(cell_poly, f))
		{
			delete_face_corners(f); delete f;
		}
	}

	//if (fabs(xc + 1.338224) < 0.001 && fabs(yc + 2.189202) < 0.001 && fabs(zc + 4.76132309)< 0.0001)
	//{
	//	plot_polyhedron(cell_poly);
	//}

	f = cell_poly->face; num = 0;
	double xc_(0.0), yc_(0.0), zc_(0.0);
	while (f != 0)
	{
		size_t num_face_corners(f->num_corners);
		FacePoint* pnt(f->corner);
		for (size_t ipnt = 0; ipnt < num_face_corners; ipnt++)
		{
			xc_+= pnt->Vpoint->x;
			yc_+= pnt->Vpoint->y;
			zc_+= pnt->Vpoint->z; num++;
		}
		f = f->nextface;
	}
	xc_ /= num; yc_ /= num; zc_ /= num;

	// check the line segment connecting the two centers c and c_:
	size_t num_crossings(0);
	for (faces_iter = cell_faces->begin(); faces_iter != cell_faces->end(); faces_iter++)
	{
		size_t iface(*faces_iter);
		size_t ii = _bfaces[iface][0];
		double ax(xc - _xb[ii]), ay(yc - _yb[ii]), az(zc - _zb[ii]);
		double bx(xc - xc_), by(yc - yc_), bz(zc - zc_);

		double dot_a(_bfaces_nx[iface] * ax + _bfaces_ny[iface] * ay + _bfaces_nz[iface] * az);
		double dot_b(_bfaces_nx[iface] * bx + _bfaces_ny[iface] * by + _bfaces_nz[iface] * bz);
		
		if (dot_a < 1E-10) continue; // both points lies on a good side of that boundary face

		// intersection point
		double u = dot_a / dot_b;
		double xx = xc + u * (xc_ - xc);
		double yy = yc + u * (yc_ - yc);
		double zz = zc + u * (zc_ - zc);

		// make sure that this point lies in the boundary face
		if (point_in_boundary_face(iface, xx, yy, zz)) num_crossings++;
	}

	delete_poly_faces(cell_poly);
	delete cell_poly;

	if (num_crossings % 2 == 1) return true;
			
	return false;
	#pragma endregion
}

inline void MeshingGenie_3d::weld_poly_corners(Polyhedron* poly)
{
	#pragma region Weld Polyhedron Cornoers:
	while (poly != 0)
	{
		PolyFace* f = poly->face;
		while (f != 0)
		{
			#pragma region Weld face corners:
			if (f->num_corners >= 3)
			{
				FacePoint* pnt  = f->corner;
				FacePoint* fpnt = pnt;
				while (true)
				{					
					weld_Voronoi_Point(pnt->Vpoint, poly, _tol);				

					if (pnt->Vpoint == pnt->prevcorner->Vpoint)
					{
						// remove that point from this face
						pnt->prevcorner->nextcorner = pnt->nextcorner;
						pnt->nextcorner->prevcorner = pnt->prevcorner;
						f->num_corners--;
						FacePoint* tmp_pnt = pnt;
						pnt = pnt->prevcorner;
						delete tmp_pnt;
					}
					pnt = pnt->nextcorner;
					if (pnt == fpnt)
					{
						if (pnt->Vpoint == pnt->prevcorner->Vpoint)
						{
							pnt->prevcorner->nextcorner = pnt->nextcorner;
							pnt->nextcorner->prevcorner = pnt->prevcorner;
							f->num_corners--;
							f->corner = pnt->prevcorner;
							delete pnt;
						}
						break;
					}
				}
			}
			f = f->nextface;
			#pragma endregion
		}
		poly = poly->nextpoly;
	}
	#pragma endregion
}

inline void MeshingGenie_3d::remove_degenerate_faces(Polyhedron* poly)
{
	#pragma region Weld Polyhedron Cornoers:
	while (poly != 0)
	{
		PolyFace* f = poly->face;
		PolyFace* flast(0); 
		while (f != 0)
		{
			PolyFace* fn = f->nextface;
			if (f->num_corners < 3)
			{				
				if (f == poly->face) poly->face = fn;
				if (flast != 0) flast->nextface = fn;
				delete f;
			}
			else flast = f;
			f = fn;
		}
		poly = poly->nextpoly;
	}
	#pragma endregion
}

inline MeshingGenie_3d::Polyhedron* MeshingGenie_3d::create_a_cell(size_t icell)
{	
	#pragma region Create A Cell:
	size_t i, j, k;
	get_cell_indices(icell, i, j, k);
	double xo(_xo + i * _s), yo(_yo + j * _s), zo(_zo + k * _s);	
	
	PolyFace* fo; PolyFace* fn;
	PolyFace* f;
	FacePoint* po; FacePoint* p1; FacePoint* p2; FacePoint* p3;
	VoronoiPoint* vpnt;
	
	// i-face
	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo; vpnt->z = zo;
	po = new FacePoint(); po->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo; vpnt->z = zo + _s;
	p1 = new FacePoint(); p1->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo + _s; vpnt->z = zo + _s;
	p2 = new FacePoint(); p2->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo + _s; vpnt->z = zo;
	p3 = new FacePoint(); p3->Vpoint = vpnt;

	po->nextcorner = p1; p1->nextcorner = p2; p2->nextcorner = p3; p3->nextcorner = po;
	po->prevcorner = p3; p1->prevcorner = po; p2->prevcorner = p1; p3->prevcorner = p2;

	f = new PolyFace(); fo = f; fn = f;
	f->corner = po; f->nextface = 0; f->num_corners = 4; 
	f->nx = -1.0; f->ny = 0.0; f->nz = 0.0;

	// ip-face
	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo; vpnt->z = zo;
	po = new FacePoint(); po->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo; vpnt->z = zo + _s;
	p1 = new FacePoint(); p1->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo + _s; vpnt->z = zo + _s;
	p2 = new FacePoint(); p2->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo + _s; vpnt->z = zo;
	p3 = new FacePoint(); p3->Vpoint = vpnt;

	po->nextcorner = p3; p1->nextcorner = po; p2->nextcorner = p1; p3->nextcorner = p2;
	po->prevcorner = p1; p1->prevcorner = p2; p2->prevcorner = p3; p3->prevcorner = po;

	f = new PolyFace(); fn->nextface = f; fn = f;
	f->corner = po; f->nextface = 0; f->num_corners = 4; 
	f->nx = 1.0; f->ny = 0.0; f->nz = 0.0;

	// j-face
	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo; vpnt->z = zo;
	po = new FacePoint(); po->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo; vpnt->z = zo;
	p1 = new FacePoint(); p1->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo; vpnt->z = zo + _s;
	p2 = new FacePoint(); p2->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo; vpnt->z = zo + _s;
	p3 = new FacePoint(); p3->Vpoint = vpnt;

	po->nextcorner = p1; p1->nextcorner = p2; p2->nextcorner = p3; p3->nextcorner = po;
	po->prevcorner = p3; p1->prevcorner = po; p2->prevcorner = p1; p3->prevcorner = p2;

	f = new PolyFace(); fn->nextface = f; fn = f;
	f->corner = po; f->nextface = 0; f->num_corners = 4; 
	f->nx = 0.0; f->ny = -1.0; f->nz = 0.0;

	// jp-face
	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo + _s; vpnt->z = zo;
	po = new FacePoint(); po->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo + _s; vpnt->z = zo;
	p1 = new FacePoint(); p1->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo + _s; vpnt->z = zo + _s;
	p2 = new FacePoint(); p2->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo + _s; vpnt->z = zo + _s;
	p3 = new FacePoint(); p3->Vpoint = vpnt;

	po->nextcorner = p3; p1->nextcorner = po; p2->nextcorner = p1; p3->nextcorner = p2;
	po->prevcorner = p1; p1->prevcorner = p2; p2->prevcorner = p3; p3->prevcorner = po;

	f = new PolyFace(); fn->nextface = f; fn = f;
	f->corner = po; f->nextface = 0; f->num_corners = 4; 
	f->nx = 0.0; f->ny = 1.0; f->nz = 0.0;

	// k-face
	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo; vpnt->z = zo;
	po = new FacePoint(); po->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo + _s; vpnt->z = zo;
	p1 = new FacePoint(); p1->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo + _s; vpnt->z = zo;
	p2 = new FacePoint(); p2->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo; vpnt->z = zo;
	p3 = new FacePoint(); p3->Vpoint = vpnt;

	po->nextcorner = p1; p1->nextcorner = p2; p2->nextcorner = p3; p3->nextcorner = po;
	po->prevcorner = p3; p1->prevcorner = po; p2->prevcorner = p1; p3->prevcorner = p2;

	f = new PolyFace(); fn->nextface = f; fn = f;
	f->corner = po; f->nextface = 0; f->num_corners = 4; 
	f->nx = 0.0; f->ny = 0.0; f->nz = -1.0;

	// kp-face
	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo; vpnt->z = zo + _s;
	po = new FacePoint(); po->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo; vpnt->y = yo + _s; vpnt->z = zo + _s;
	p1 = new FacePoint(); p1->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo + _s; vpnt->z = zo + _s;
	p2 = new FacePoint(); p2->Vpoint = vpnt;

	vpnt = new VoronoiPoint(); vpnt->x = xo + _s; vpnt->y = yo; vpnt->z = zo + _s;
	p3 = new FacePoint(); p3->Vpoint = vpnt;

	po->nextcorner = p3; p1->nextcorner = po; p2->nextcorner = p1; p3->nextcorner = p2;
	po->prevcorner = p1; p1->prevcorner = p2; p2->prevcorner = p3; p3->prevcorner = po;

	f = new PolyFace(); fn->nextface = f; fn = f;
	f->corner = po; f->nextface = 0; f->num_corners = 4; 
	f->nx = 0.0; f->ny = 0.0; f->nz = 1.0;

	Polyhedron* poly = new Polyhedron();
	poly->face = fo; poly->nextpoly = 0;

	return poly;
	#pragma endregion
}

inline bool MeshingGenie_3d::add_face_to_poly(Polyhedron* poly, PolyFace* face)
{
	#pragma region Add a face to poly:
	PolyFace* f = poly->face;
	while (f != 0)
	{
		// trim face with internal faces only listed first
		double rx, ry, rz, sx, sy, sz;
		if (face->num_corners > 0 && f->num_corners > 0) 
		{
			trim_face(face, f->corner->Vpoint->x, f->corner->Vpoint->y, f->corner->Vpoint->z, f->nx, f->ny, f->nz,
						rx, ry, rz, sx, sy, sz);
			if (face->num_corners < 3) return false;
		}
		f = f->nextface;
	}
	if (face->num_corners < 3) return false;
	f = poly->face;
	PolyFace* fn;
	while (f != 0)
	{
		// trim old faces using new face
		double rx, ry, rz, sx, sy, sz;
		if (f->num_corners > 2) 
		{						
			trim_face(f, face->corner->Vpoint->x, face->corner->Vpoint->y, face->corner->Vpoint->z, face->nx, face->ny, face->nz,
							rx, ry, rz, sx, sy, sz);
			if (face->num_corners < 3) break;
		}
		fn = f;
		f = f->nextface;
	}
	fn->nextface = face;
	return true;
	#pragma endregion
}

inline MeshingGenie_3d::Polyhedron* MeshingGenie_3d::copy_polyhedron(Polyhedron* poly)
{
	#pragma region Create a copy of a polyhedron:
	Polyhedron* newpoly = new Polyhedron();
	newpoly->nextpoly = 0; 
	PolyFace* f = poly->face;
	PolyFace* newf;  PolyFace* fn(0);
	while (f != 0)
	{		
		newf = copy_polyface(f, false);
		if (fn == 0)
		{
			newpoly->face = newf;
			fn = newf;
		}
		else
		{
			fn->nextface = newf;
			fn = newf;
		}
		f = f->nextface;
	}	
	return newpoly;
	#pragma endregion
}

inline MeshingGenie_3d::PolyFace* MeshingGenie_3d::copy_polyface(PolyFace* f, bool invert_normals)
{
	#pragma region Create A copy of a PolyFace:
	PolyFace* newf = new PolyFace();

	newf->q = f->q;
	newf->preserved = f->preserved;
	newf->on_boundary = f->on_boundary;
	newf->nx = f->nx;
	newf->ny = f->ny;
	newf->nz = f->nz;
	newf->num_corners = f->num_corners;
	newf->inverted = false;
	if (invert_normals) newf->inverted = true;
	newf->nextface = 0;
				
	if (invert_normals)
	{
		newf->nx = -newf->nx;
		newf->ny = -newf->ny;
		newf->nz = -newf->nz;				
	}

	size_t num_corners(f->num_corners);

	if (num_corners == 0) 
	{
		newf->corner = 0;
		return newf;
	}

	std::vector<FacePoint*> pnts(num_corners);

	FacePoint* pnt(f->corner);
	for (size_t i = 0; i < num_corners; i++)
	{
		pnts[i] = new FacePoint();
		pnts[i]->Vpoint = new VoronoiPoint();
		pnts[i]->Vpoint->x = pnt->Vpoint->x;
		pnts[i]->Vpoint->y = pnt->Vpoint->y;
		pnts[i]->Vpoint->z = pnt->Vpoint->z;
		pnts[i]->Vpoint->on_boundary = pnt->Vpoint->on_boundary;
		pnt = pnt->nextcorner;
	}

	for (size_t i = 0; i < num_corners; i++)
	{
		int im, ip(i + 1);
		if (i == 0) im = num_corners - 1;
		else im = i - 1;
		if (ip == num_corners) ip = 0;

		if (invert_normals)
		{
			pnts[i]->nextcorner = pnts[im];
			pnts[i]->prevcorner = pnts[ip];
		}
		else
		{
			pnts[i]->nextcorner = pnts[ip];
			pnts[i]->prevcorner = pnts[im];
		}
	}
	newf->corner = pnts[0];
	return newf; 
	#pragma endregion
}

inline MeshingGenie_3d::PolyFace* MeshingGenie_3d::construct_super_triangle(
	                                               double xmid, double ymid, double zmid, double nx, double ny, double nz,
				                                   double xc, double yc, double zc, Point* q) // center of a cell
{
	#pragma region Construct Super Triangle:
	double mx, my, mz, lx, ly, lz;

	double x1, y1, z1, x2, y2, z2, x3, y3, z3;

	if (fabs(nx) > fabs(ny) && fabs(nx) > fabs(nz))
	{
		// plane is close to yz
		mx = -ny; my = nx; mz = 0;
	}
	else if (fabs(ny) > fabs(nz))
	{
		// plane is close to zx
		mx = 0; my = -nz; mz = ny;
	}
	else
	{
		// plane is close to xy
		mx = nz; my = 0; mz = -nx;
	}

	// normalize n and m and l
	double m = mx * mx + my * my + mz * mz;
	double m_inv = 1.0 / sqrt(m);
	mx *= m_inv; my *= m_inv; mz *= m_inv;

	double n = nx * nx + ny * ny + nz * nz;
	double n_inv = 1.0 / sqrt(n);
	nx *= n_inv; ny *= n_inv; nz *= n_inv;			

	lx = ny * mz - nz * my;
	ly = nz * mx - nx * mz;
	lz = nx * my - ny * mx;

	double l = lx * lx + ly * ly + lz * lz;
	double l_inv = 1.0 / sqrt(l);
	lx *= l_inv; ly *= l_inv; lz *= l_inv;	

	x1 = xmid + 100 * mx + 100 * lx;
	y1 = ymid + 100 * my + 100 * ly;
	z1 = zmid + 100 * mz + 100 * lz;			

	x3 = xmid + 100 * mx - 100 * lx;
	y3 = ymid + 100 * my - 100 * ly;
	z3 = zmid + 100 * mz - 100 * lz;

	x2 = xmid - 100 * mx;
	y2 = ymid - 100 * my;
	z2 = zmid - 100 * mz;

	// first point
	FacePoint* pnt = new FacePoint();
	VoronoiPoint* vpnt = new VoronoiPoint(); vpnt->on_boundary = false;
	vpnt->x = x1; vpnt->y = y1; vpnt->z = z1;
	pnt->Vpoint = vpnt; FacePoint* fpnt = pnt;

	// second point
	pnt = new FacePoint();
	vpnt = new VoronoiPoint(); vpnt->on_boundary = false;			
	vpnt->x = x2; vpnt->y = y2; vpnt->z = z2;
	pnt->Vpoint = vpnt; pnt->prevcorner = fpnt; fpnt->nextcorner = pnt; 

	// third point
	pnt = new FacePoint();
	vpnt = new VoronoiPoint(); vpnt->on_boundary = false;			
	vpnt->x = x3; vpnt->y = y3; vpnt->z = z3;
	pnt->Vpoint = vpnt; pnt->nextcorner = fpnt; fpnt->prevcorner = pnt;
	pnt->prevcorner = fpnt->nextcorner; fpnt->nextcorner->nextcorner = pnt; 

	PolyFace* face = new PolyFace(); face->on_boundary = false; face->inverted = false;
	face->corner = fpnt; face->nextface = 0;
	face->nx = nx; face->ny = ny; face->nz = nz;
	face->num_corners = 3; face->q = q;
	face->preserved = false;

	double S = 4.5 * _s;
			
	double px, py, pz, qx, qy, qz;
	trim_face(face, xc + S, yc, zc,  1.0, 0.0, 0.0, px, py, pz, qx, qy, qz);
	if (face->num_corners == 0) return face;

	trim_face(face, xc - S, yc, zc, -1.0, 0.0, 0.0, px, py, pz, qx, qy, qz);
	if (face->num_corners == 0) return face;

	trim_face(face, xc, yc + S, zc,  0.0,  1.0, 0.0, px, py, pz, qx, qy, qz);
	if (face->num_corners == 0) return face;

	trim_face(face, xc, yc - S, zc,  0.0, -1.0, 0.0, px, py, pz, qx, qy, qz);
	if (face->num_corners == 0) return face;

	trim_face(face, xc, yc, zc + S,  0.0, 0.0,  1.0, px, py, pz, qx, qy, qz);
	if (face->num_corners == 0) return face;

	trim_face(face, xc, yc, zc - S,  0.0, 0.0, -1.0, px, py, pz, qx, qy, qz);
	if (face->num_corners == 0) return face;

	return face;
	#pragma endregion
}


inline bool MeshingGenie_3d::point_in_boundary_face(size_t iface, double px, double py, double pz)
{
	#pragma region point in boundary face check:
	double nx(_bfaces_nx[iface]), ny(_bfaces_ny[iface]), nz(_bfaces_nz[iface]);
	size_t num_face_corners(_bfaces[iface].size());
	for (size_t j = 0; j < num_face_corners; j++)
	{
		size_t jp(j + 1);
		if (jp == num_face_corners) jp = 0;
					
		size_t node_i = _bfaces[iface][j];
		size_t node_j = _bfaces[iface][jp];

		double xi(_xb[node_i]), yi(_yb[node_i]), zi(_zb[node_i]);
		double xj(_xb[node_j]), yj(_yb[node_j]), zj(_zb[node_j]);

		double ax(xj - xi), ay(yj - yi), az(zj - zi);
		double mx(ay * nz - az * ny);
		double my(az * nx - ax * nz);
		double mz(ax * ny - ay * nx);

		double dot(mx * (px - xi) + my * (py - yi) + mz * (pz - zi));
		if (dot > 1E-10)
		{
			return false; 
		}
	}
	return true;
	#pragma endregion
}


inline void MeshingGenie_3d::delete_face_corners(PolyFace* f)
{
	#pragma region Delete  Face Corners:
	size_t num_face_corners(f->num_corners);
	FacePoint* pnt;
	for (size_t icorner = 0; icorner < num_face_corners; icorner++)
	{
		pnt = f->corner->nextcorner;
		delete f->corner->Vpoint;
		delete f->corner;
		f->corner = pnt;
	}
	f->num_corners = 0;
	f->corner = 0;
	#pragma endregion
}


inline MeshingGenie_3d::PolyFace* MeshingGenie_3d::create_boundary_face(size_t iface)
{
	#pragma region Create an instance of a boundary Face:
	size_t num_face_corners(_bfaces[iface].size());
	FacePoint** pnts= new FacePoint*[num_face_corners];
	for (size_t i = 0; i < num_face_corners; i++)
	{
		size_t node_i = _bfaces[iface][i];
		double xi(_xb[node_i]), yi(_yb[node_i]), zi(_zb[node_i]);

		pnts[i] = new FacePoint();
		VoronoiPoint* vpnt = new VoronoiPoint();
		vpnt->x = xi; vpnt->y = yi; vpnt->z = zi;
		pnts[i]->Vpoint = vpnt;
	}

	for (size_t i = 0; i < num_face_corners; i++)
	{
		size_t ip = i + 1;
		if (ip == num_face_corners) ip = 0;
		pnts[i]->nextcorner = pnts[ip];
		pnts[ip]->prevcorner = pnts[i];
	}

	PolyFace* face = new PolyFace();
	face->nx = _bfaces_nx[iface]; 
	face->ny = _bfaces_ny[iface]; 
	face->nz = _bfaces_nz[iface]; 
	face->corner = pnts[0];
	face->num_corners = num_face_corners;
	face->on_boundary = true;
	face->preserved = false;
	face->q = 0;
	face->nextface = 0;
	face->inverted = false;
	return face;
	#pragma endregion
}

inline MeshingGenie_3d::PolyFace* MeshingGenie_3d::create_internal_face(size_t iface)
{
	#pragma region Create an instance of an internal Face:
	size_t num_face_corners(_xf[iface].size());
	FacePoint** pnts= new FacePoint*[num_face_corners];
	for (size_t i = 0; i < num_face_corners; i++)
	{
		double xi(_xf[iface][i]), yi(_yf[iface][i]), zi(_zf[iface][i]);

		pnts[i] = new FacePoint();
		VoronoiPoint* vpnt = new VoronoiPoint();
		vpnt->x = xi; vpnt->y = yi; vpnt->z = zi;
		pnts[i]->Vpoint = vpnt;
	}

	for (size_t i = 0; i < num_face_corners; i++)
	{
		size_t ip = i + 1;
		if (ip == num_face_corners) ip = 0;
		pnts[i]->nextcorner = pnts[ip];
		pnts[ip]->prevcorner = pnts[i];
	}

	PolyFace* face = new PolyFace();
	face->nx = _nfx[iface]; 
	face->ny = _nfy[iface]; 
	face->nz = _nfz[iface]; 
	face->corner = pnts[0];
	face->num_corners = num_face_corners;
	face->on_boundary = true;
	face->preserved = false;
	face->q = 0;
	face->nextface = 0;
	face->inverted = false;
	return face;
	#pragma endregion
}

void MeshingGenie_3d::sample_face(size_t iface)
{
	#pragma region Sample Face:

	size_t num_face_cracks(_xfc[iface].size());
	size_t num_face_corners(_xf[iface].size());

	double nx(_nfx[iface]), ny(_nfy[iface]), nz(_nfz[iface]);

	double nxy = sqrt(nx * nx + ny * ny);
	double cos_phi(1.0); double sin_phi(0.0);
	if (fabs(nxy) > 1E-10) 
	{
		cos_phi = nx / nxy; sin_phi = ny / nxy;
	}
	

	double mx = nx * cos_phi + ny * sin_phi;
	double my = -nx * sin_phi + ny * cos_phi;
	double mz = nz;

	for (size_t inode = 0; inode < num_face_corners; inode++)
	{
		// rotate around the z-axis -phi
		double xx(_xf[iface][inode]), yy(_yf[iface][inode]), zz(_zf[iface][inode]);
		_xf[iface][inode] =  xx * cos_phi + yy * sin_phi;
		_yf[iface][inode] = -xx * sin_phi + yy * cos_phi;
	}

	for (size_t icrack = 0; icrack < num_face_cracks; icrack++)
	{
		size_t num_crack_nodes(_xfc[iface][icrack].size());
		for (size_t inode = 0; inode < num_crack_nodes; inode++)
		{
			// rotate around the z-axis -phi
			double xx(_xfc[iface][icrack][inode]), yy(_yfc[iface][icrack][inode]), zz(_zfc[iface][icrack][inode]);
			_xfc[iface][icrack][inode] =  xx * cos_phi + yy * sin_phi;
			_yfc[iface][icrack][inode] = -xx * sin_phi + yy * cos_phi;
		}
	}

	double mxz = sqrt(mx * mx + mz * mz);
	double cos_theta(mz / mxz), sin_theta(mx / mxz);

	double oz =  mz * cos_theta + mx * sin_theta;
	double ox = -mz * sin_theta + mx * cos_theta;
	
	for (size_t inode = 0; inode < num_face_corners; inode++)
	{
		// rotate around the y-axis theta
		double xx(_xf[iface][inode]), yy(_yf[iface][inode]), zz(_zf[iface][inode]);
		_zf[iface][inode] =  zz * cos_theta + xx * sin_theta;
		_xf[iface][inode] = -zz * sin_theta + xx * cos_theta;
	}

	for (size_t icrack = 0; icrack < num_face_cracks; icrack++)
	{
		size_t num_crack_nodes(_xfc[iface][icrack].size());
		for (size_t inode = 0; inode < num_crack_nodes; inode++)
		{
			// rotate around the z-axis -phi
			double xx(_xfc[iface][icrack][inode]), yy(_yfc[iface][icrack][inode]), zz(_zfc[iface][icrack][inode]);
			_zfc[iface][icrack][inode] =  zz * cos_theta + xx * sin_theta;
			_xfc[iface][icrack][inode] = -zz * sin_theta + xx * cos_theta;
		}
	}

	double DZ = _zf[iface][0];

	std::vector<double> ExtBound; ExtBound.reserve(num_face_corners * 2);
	std::vector< std::vector<double> > Holes;
	std::vector< std::vector<double> > Cracks(num_face_cracks);
		
	for (size_t inode = 0; inode < num_face_corners; inode++)
	{
		ExtBound.push_back(_xf[iface][inode]);
		ExtBound.push_back(_yf[iface][inode]);
	}

	for (size_t icrack = 0; icrack < num_face_cracks; icrack++)
	{
		size_t num_crack_nodes(_xfc[iface][icrack].size());
		for (size_t inode = 0; inode < num_crack_nodes; inode++)
		{
			Cracks[icrack].push_back(_xfc[iface][icrack][inode]);
			Cracks[icrack].push_back(_yfc[iface][icrack][inode]);
		}
	}


	MeshingGenie_2d genie = MeshingGenie_2d(_dm / sqrt(2.0), ExtBound, Holes, Cracks, 0, false); 	
	
	//genie.use_fixed_seed(atoi(argv[2])); //nonzero to use a fixed seed
	genie.use_fixed_seed(_fixed_seed); //nonzero to use a fixed seed
	genie.execute();

	std::vector<double> xs; std::vector<double> ys; std::vector<double> zs;

	genie.get_point_cloud(xs, ys);

	// move samples back to the 3D domain
	size_t number_samples(xs.size()); zs.resize(number_samples);
	for (size_t inode = 0; inode < number_samples; inode++)
	{
		// rotate around the y-axis theta
		double xx(xs[inode]), yy(ys[inode]), zz(DZ);
		zs[inode] =  zz * cos_theta - xx * sin_theta;
		xs[inode] =  zz * sin_theta + xx * cos_theta;
	}
	for (size_t inode = 0; inode < num_face_corners; inode++)
	{
		// rotate around the y-axis theta
		double xx(_xf[iface][inode]), yy(_yf[iface][inode]), zz(_zf[iface][inode]);
		_zf[iface][inode] = zz * cos_theta - xx * sin_theta;
		_xf[iface][inode] = zz * sin_theta + xx * cos_theta;
	}
	for (size_t icrack = 0; icrack < num_face_cracks; icrack++)
	{
		size_t num_crack_nodes(_xfc[iface][icrack].size());
		for (size_t inode = 0; inode < num_crack_nodes; inode++)
		{
			// rotate around the z-axis -phi
			double xx(_xfc[iface][icrack][inode]), yy(_yfc[iface][icrack][inode]), zz(_zfc[iface][icrack][inode]);
			_zfc[iface][icrack][inode] =  zz * cos_theta - xx * sin_theta;
			_xfc[iface][icrack][inode] =  zz * sin_theta + xx * cos_theta;
		}
	}

	for (size_t inode = 0; inode < number_samples; inode++)
	{
		// rotate around the z-axis -phi
		double xx(xs[inode]), yy(ys[inode]), zz(zs[inode]);
		xs[inode] =  xx * cos_phi - yy * sin_phi;
		ys[inode] =  xx * sin_phi + yy * cos_phi;

		//std::cout << xs[inode] << " " << ys[inode] << " " << zs[inode] << std::endl;
	}

	for (size_t inode = 0; inode < num_face_corners; inode++)
	{
		// rotate around the z-axis -phi
		double xx(_xf[iface][inode]), yy(_yf[iface][inode]), zz(_zf[iface][inode]);
		_xf[iface][inode] =  xx * cos_phi - yy * sin_phi;
		_yf[iface][inode] =  xx * sin_phi + yy * cos_phi;
	}

	for (size_t icrack = 0; icrack < num_face_cracks; icrack++)
	{
		size_t num_crack_nodes(_xfc[iface][icrack].size());
		for (size_t inode = 0; inode < num_crack_nodes; inode++)
		{
			// rotate around the z-axis -phi
			double xx(_xfc[iface][icrack][inode]), yy(_yfc[iface][icrack][inode]), zz(_zfc[iface][icrack][inode]);
			_xfc[iface][icrack][inode] =  xx * cos_phi - yy * sin_phi;
			_yfc[iface][icrack][inode] =  xx * sin_phi + yy * cos_phi;
		}
	}

	for (size_t inode = 0; inode < number_samples; inode++)
	{
		Point* newpoint = closest_point(xs[inode], ys[inode], zs[inode], 1E-10);
		if (newpoint == 0)
		{
			newpoint = new Point();
			newpoint->x = xs[inode];
			newpoint->y = ys[inode];
			newpoint->z = zs[inode];
			newpoint->_on_reflex_corner = false;
			newpoint->_on_internal_boundaries = true;
			newpoint->next = 0;
			newpoint->poly = 0;

			size_t i, j, k, icell;
			i = size_t((newpoint->x - _xo) / _s);
			j = size_t((newpoint->y - _yo) / _s);
			k = size_t((newpoint->z - _zo) / _s);
			get_cell_index(i, j, k, icell);

			if (_cell_points[icell] == 0) _cell_points[icell] = newpoint;
			else
			{
				Point* q = _cell_points[icell];
				while (q->next != 0) q = q->next;
				q->next = newpoint;				
			}
			_num_inserted_points++;

			if (!_invalid_cells[icell]) 
			{
				_invalid_cells[icell] = true; _num_valid--;
			}

			_cell_faces_iter = _cell_internal_faces.find(icell);
			if (_cell_faces_iter == _cell_internal_faces.end())
			{
				std::set<size_t> newset;
				newset.insert(iface);
				_cell_internal_faces[icell] = newset;
			}
			else _cell_faces_iter->second.insert(iface);
		}
	}
	#pragma endregion
}

void MeshingGenie_3d::detect_internal_cracks_on_internal_faces()
{
	#pragma region Detect Internal Cracks on Internal Faces:
	size_t num_internal_faces(_xf.size());
	_xfc.clear(); _yfc.clear(); _zfc.clear();
	_xfc.resize(num_internal_faces);
	_yfc.resize(num_internal_faces);
	_zfc.resize(num_internal_faces);
	for (size_t iface = 0; iface < num_internal_faces; iface++)
	{
		PolyFace* fi = create_internal_face(iface);

		for (size_t jface = iface + 1; jface < num_internal_faces; jface++)
		{
			PolyFace* fj = create_internal_face(jface);

			double px, py, pz, qx, qy, qz, rx, ry, rz, sx, sy, sz;

			if (get_trimming_points(fi, fj->corner->Vpoint->x, fj->corner->Vpoint->y, fj->corner->Vpoint->z, fj->nx, fj->ny, fj->nz,
				                    px, py, pz, qx, qy, qz) && 
				get_trimming_points(fj, fi->corner->Vpoint->x, fi->corner->Vpoint->y, fi->corner->Vpoint->z, fi->nx, fi->ny, fi->nz,
				                    rx, ry, rz, sx, sy, sz))
			{
				
				bool use_x(false), use_y(false), use_z(false);
				double pq_x(fabs(px - qx)), pq_y(fabs(py - qy)), pq_z(fabs(pz - qz));
				if (pq_x > pq_y && pq_x > pq_z) use_x = true;
				else if (pq_y > pq_z && pq_y > pq_x) use_y = true;
				else use_z = true;

				bool p_is_r(false), p_is_s(false), q_is_r(false), q_is_s(false);
				if (fabs(px - rx) < 1E-10 && fabs(py - ry) < 1E-10 && fabs(pz - rz) < 1E-10) p_is_r = true;
				if (fabs(px - sx) < 1E-10 && fabs(py - sy) < 1E-10 && fabs(pz - sz) < 1E-10) p_is_s = true;
				if (fabs(qx - rx) < 1E-10 && fabs(qy - ry) < 1E-10 && fabs(qz - rz) < 1E-10) q_is_r = true;
				if (fabs(qx - sx) < 1E-10 && fabs(qy - sy) < 1E-10 && fabs(qz - sz) < 1E-10) q_is_s = true;

				bool p_in_rs(false), q_in_rs(false), r_in_pq(false), s_in_pq(false);
				
				if (use_x)
				{
					if (!p_is_r && !p_is_s && ((px - rx >= -1E-10 && sx - px >= -1E-10) || (rx - px >= -1E-10 && px - sx >= -1E-10))) p_in_rs = true;
					if (!q_is_r && !q_is_s && ((qx - rx >= -1E-10 && sx - qx >= -1E-10) || (rx - qx >= -1E-10 && qx - sx >= -1E-10))) q_in_rs = true;
					if (!p_is_r && !q_is_r && ((px - rx >= -1E-10 && rx - qx >= -1E-10) || (rx - px >= -1E-10 && qx - rx >= -1E-10))) r_in_pq = true;
					if (!p_is_s && !q_is_s && ((px - sx >= -1E-10 && sx - qx >= -1E-10) || (sx - px >= -1E-10 && qx - sx >= -1E-10))) s_in_pq = true;
				}
				else if (use_y)
				{
					if (!p_is_r && !p_is_s && ((py - ry >= -1E-10 && sy - py >= -1E-10) || (ry - py >= -1E-10 && py - sy >= -1E-10))) p_in_rs = true;
					if (!q_is_r && !q_is_s && ((qy - ry >= -1E-10 && sy - qy >= -1E-10) || (ry - qy >= -1E-10 && qy - sy >= -1E-10))) q_in_rs = true;
					if (!p_is_r && !q_is_r && ((py - ry >= -1E-10 && ry - qy >= -1E-10) || (ry - py >= -1E-10 && qy - ry >= -1E-10))) r_in_pq = true;
					if (!p_is_s && !q_is_s && ((py - sy >= -1E-10 && sy - qy >= -1E-10) || (sy - py >= -1E-10 && qy - sy >= -1E-10))) s_in_pq = true;				
				}
				else
				{				
					if (!p_is_r && !p_is_s && ((pz - rz >= -1E-10 && sz - pz >= -1E-10) || (rz - pz >= -1E-10 && pz - sz >= -1E-10))) p_in_rs = true;
					if (!q_is_r && !q_is_s && ((qz - rz >= -1E-10 && sz - qz >= -1E-10) || (rz - qz >= -1E-10 && qz - sz >= -1E-10))) q_in_rs = true;
					if (!p_is_r && !q_is_r && ((pz - rz >= -1E-10 && rz - qz >= -1E-10) || (rz - pz >= -1E-10 && qz - rz >= -1E-10))) r_in_pq = true;
					if (!p_is_s && !q_is_s && ((pz - sz >= -1E-10 && sz - qz >= -1E-10) || (sz - pz >= -1E-10 && qz - sz >= -1E-10))) s_in_pq = true;
				}

				std::vector< double > ex; std::vector< double > ey; std::vector< double > ez;
				if      (p_in_rs && q_in_rs) sample_edge(px, py, pz, qx, qy, qz, 0.5 * _dm, ex, ey, ez);
				else if (p_in_rs && r_in_pq) sample_edge(px, py, pz, rx, ry, rz, 0.5 * _dm, ex, ey, ez);
				else if (p_in_rs && s_in_pq) sample_edge(px, py, pz, sx, sy, sz, 0.5 * _dm, ex, ey, ez);
				else if (q_in_rs && r_in_pq) sample_edge(qx, qy, qz, rx, ry, rz, 0.5 * _dm, ex, ey, ez);
				else if (q_in_rs && s_in_pq) sample_edge(qx, qy, qz, sx, sy, sz, 0.5 * _dm, ex, ey, ez);
				else if (r_in_pq && s_in_pq) sample_edge(rx, ry, rz, sx, sy, sz, 0.5 * _dm, ex, ey, ez);
				
				if (ex.size() > 0)
				{
					_xfc[iface].push_back(ex); _yfc[iface].push_back(ey); _zfc[iface].push_back(ez);
					_xfc[jface].push_back(ex); _yfc[jface].push_back(ey); _zfc[jface].push_back(ez);					
				}
			}			
		}
	}
	#pragma endregion
}
