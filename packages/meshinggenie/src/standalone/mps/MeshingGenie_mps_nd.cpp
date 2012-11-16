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
 * Author: Mohamed S. Ebeida (msebeid@sandia.gov)
 * Description: 
      This class is based on the method presented in the following article
      M. S. Ebeida, S. A. Mitchell, A. Patney, A. A. Davidson, and J. D. Owens, 
      "A simple algorithm for maximal Poisson-disk sampling in high dimensions", 
	  Computer Graphics Forum (Eurographics 2012), 31(2), May 2012.
   
   input: a set of conforming simplices defining the boundaries of a (non-convex) domain
          + a ditribution radius (greater than the smallest feature size)
   output: a maximal Poisson-disk sample

 * The Random Number Generator is provided by George Marsaglia available at
   http://www.velocityreviews.com/forums/t720512-re-rngs-a-double-kiss.html

 * Last modified: 11/16/2012
********************************************************************************/

// Todo List:
// 1. Add active cells from boundary cells
// 2. Improve boundary face sampling by refining boundary faces and sampling non-covered only
// 3. Prepare output, nearest neighbor queries
// 4. Perform timing results


#include "MeshingGenie_mps_nd.h"

#include "MeshingGenie_plotter_2d.h"

MeshingGenie_plotter_2d _plotter;

int MeshingGenie_mps_nd::solve_mps(size_t ndim, double r, size_t random_seed, 
			                       std::vector< std::vector<double> > &boundary_points,
			                       std::vector< std::vector<size_t> > &boundary_faces, 
				                   std::vector< std::vector<double> > &sample_points)
{
	_ndim = ndim; _r = r; _rsq = _r * _r; _random_seed = random_seed;
	_boundary_points = &boundary_points;
	_boundary_faces = &boundary_faces;
	_sample_points = &sample_points;

	_num_active_cells = 0; _num_inserted_points = 0; _num_darts = 0;
	_current_dart = 0;

	clock_t start_time, end_time; double cpu_time;	 
	start_time = clock();

	init();
	
	initiate_random_generator(random_seed);

	generate_bg_grid();

	create_generic_neighbor_list();

	initiate_active_pool();

	identify_boundary_cells();

	identify_exterior_and_interior_cells();

	_boundary_cell_points.resize(_boundary_cells.size());
	_interior_cell_points.resize(_interior_cells.size());

	for (size_t iref = 0; iref < 30; iref++)
	{
		throw_darts(iref);
		if (iref == 0 && _num_active_cells == 0) 
		{			
			break;
		}

		if (_num_active_cells == 0)
		{
			// destroy parents of the previous level
			size_t num_active_cells = _active_cells.size();
			for (size_t icell = 0; icell < num_active_cells; icell++)
			{
				delete[] _active_cells[icell];
			}
			break;
		}
	}
	end_time = clock();

	// Report Number of points and time consumed

	cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

	std::cout<< "======================= MeshingGenie =============================" << std::endl;
	std::cout<< "Distribution Radius = " << _r << std::endl;
	std::cout<< "Number of inserted points = " << _num_inserted_points << std::endl;
	std::cout<< "Number of thrown darts = " << _num_darts << std::endl;
	std::cout<< "Execution Time = " << cpu_time << " seconds." << std::endl;
	std::cout<< "==================================================================" << std::endl;


	return 0;
}

void MeshingGenie_mps_nd::init()
{
	#pragma region Initialize Arrays:
	// bavkground grid data
	_xmin = new double[_ndim];
	_xmax = new double[_ndim];
	_nc = new size_t [_ndim];
	#pragma endregion
}

void MeshingGenie_mps_nd::destroy_mps_data()
{
	#pragma region Clean up dynamic arrays:
	delete[] _xmin;	delete[] _xmax;	delete[] _nc;

	clear_cell_pool(_boundary_cells);
	clear_cell_pool(_exterior_cells);
	clear_cell_pool(_interior_cells);
	clear_cell_pool(_covered_cells);

	clear_cell_points(_boundary_cell_points);
	clear_cell_points(_interior_cell_points);

	for (size_t icell = 0; icell< _num_neighbor_cells; icell++)
	{
		delete[] _generic_neighbor_cells[icell];
	}
	delete[] _generic_neighbor_cells;

	_active_cells.clear();
	#pragma endregion
}

inline void MeshingGenie_mps_nd::clear_cell_pool(std::vector<size_t*> &pool)
{
	#pragma region Clear Cell Pool:
	size_t num_cells(pool.size());
	for (size_t icell = 0; icell < num_cells; icell++)
	{
		delete[] pool[icell];
	}
	pool.clear();
	#pragma endregion
}

inline void MeshingGenie_mps_nd::clear_cell_points(std::vector<double*> &cell_points)
{
	#pragma region Clear Cell points:
	size_t num_cells(cell_points.size());
	for (size_t icell = 0; icell < num_cells; icell++)
	{
		if (cell_points[icell] == 0) continue;
		delete[] cell_points[icell];
	}
	cell_points.clear();
	#pragma endregion
}


void MeshingGenie_mps_nd::initiate_random_generator(unsigned long x)
{
	#pragma region initiate The Random generator:
	//assert(sizeof (double) >= 54) ;

	cc = 1.0 / 9007199254740992.0; // inverse of 2^53rd power
	size_t i;
	qlen = indx = sizeof Q / sizeof Q[0];
	for (i = 0; i < qlen; i++) Q[i] = 0;

	c = 0.0, zc = 0.0,	/* current CSWB and SWB `borrow` */
	zx = 5212886298506819.0 / 9007199254740992.0,	/* SWB seed1 */
	zy = 2020898595989513.0 / 9007199254740992.0;	/* SWB seed2 */
	
	size_t j;
	double s, t;	 /* Choose 32 bits for x, 32 for y */
	if (x == 0) x = (size_t) time(NULL); //123456789; /* default seeds */
	unsigned long y = 362436069; /* default seeds */

	/* Next, seed each Q[i], one bit at a time, */
	for (i = 0; i < qlen; i++)
	{ /* using 9th bit from Cong+Xorshift */
		s = 0.0;
		t = 1.0;
		for (j = 0; j < 52; j++)
		{
			t = 0.5 * t; /* make t=.5/2^j */
			x = 69069 * x + 123;
			y ^= (y << 13);
			y ^= (y >> 17);
			y ^= (y << 5);
			if (((x + y) >> 23) & 1) s = s + t; /* change bit of s, maybe */
		}	 /* end j loop */
		Q[i] = s;
	} /* end i seed loop, Now generate 10^9 dUNI's: */
	#pragma endregion
}

double MeshingGenie_mps_nd::generate_a_random_number()
{ 
	#pragma region generate a random number:
	/* Takes 14 nanosecs, Intel Q6600,2.40GHz */
	int i, j;
	double t; /* t: first temp, then next CSWB value */
	/* First get zy as next lag-2 SWB */
	t = zx - zy - zc;
	zx = zy;
	if (t < 0)
	{
		zy = t + 1.0;
		zc = cc;
	}
	else
	{
		zy = t;
		zc = 0.0;
	}

	/* Then get t as the next lag-1220 CSWB value */
	if (indx < 1220)
		t = Q[indx++];
	else
	{ /* refill Q[n] via Q[n-1220]-Q[n-1190]-c, */
		for (i = 0; i < 1220; i++)
		{
			j = (i < 30) ? i + 1190 : i - 30;
			t = Q[j] - Q[i] + c; /* Get next CSWB element */
			if (t > 0)
			{
				t = t - cc;
				c = cc;
			}
			else
			{
				t = t - cc + 1.0;
				c = 0.0;
			}
			Q[i] = t;
		}	 /* end i loop */
		indx = 1;
		t = Q[0]; /* set indx, exit 'else' with t=Q[0] */
	} /* end else segment; return t-zy mod 1 */
	
	return ((t < zy) ? 1.0 + (t - zy) : t - zy);	
	#pragma endregion
}


void MeshingGenie_mps_nd::generate_bg_grid()
{
	#pragma region Generate Background Grid:
	_num_ghost_layers = (size_t) ceil(sqrt(double(_ndim)));
	_s = _r / sqrt(double(_ndim)); _ss = _s;
	size_t num_bp = _boundary_points->size();
	for (size_t idim = 0; idim < _ndim; idim++)
	{
		_xmin[idim] = (*_boundary_points)[0][idim];
		_xmax[idim] = (*_boundary_points)[0][idim];
	}
	for (size_t ipnt = 1; ipnt < num_bp; ipnt++)
	{
		for (size_t idim = 0; idim < _ndim; idim++)
		{
			if  (_xmin[idim] > (*_boundary_points)[ipnt][idim]) _xmin[idim] = (*_boundary_points)[ipnt][idim];
			if  (_xmax[idim] < (*_boundary_points)[ipnt][idim]) _xmax[idim] = (*_boundary_points)[ipnt][idim];
		}
	}
	for (size_t idim = 0; idim < _ndim; idim++) _xmin[idim] -= _num_ghost_layers * _s; // taking into account ghost layers
	for (size_t idim = 0; idim < _ndim; idim++) _xmax[idim] += (_num_ghost_layers) * _s; // taking into account ghost layers
	for (size_t idim = 0; idim < _ndim; idim++) _nc[idim] = (size_t) ceil((_xmax[idim] - _xmin[idim]) / _s);
	#pragma endregion
}


void MeshingGenie_mps_nd::create_generic_neighbor_list()
{
	#pragma region Create a generic neighbor list sorted based on distance beween cells centers:
	size_t* jcent = new size_t[_ndim];
	size_t* j = new size_t[_ndim];
	int* dj = new int[_ndim];

	for (size_t idim = 0; idim < _ndim; idim++) jcent[idim] = _num_ghost_layers;
	for (size_t idim = 0; idim < _ndim; idim++) j[idim] = 0;

	size_t k_dim(_ndim - 1), jmax = 2 * _num_ghost_layers;

	_num_neighbor_cells = 0;
	while (true)
	{
		#pragma region Count Neighbor Cells:
		while (j[k_dim] <= jmax)
		{
			// get furthest corner of that cell from xcent
			if (conflicting_cells(j, jcent) && !cells_equal(j, jcent))  _num_neighbor_cells++;
			j[k_dim]++;
		}
		size_t kk_dim(k_dim - 1);
		bool done(false);
		while (true)
		{
			j[kk_dim]++;
			if (j[kk_dim] > jmax)
			{
				j[kk_dim] = 0;
				if (kk_dim == 0)
				{
					done = true;
					break;
				}
				kk_dim--;
			}
			else break;
		}
		if (done) break;
		j[k_dim] = 0;
		#pragma endregion
	}
	_generic_neighbor_cells = new int*[_num_neighbor_cells];

	k_dim = _ndim - 1; _num_neighbor_cells = 0;
	for (size_t idim = 0; idim < _ndim; idim++) j[idim] = 0;
	while (true)
	{
		#pragma region Fill-in Neighbor list:
		while (j[k_dim] <= jmax)
		{
			// get furthest corner of that cell from xcent
			if (conflicting_cells(j, jcent) && !cells_equal(j, jcent))  
			{
				for (size_t idim = 0; idim < _ndim; idim++) dj[idim] = int(j[idim]) - int(jcent[idim]);
				_generic_neighbor_cells[_num_neighbor_cells] = dj; _num_neighbor_cells++;
				dj = new int[_ndim];
			}
			j[k_dim]++;
		}
		size_t kk_dim(k_dim - 1);
		bool done(false);
		while (true)
		{
			j[kk_dim]++;
			if (j[kk_dim] > jmax)
			{
				j[kk_dim] = 0;
				if (kk_dim == 0)
				{
					done = true;
					break;
				}
				kk_dim--;
			}
			else break;
		}
		if (done) break;
		j[k_dim] = 0;
		#pragma endregion
	}

	// Lazy quadratc sorting should not affect performance done only once!
	for (size_t i = 0; i < _num_neighbor_cells; i++)
	{
		for (size_t idim = 0; idim < _ndim; idim++) j[idim] = jcent[idim] + _generic_neighbor_cells[i][idim];
		double dst_i = get_cells_sq_distance(jcent, j);

		for (size_t k = i + 1; k < _num_neighbor_cells; k++)
		{
			for (size_t idim = 0; idim < _ndim; idim++) j[idim] = jcent[idim] + _generic_neighbor_cells[k][idim];
			double dst_k = get_cells_sq_distance(jcent, j);
			if (dst_k < dst_i)
			{
				int* tmp = _generic_neighbor_cells[i];
				_generic_neighbor_cells[i] = _generic_neighbor_cells[k];
				_generic_neighbor_cells[k] = tmp;
				dst_i = dst_k;
			}
		}
	}
	delete[] j;
	delete[] jcent;
	delete[] dj;
	#pragma endregion
}

void MeshingGenie_mps_nd::initiate_active_pool()
{
	#pragma region initiate Active pool with all cell in the background grid:
	size_t* icell = new size_t[_ndim];
	for (size_t idim = 0; idim< _ndim; idim++) icell[idim] = 0;

	size_t k_dim(_ndim - 1);
	while (true)
	{
		while (icell[k_dim] < _nc[k_dim])
		{
			size_t* jcell = new size_t[_ndim];
			for (size_t idim = 0; idim< _ndim; idim++) jcell[idim] = icell[idim];
			add_cell(jcell, _active_cells);
			icell[k_dim]++;
		}

		size_t kk_dim(k_dim - 1);
		bool done(false);
		while (true)
		{
			icell[kk_dim]++;
			if (icell[kk_dim] >= _nc[kk_dim])
			{
				icell[kk_dim] = 0;
				if (kk_dim == 0)
				{
					done = true;
					break;
				}
				kk_dim--;
			}
			else break;
		}
		if (done) break;
		icell[k_dim] = 0;
	}
	delete [] icell;
	#pragma endregion
}

void MeshingGenie_mps_nd::identify_boundary_cells()
{
	#pragma region Identifying Boundary Cells via Random Sampling of boundary faces:
	size_t num_bf(_boundary_faces->size());

	size_t* icell = new size_t[_ndim];
	size_t* iparent = new size_t[_ndim];

	double* u = new double[_ndim];
	double* x = new double[_ndim];
	double* xfmin = new double[_ndim];
	double* xfmax = new double[_ndim];

	for (size_t iface = 0; iface < num_bf; iface++)
	{
		// retrieve bounding box of that face
		for (size_t idim = 0; idim < _ndim; idim++)
		{
			size_t face_corner = (*_boundary_faces)[iface][0];
			xfmin[idim] = (*_boundary_points)[face_corner][idim];
			xfmax[idim] = (*_boundary_points)[face_corner][idim];
		}

		size_t num_face_corners = (*_boundary_faces)[iface].size();
		for (size_t ipnt = 1; ipnt < num_face_corners; ipnt++)
		{
			size_t face_corner = (*_boundary_faces)[iface][ipnt];
			for (size_t idim = 0; idim < _ndim; idim++)
			{
				if  (xfmin[idim] > (*_boundary_points)[face_corner][idim]) xfmin[idim] = (*_boundary_points)[face_corner][idim];
				if  (xfmax[idim] < (*_boundary_points)[face_corner][idim]) xfmax[idim] = (*_boundary_points)[face_corner][idim];
			}
		}

		// Estimating the number of cells that this face extends through
		size_t num_sample_points(100);
		for (size_t idim = 0; idim < _ndim; idim++)
		{
			size_t num_oneD_cells = (size_t) ceil((xfmax[idim]- xfmin[idim]) / _ss);
			num_oneD_cells++;
			num_sample_points *= num_oneD_cells;
		}
		
		for (size_t isample = 0; isample < num_sample_points; isample++)
		{
			// generate random barycentric point
			double u_norm(0.0);
			for (size_t idim = 0; idim < _ndim; idim++)
			{
				u[idim] = generate_a_random_number();
				u_norm += u[idim];
			}
			u_norm = 1.0 / u_norm;
			for (size_t idim = 0; idim < _ndim; idim++) u[idim] *= u_norm;
			
			for (size_t idim = 0; idim < _ndim; idim++) 
			{
				x[idim] = 0.0;

				size_t num_face_corners = (*_boundary_faces)[iface].size();
				for (size_t ipnt = 0; ipnt < num_face_corners; ipnt++)
				{
					size_t face_corner = (*_boundary_faces)[iface][ipnt];
					x[idim] += u[ipnt] * (*_boundary_points)[face_corner][idim];
				}
				// Retrieve Cell indices
				icell[idim] = size_t(floor((x[idim] - _xmin[idim]) / _ss));
				iparent[idim] = size_t(floor((x[idim] - _xmin[idim]) / _s));
			}

			if (!cell_in_vec(iparent, _covered_cells))
			{
				// parent cell is not covered so add this cell to _active_boundary_cells 
				if (add_boundary_cell(icell)) 
				{
					icell = new size_t[_ndim];
				}
			}
		}
	}
	delete[] icell;	delete[] iparent;
	delete[] u;	delete[] x;
	delete[] xfmin;	delete[] xfmax;
	#pragma endregion
}

void MeshingGenie_mps_nd::identify_exterior_and_interior_cells()
{
	#pragma region Identify Exterior Cells Assuming One Material only, may contain holes:
	size_t num_face_neighbors(2 * _ndim);
	
	// Assume that _r is smaller than feature size, nonconvext domains with holes
	if (fabs(_ss / _s - 1.0) < 1E-10)
	{
		// Level zero	
		size_t* icell = new size_t[_ndim];
		for (size_t idim = 0; idim< _ndim; idim++) icell[idim] = 0;
		add_cell(icell, _exterior_cells);
	}

	size_t num_active_cells(_active_cells.size());
	size_t* jcell = new size_t[_ndim];
	while (true)
	{
		#pragma region propagating exterior cells thorugh active cells until we hit a boundary cell:
		bool done(true);

		for (size_t iactive = 0; iactive < num_active_cells; iactive++)
		{
			if (!cell_in_vec(_active_cells[iactive], _exterior_cells)) continue;
			
			// _active_cells[iactive] is an exterior cell, try to propagate that to its face neighbors
			for (size_t i_neighbor = 0; i_neighbor < num_face_neighbors; i_neighbor++)
			{
				bool valid_cell(true);
				for (size_t idim = 0; idim < _ndim; idim++)
				{
					if ((int)_active_cells[iactive][idim] + _generic_neighbor_cells[i_neighbor][idim] < 0) {valid_cell = false; break;}
					if ((int)_active_cells[iactive][idim] + _generic_neighbor_cells[i_neighbor][idim] >= (int)_nc[idim]) {valid_cell = false; break;}

					jcell[idim] = _active_cells[iactive][idim] + _generic_neighbor_cells[i_neighbor][idim];
				}
				if (!valid_cell) continue;

				if (cell_in_vec(jcell, _boundary_cells)) continue;

				if (cell_in_vec(jcell, _exterior_cells)) continue;

				add_cell(jcell, _exterior_cells);
				jcell = new size_t[_ndim];
				done = false;
			}
		}
		if (done) break;
		#pragma endregion
	}

	#pragma region Propagating Exterior cells to interior cells through active boundary cells:
	size_t num_bcells(_boundary_cells.size());
	for (size_t ibactive = 0; ibactive < num_bcells; ibactive++)
	{
		bool touching_exterior(false);
		for (size_t i_neighbor = 0; i_neighbor < num_face_neighbors; i_neighbor++)
		{
			bool valid_cell(true);
			for (size_t idim = 0; idim < _ndim; idim++)
			{
				if ((int)_boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim] < 0) {valid_cell = false; break;}
				if ((int)_boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim] >= (int)_nc[idim]) {valid_cell = false; break;}

				jcell[idim] = _boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim];
			}
			if (!valid_cell) continue;
			if (cell_in_vec(jcell, _exterior_cells)) {touching_exterior = true; break;} // A boundary cell touching the exterior
		}
		if (touching_exterior)
		{
			for (size_t i_neighbor = 0; i_neighbor < num_face_neighbors; i_neighbor++)
			{
				bool valid_cell(true);
				for (size_t idim = 0; idim < _ndim; idim++)
				{
					if ((int)_boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim] < 0) {valid_cell = false; break;}
					if ((int)_boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim] >= (int)_nc[idim]) {valid_cell = false; break;}

					jcell[idim] = _boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim];
				}
				if (!valid_cell) continue;
				if (!cell_in_vec(jcell, _exterior_cells) && !cell_in_vec(jcell, _boundary_cells) && !cell_in_vec(jcell, _interior_cells))
				{ 
					// jcell is a non-classified active cell, classify it as an interior cell
					add_cell(jcell, _interior_cells);
					jcell = new size_t[_ndim];
				}
			}
		}
	}
	#pragma endregion

	while (true)
	{
		#pragma region propagating interior cells thorugh active cells until we hit a boundary cell:
		bool done(true);

		for (size_t iactive = 0; iactive < num_active_cells; iactive++)
		{
			if (!cell_in_vec(_active_cells[iactive], _interior_cells)) continue;
			
			// _active_cells[iactive] is an exterior cell, try to propagate that to its face neighbors
			for (size_t i_neighbor = 0; i_neighbor < num_face_neighbors; i_neighbor++)
			{
				bool valid_cell(true);
				for (size_t idim = 0; idim < _ndim; idim++)
				{
					if ((int)_active_cells[iactive][idim] + _generic_neighbor_cells[i_neighbor][idim] < 0) {valid_cell = false; break;}
					if ((int)_active_cells[iactive][idim] + _generic_neighbor_cells[i_neighbor][idim] >= (int)_nc[idim]) {valid_cell = false; break;}

					jcell[idim] = _active_cells[iactive][idim] + _generic_neighbor_cells[i_neighbor][idim];
				}
				if (!valid_cell) continue;

				if (cell_in_vec(jcell, _boundary_cells)) continue;

				if (cell_in_vec(jcell, _interior_cells)) continue;

				add_cell(jcell, _interior_cells);
				jcell = new size_t[_ndim];
				done = false;
			}
		}
		if (done) break;
		#pragma endregion
	}

	#pragma region Propagating Interior cells to Exterior cells through active boundary cells for domain with holes:
	for (size_t ibactive = 0; ibactive < num_bcells; ibactive++)
	{
		bool touching_interior(false);
		for (size_t i_neighbor = 0; i_neighbor < num_face_neighbors; i_neighbor++)
		{
			bool valid_cell(true);
			for (size_t idim = 0; idim < _ndim; idim++)
			{
				if ((int)_boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim] < 0) {valid_cell = false; break;}
				if ((int)_boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim] >= (int)_nc[idim]) {valid_cell = false; break;}

				jcell[idim] = _boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim];
			}
			if (!valid_cell) continue;
			if (cell_in_vec(jcell, _interior_cells)) {touching_interior = true; break;} // A boundary cell touching the exterior
		}
		if (touching_interior)
		{
			for (size_t i_neighbor = 0; i_neighbor < num_face_neighbors; i_neighbor++)
			{
				bool valid_cell(true);
				for (size_t idim = 0; idim < _ndim; idim++)
				{
					if ((int)_boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim] < 0) {valid_cell = false; break;}
					if ((int)_boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim] >= (int)_nc[idim]) {valid_cell = false; break;}

					jcell[idim] = _boundary_cells[ibactive][idim] + _generic_neighbor_cells[i_neighbor][idim];
				}
				if (!valid_cell) continue;
				if (!cell_in_vec(jcell, _exterior_cells) && !cell_in_vec(jcell, _boundary_cells) && !cell_in_vec(jcell, _interior_cells))
				{ 
					// jcell is a non-classified active cell, classify it as an interior cell
					add_cell(jcell, _exterior_cells);
					jcell = new size_t[_ndim];
				}
			}
		}
	}
	#pragma endregion

	while (true)
	{
		#pragma region propagating exterior cells thorugh active cells until we hit a boundary cell:
		bool done(true);

		for (size_t iactive = 0; iactive < num_active_cells; iactive++)
		{
			if (!cell_in_vec(_active_cells[iactive], _exterior_cells)) continue;
			
			// _active_cells[iactive] is an exterior cell, try to propagate that to its face neighbors
			for (size_t i_neighbor = 0; i_neighbor < num_face_neighbors; i_neighbor++)
			{
				bool valid_cell(true);
				for (size_t idim = 0; idim < _ndim; idim++)
				{
					if ((int)_active_cells[iactive][idim] + _generic_neighbor_cells[i_neighbor][idim] < 0) {valid_cell = false; break;}
					if ((int)_active_cells[iactive][idim] + _generic_neighbor_cells[i_neighbor][idim] >= (int)_nc[idim]) {valid_cell = false; break;}

					jcell[idim] = _active_cells[iactive][idim] + _generic_neighbor_cells[i_neighbor][idim];
				}
				if (!valid_cell) continue;

				if (cell_in_vec(jcell, _boundary_cells)) continue;

				if (cell_in_vec(jcell, _exterior_cells)) continue;

				add_cell(jcell, _exterior_cells);
				jcell = new size_t[_ndim];
				done = false;
			}
		}
		if (done) break;
		#pragma endregion
	}

	delete[] jcell;
	#pragma endregion
}

void MeshingGenie_mps_nd::update_active_pool(size_t refLevel)
{
	#pragma region update active pool - Still need to consider cells from boundary cells:
	if (refLevel == 0)
	{
		size_t num_active_cells = _active_cells.size();
		for (size_t icell = 0; icell < num_active_cells; icell++)
		{
			delete[] _active_cells[icell];
		}
		_active_cells.clear();
		num_active_cells = _interior_cells.size();
		_active_cells.resize(num_active_cells);
		for (size_t icell = 0; icell < num_active_cells; icell++)
		{
			_active_cells[icell] = _interior_cells[icell];
		}
	}
	else
	{
		// refine active cells, add noncovered
		std::vector<size_t*> active_children;		
		for (size_t icell = 0; icell < _num_active_cells; icell++)
		{
			get_uncovered_children(_active_cells[icell], refLevel, active_children);
		}
		if (refLevel > 1)
		{
			// destroy parents of the previous level
			size_t num_active_cells = _active_cells.size();
			for (size_t icell = 0; icell < num_active_cells; icell++)
			{
				delete[] _active_cells[icell];
			}
		}
		_active_cells.clear();
		size_t num_active_children = active_children.size();
		_active_cells.resize(num_active_children);
		for (size_t icell = 0; icell < num_active_children; icell++)
		{
			_active_cells[icell] = active_children[icell];
		}
		_ss *= 0.5;
	}
	//_plotter.plot_active_pool(this);
	#pragma endregion
}

void MeshingGenie_mps_nd::throw_darts(size_t refLevel)
{
	#pragma region Throw darts:
	update_active_pool(refLevel);

	size_t* iparent = new size_t[_ndim];
	double* x = new double[_ndim];
	_num_active_cells = _active_cells.size();
	size_t num_darts(_num_active_cells);
	for (size_t idart = 0; idart < num_darts; idart++)
	{
		// pick a random cell
		if (_num_active_cells == 0) break;
		double u = generate_a_random_number();
		size_t icell = size_t(u * _num_active_cells);
		if (icell == _num_active_cells) icell--; // since u could be 1

		// get parent cell
		for (size_t idim = 0; idim < _ndim; idim++) iparent[idim] = _active_cells[icell][idim] / ipow(2, refLevel);

		if (cell_in_vec(iparent, _covered_cells))
		{
			// deactivate icell cause its parent is covered
			_num_active_cells--;
			size_t* tmp = _active_cells[icell];
			_active_cells[icell] = _active_cells[_num_active_cells];
			_active_cells[_num_active_cells] = tmp;
			continue;
		}

		_num_darts++;

		// pick a random point from icell
		for (size_t idim = 0; idim < _ndim; idim++)
		{
			u = generate_a_random_number();
			x[idim] = _xmin[idim] + (_active_cells[icell][idim] + u) * _ss;
		}

		_current_dart = x;
		//_plotter.plot_active_pool(this);

		if (valid_dart(x, iparent))
		{
			if (cell_in_vec(iparent, _boundary_cells))
			{
				size_t index = find_cell_binary(iparent, _boundary_cells);
				_boundary_cell_points[index] = x; _num_inserted_points++;
				x = new double[_ndim];
			}
			else if (cell_in_vec(iparent, _interior_cells))
			{
				size_t index = find_cell_binary(iparent, _interior_cells);
				_interior_cell_points[index] = x; _num_inserted_points++;
				x = new double[_ndim];
			}

			// deactivate icell and cover its parent
			_num_active_cells--;
			size_t* tmp = _active_cells[icell];
			_active_cells[icell] = _active_cells[_num_active_cells];
			_active_cells[_num_active_cells] = tmp;
			add_cell(iparent, _covered_cells);
			iparent = new size_t[_ndim];
			continue;
		}
	}

	_current_dart = 0;
	delete[] iparent;
	delete[] x;
	#pragma endregion
}

inline void MeshingGenie_mps_nd::get_uncovered_children(size_t* active_cell, size_t refLevel, std::vector<size_t*> &children)
{
	#pragma region Get uncovered children:
	size_t* icell = new size_t[_ndim];
	for (size_t idim = 0; idim< _ndim; idim++) icell[idim] = 0;

	size_t* jcell = new size_t[_ndim];

	size_t k_dim(_ndim - 1);
	while (true)
	{
		while (icell[k_dim] < 2)
		{
			for (size_t idim = 0; idim< _ndim; idim++) jcell[idim] = active_cell[idim] * 2 + icell[idim];

			if (!covered_cell(jcell, refLevel))
			{
				children.push_back(jcell);
				jcell = new size_t[_ndim];
			}
			icell[k_dim]++;
		}

		size_t kk_dim(k_dim - 1);
		bool done(false);
		while (true)
		{
			icell[kk_dim]++;
			if (icell[kk_dim] == 2)
			{
				icell[kk_dim] = 0;
				if (kk_dim == 0)
				{
					done = true;
					break;
				}
				kk_dim--;
			}
			else break;
		}
		if (done) break;
		icell[k_dim] = 0;
	}
	delete [] icell;
	delete [] jcell;
    #pragma endregion
}


inline bool MeshingGenie_mps_nd::valid_dart(double* dart, size_t* dart_parent_cell)
{
	#pragma region Check if dart is Valid:
	size_t* parent_neighbor = new size_t[_ndim];
	
	for (size_t i_neighbor = 0; i_neighbor < _num_neighbor_cells; i_neighbor++)
	{
		bool valid_cell(true);
		for (size_t idim = 0; idim < _ndim; idim++)
		{
			if ((int)dart_parent_cell[idim] + _generic_neighbor_cells[i_neighbor][idim] < 0) {valid_cell = false; break;}

			if ((int)dart_parent_cell[idim] + _generic_neighbor_cells[i_neighbor][idim] >= (int)_nc[idim]) {valid_cell = false; break;}

			parent_neighbor[idim] = dart_parent_cell[idim] + _generic_neighbor_cells[i_neighbor][idim];
		}
		if (!valid_cell) continue;

		double* x = 0;
		if (cell_in_vec(parent_neighbor, _boundary_cells))
		{
			size_t index = find_cell_binary(parent_neighbor, _boundary_cells);
			x = _boundary_cell_points[index];
		}
		else if (cell_in_vec(parent_neighbor, _interior_cells))
		{
			size_t index = find_cell_binary(parent_neighbor, _interior_cells);
			x = _interior_cell_points[index];
		}
		if (x != 0)
		{
			// check if disk at x covers icell
			double dd(0.0);
			for (size_t idim = 0; idim < _ndim; idim++)
			{
				double dx = dart[idim] - x[idim];
				dd+= dx * dx;
			}
			if (dd < _rsq) 
			{
				delete[] parent_neighbor;
				return false;
			}
		}
	}
	delete[] parent_neighbor;
	return true;
	#pragma endregion
}

inline bool MeshingGenie_mps_nd::covered_cell(size_t* icell, size_t refLevel)
{
	#pragma region Check if a cell is covered:
	double ss = _s / ipow(2, refLevel);
	
	size_t* iparent = new size_t[_ndim];
	for (size_t idim = 0; idim < _ndim; idim++) iparent[idim] = icell[idim] / ipow(2, refLevel);

	if (cell_in_vec(iparent, _covered_cells))
	{
		delete[] iparent;
		return true;
	}

	size_t* parent_neighbor = new size_t[_ndim];
	for (size_t i_neighbor = 0; i_neighbor < _num_neighbor_cells; i_neighbor++)
	{
		bool valid_cell(true);
		for (size_t idim = 0; idim < _ndim; idim++)
		{
			if ((int)iparent[idim] + _generic_neighbor_cells[i_neighbor][idim] < 0) {valid_cell = false; break;}

			if ((int)iparent[idim] + _generic_neighbor_cells[i_neighbor][idim] >= (int)_nc[idim]) {valid_cell = false; break;}

			parent_neighbor[idim] = iparent[idim] + _generic_neighbor_cells[i_neighbor][idim];
		}
		if (!valid_cell) continue;

		double* x = 0;
		if (cell_in_vec(parent_neighbor, _boundary_cells))
		{
			size_t index = find_cell_binary(parent_neighbor, _boundary_cells);
			x = _boundary_cell_points[index];
		}
		else if (cell_in_vec(parent_neighbor, _interior_cells))
		{
			size_t index = find_cell_binary(parent_neighbor, _interior_cells);
			x = _interior_cell_points[index];
		}
		if (x != 0)
		{
			// check if disk at x covers icell
			double dd(0.0);
			for (size_t idim = 0; idim < _ndim; idim++)
			{
				double xx = _xmin[idim] + icell[idim] * ss;
				if (fabs(x[idim] - xx) < fabs(x[idim] - xx - ss)) xx += ss;
				double dx = xx - x[idim];
				dd += dx * dx;
			}
			if (dd < _rsq) 
			{
				delete[] iparent;
				delete[] parent_neighbor;
				return true;
			}
		}
	}
	delete[] iparent;
	delete[] parent_neighbor;
	return false;
	#pragma endregion
}

inline bool MeshingGenie_mps_nd::cell_in_vec(size_t* icell, std::vector<size_t*> &cell_vec)
{
	#pragma region Use binary search to check if icell exists in cell_vec:
	vec_range loc = check_location(icell, cell_vec);
	if (loc == last) return true;
	if (loc == within_range)
	{
		size_t index = find_cell_binary(icell, cell_vec);
		if (cells_equal(icell, cell_vec[index])) return true; // icell exists in cell_vec
	}
	return false;
	#pragma endregion
}

inline bool MeshingGenie_mps_nd::add_cell(size_t* icell, std::vector<size_t*> &cell_vec)
{
	#pragma region Use binary search to add icell to cell_vec:
	vec_range loc = check_location(icell, cell_vec);
	if (loc == empty || loc == back)
	{
		// active_boundary_cell is empty
		cell_vec.push_back(icell);
		return true;
	}
	if (loc == front)
	{
		// if active_boundary_cell < _active_boundary_cells[0]
		insert_cell(icell, cell_vec, 0);
		return true;
	}

	if (loc == last)
	{
		// active_boundary_cell == _active_boundary_cells[num_active_bcells - 1]
		return false;
	}

	size_t index = find_cell_binary(icell, cell_vec);       

	if (cells_equal(icell, cell_vec[index])) return false; // active_boundary_cell exists in _active_boundary_cells
	
	// insert cell at index + 1
	insert_cell(icell, cell_vec, index + 1);
	return true;
	#pragma endregion
}

inline bool MeshingGenie_mps_nd::add_boundary_cell(size_t* boundary_cell)
{
	return add_cell(boundary_cell, _boundary_cells);
}


inline MeshingGenie_mps_nd::vec_range MeshingGenie_mps_nd::check_location(size_t* icell,std::vector<size_t*> &cell_vec)
{
	#pragma region Check the location of icell with regard to the range of cell_vec
	size_t num_cells(cell_vec.size());
	if (num_cells == 0)                                   return empty;   // cell_vec is empty
	if (cells_less_than(icell, cell_vec[0]))              return front;   // icell < cell_vec[0]
	if (cells_less_than( cell_vec[num_cells - 1], icell)) return back;    // icell > cell_vec[num_cells - 1]
	if (cells_equal(cell_vec[num_cells - 1], icell))      return last;    // icell == cell_vec[num_cells - 1]
	return within_range;                                                  // cell_vec[0] < icell < cell_vec[num_cells - 1]
	#pragma endregion
}

inline size_t MeshingGenie_mps_nd::find_cell_binary(size_t* icell,std::vector<size_t*> &cell_vec)
{
	#pragma region Find an index such that cell_vec[ist = index] <= icell < cell_vec[index + 1]
	size_t num_cells(cell_vec.size());
	size_t ist(0), iend(num_cells);
	while (iend - ist > 1)
	{
		size_t imid = (ist + iend) / 2;		
		if (cells_less_than(icell, cell_vec[imid]))   
		{
			//  ist <= active_boundary_cell < imid
			iend = imid;
		}
		else
		{
			//  imid <= active_boundary_cell < iend
			ist = imid;                                                          
		}
	}
	return ist;
	#pragma endregion
}

inline void MeshingGenie_mps_nd::insert_cell(size_t* icell,std::vector<size_t*> &cell_vec, size_t index)
{
	#pragma region Insert icell at location index in Vector cell_vec:
	size_t num_cells(cell_vec.size());
	cell_vec.push_back(icell);
	for (size_t i = num_cells; i > index; i--)
	{
		cell_vec[i] = cell_vec[i-1];
	}
	cell_vec[index] = icell;
	#pragma endregion
}


inline bool MeshingGenie_mps_nd::cells_less_than(size_t* icell, size_t* jcell)
{
	#pragma region Check if icell is less than jcell:
	for (size_t i = 0; i < _ndim; i++)
	{
		if (icell[i] < jcell[i]) return true;       // icell < jcell
		else if (icell[i] > jcell[i]) return false; // icell > jcell
	}
	return false; // icell == jcell
	#pragma endregion
}

inline bool MeshingGenie_mps_nd::cells_equal(size_t* icell, size_t* jcell)
{
	#pragma region Check if icell equals jcell:
	for (size_t i = 0; i < _ndim; i++)
	{
		if (icell[i] != jcell[i]) return false;
	}
	return true;
	#pragma endregion
}

inline bool MeshingGenie_mps_nd::conflicting_cells(size_t* icell, size_t* jcell)
{
	#pragma region Check if there is a chance that disks in the two cells may intersect:
	double hh(0.0); // closest distance between the two cells
	for (size_t idim = 0; idim < _ndim; idim++)
	{
		if (abs(int(icell[idim]) - int(jcell[idim])) > int(_num_ghost_layers)) return false;
		
		double xi = _xmin[idim] + icell[idim] * _s;
		double xj = _xmin[idim] + jcell[idim] * _s;
		double dx;
		if (icell[idim] > jcell[idim]) dx = xi - (xj + _s);
		else if (icell[idim] < jcell[idim]) dx = xj - (xi + _s);
		else dx = 0.0;
		hh += dx * dx;
	}
	if (hh > _rsq - 1E-10) return false;
	return true;
	#pragma endregion
}

inline double MeshingGenie_mps_nd::get_cells_sq_distance(size_t* icell, size_t* jcell)
{
	#pragma region Get distance squared between the centers of two cells:
	double hh(0.0); // closest distance between the two cells
	for (size_t idim = 0; idim < _ndim; idim++)
	{
		double dx = (int(icell[idim]) - int(jcell[idim])) * _s;
		hh += dx * dx;
	}
	return hh;
	#pragma endregion
}

inline size_t MeshingGenie_mps_nd::ipow(size_t i, size_t p)
{
	size_t ans(1);
	for (size_t j = 0; j < p; j++) ans*= i;
	return ans;
}

