
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

 * Last modified: 10/29/2012
********************************************************************************/


#ifndef MESHING_GENIE_mps_nd_H
#define MESHING_GENIE_mps_nd_H

#include "MeshingGenie_defs.h"

class MeshingGenie_mps_nd   
{	
    public:
        //! constructor
		MeshingGenie_mps_nd(){ };
   

        //! Destructor                
        ~MeshingGenie_mps_nd(){ };


		int solve_mps(size_t ndim, double r, size_t random_seed, 
			          std::vector< std::vector<double> > &boundary_points,
			          std::vector< std::vector<size_t> > &boundary_faces, 
				      std::vector< std::vector<double> > &sample_points);
	private:
		enum vec_range{empty,                     // vector is empty
			           front,                     // entity < vec[first]
			           back,                      // entity > vec[last]
					   last,                      // entity = vec[last]
					   within_range,              // vec[first] <= entity < vec[last]
					  };

	private:
		void init();
		void cleanup();

		void initiate_random_generator(unsigned long x);
		double generate_a_random_number();

		void generate_bg_grid();
		void create_generic_neighbor_list();
		

		// Identifying Active Boundary Cells via Random Sampling of boundary faces
		void initiate_active_pool();
		void identify_boundary_cells();
		void identify_exterior_and_interior_cells();
		void update_active_pool(size_t refLevel);
		void throw_darts(size_t refLevel);

		// retrieve uncovered children of an active cell
		inline void get_uncovered_children(size_t* icell, size_t refLevel, std::vector<size_t*> &children);

		inline bool valid_dart(double* dart, size_t* dart_parent_cell);
		inline bool covered_cell(size_t icell, size_t refLevel);

		// use binary search to check if icell exists in cell_vec
		inline bool cell_in_vec(size_t* icell, std::vector<size_t*> &cell_vec);
		
		// use binary search to add icell to cell_vec, false == icell already exists in cell_vec
		inline bool add_cell(size_t* icell, std::vector<size_t*> &cell_vec);
		inline bool add_boundary_cell(size_t* active_boundary_cell);

		// Check the location of icell with regard to the range of cell_vec
		inline vec_range check_location(size_t* icell,std::vector<size_t*> &cell_vec);
		// returns a number index from 0 to cell_vec.siz(0) - 1 s.t. cell_vec[index] <= icell < cell_vec[index+1]
		inline size_t find_cell_binary(size_t* icell,std::vector<size_t*> &cell_vec);
		// insert icell in the location at index
		inline void insert_cell(size_t* icell,std::vector<size_t*> &cell_vec, size_t index);

		// Check if icell < jcell
		inline bool cells_less_than(size_t* icell, size_t* jcell);
		// Check if icell == jcell
		inline bool cells_equal(size_t* icell, size_t* jcell);
		// Check if there is a chance that disks in the two cells may intersect:
		inline bool conflicting_cells(size_t* icell, size_t* jcell);
		// Get distance squared between the centers of two cells:
		inline double get_cells_sq_distance(size_t* icell, size_t* jcell);
		
		// returns the i^p
		inline size_t ipow(size_t i, size_t p);

    private:
		

		size_t _ndim;                                         // number of dimensions
		double _r, _rsq;                                      // distribution radius
		size_t _random_seed;								  // random seed
		std::vector< std::vector<double> > *_boundary_points; // pointer to boundary points
		std::vector< std::vector<size_t> > *_boundary_faces;  // pointer to boundary faces
		std::vector< std::vector<double> > *_sample_points;   // pointer to sample points	

		size_t _expected_number_of_points;

		// variables for Random number generator
		double Q[1220];
		int indx;
		double cc;
		double c; /* current CSWB */
		double zc;	/* current SWB `borrow` */
		double zx;	/* SWB seed1 */
		double zy;	/* SWB seed2 */
		size_t qlen;/* length of Q array */	


		// background grids
		size_t _num_ghost_layers; 
		size_t _num_valid_cells, _num_inserted_points;

		double _s;      // parent grid spacing
		double _ss;     // active grid spacing
		size_t* _nc;    // Number of cells in each direction
		double* _xmin;  // lower left corner of bg grid
		double* _xmax;  // upper right corner of bg grid

		size_t _num_neighbor_cells;                // number of neighbor cells
		int**  _generic_neighbor_cells;            // list of generic neigbor cells

		
		std::vector< size_t* > _active_cells;               // active Cells


		std::vector< size_t* > _boundary_cells;             // boundary cells
		std::vector< size_t* > _exterior_cells;             // cells that lie outside the domain
		std::vector< size_t* > _interior_cells;             // cells that lie inside the domain
		std::vector< size_t* > _covered_cells;              //  cells that are already covered
		
		std::vector< double* > _boundary_cell_points;       // sample points in boundary cells
		std::vector< double* > _interior_cell_points;       // sample points in interior cells		
		
};
                                
#endif	



