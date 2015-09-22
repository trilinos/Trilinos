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
      This class helps debugging the method presented in the following article
      M. S. Ebeida, S. A. Mitchell, A. Patney, A. A. Davidson, and J. D. Owens, 
      "A simple algorithm for maximal Poisson-disk sampling in high dimensions", 
	  Computer Graphics Forum (Eurographics 2012), 31(2), May 2012.
   
   input: mesh data (Active pool, boundary cells, ... etc)

   output: a postscript file that shows various details in 2D domains

 * Last modified: 11/21/2012
********************************************************************************/

#include "MeshingGenie_plotter_2d.h"

int MeshingGenie_plotter_2d::plot_active_pool(MeshingGenie_mps_nd* mps_solver)
{
    #pragma region Plot Background grid:
	std::fstream file("results.ps", std::ios::out);
    file << "%!PS-Adobe-3.0" << std::endl;
    file << "72 72 scale     % one unit = one inch" << std::endl;

    #pragma region Retrieve bounding box, scale, and translate:
    double xmin(mps_solver->_xmin[0]);
	double ymin(mps_solver->_xmin[1]);
	double Lx(mps_solver->_xmax[0] - mps_solver->_xmin[0]);
	double Ly(mps_solver->_xmax[1] - mps_solver->_xmin[1]);        
    
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

	file << "/seg_blue      % stack: x1 y1 x2 y2" << std::endl;
    file << "{newpath" << std::endl; 
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
	file << " gsave" << std::endl;       // saves saved color and line width
	file << "0 0 1 setrgbcolor" << std::endl; 
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
	file << " grestore" << std::endl;    // stores saved color and line width
    file << "} def" << std::endl;

    file << "/seg_bold      % stack: x1 y1 x2 y2" << std::endl;
    file << "{newpath" << std::endl; 
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
    file << "} def" << std::endl;

    file << "/seg      % stack: x1 y1 x2 y2" << std::endl;
    file << "{newpath" << std::endl; 
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " 0.008 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
    file << "} def" << std::endl;

    file << "/quad      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
    file << "{newpath" << std::endl; 
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " 0.0 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
    file << "} def" << std::endl;

	file << "/quad_bold      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
    file << "{newpath" << std::endl; 
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
    file << "} def" << std::endl;

	file << "/quad_dark      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
    file << "{newpath" << std::endl; 
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
	file << " gsave" << std::endl;
    file << " 0.2 setgray fill" << std::endl;
	file << " grestore" << std::endl;
    file << " 0.002 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
    file << "} def" << std::endl;

	file << "/quad_light      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
    file << "{newpath" << std::endl; 
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
	file << " gsave" << std::endl;
    file << " 0.9 setgray fill" << std::endl;
	file << " grestore" << std::endl;
    file << " 0.002 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
    file << "} def" << std::endl;

	file << "/quad_white      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
    file << "{newpath" << std::endl; 
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
	file << " gsave" << std::endl;
    file << " 1.0 setgray fill" << std::endl;
	file << " grestore" << std::endl;   
    file << "} def" << std::endl;
    #pragma endregion
    
    double s(mps_solver->_r * 0.05);
   
	// plot filled circles
	size_t num_interior_cells(mps_solver->_interior_cells.size());
	for (size_t icell = 0; icell < num_interior_cells; icell++)
	{
		if (mps_solver->_interior_cell_points[icell] == 0) continue;
		double* x = mps_solver->_interior_cell_points[icell];
		file << x[0] * scale << "  " << x[1] * scale << "  " << mps_solver->_r * scale << "  ";
		file << "fcirc"     << std::endl;			
	}
	
	// plot active cells
	size_t num_active_cells = mps_solver->_active_cells.size();
	for (size_t icell = 0; icell < num_active_cells; icell++)
	{
		size_t* i = mps_solver->_active_cells[icell];

		double xo = mps_solver->_xmin[0] + i[0] * mps_solver->_ss;
		double yo = mps_solver->_xmin[1] + i[1] * mps_solver->_ss;
		double xn = xo + mps_solver->_ss;
		double yn = yo + mps_solver->_ss;

		// plot cell
		file << xo * scale << "  " << yo * scale << "  ";
		file << xn * scale << "  " << yo * scale << "  ";          
		file << xn * scale << "  " << yn * scale << "  ";          
		file << xo * scale << "  " << yn * scale << "  "; 
		file << "quad_light"     << std::endl;			
	}

	size_t num_bcells(mps_solver->_boundary_cells.size());
	for (size_t icell = 0; icell < num_bcells; icell++)
	{
		size_t* i = mps_solver->_boundary_cells[icell];

		double xo = mps_solver->_xmin[0] + i[0] * mps_solver->_s;
		double yo = mps_solver->_xmin[1] + i[1] * mps_solver->_s;
		double xn = xo + mps_solver->_s;
		double yn = yo + mps_solver->_s;

		// plot cell
		file << xo * scale << "  " << yo * scale << "  ";
		file << xn * scale << "  " << yo * scale << "  ";          
		file << xn * scale << "  " << yn * scale << "  ";          
		file << xo * scale << "  " << yn * scale << "  "; 
		file << "quad_dark"     << std::endl;			
	}


	// plot discs boundaries
	for (size_t icell = 0; icell < num_interior_cells; icell++)
	{
		if (mps_solver->_interior_cell_points[icell] == 0) continue;
		double* x = mps_solver->_interior_cell_points[icell];
		file << x[0] * scale << "  " << x[1] * scale << "  " << mps_solver->_r * scale << "  ";
		file << "circ"     << std::endl;			
	}

	for (size_t icell = 0; icell < num_interior_cells; icell++)
	{
		if (mps_solver->_interior_cell_points[icell] == 0) continue;
		double* x = mps_solver->_interior_cell_points[icell];

		// plot vertex
        file << (x[0] - s) * scale << "  " << (x[1] - s) * scale << "  ";
        file << (x[0] + s) * scale << "  " << (x[1] + s) * scale << "  ";
        file << "seg_bold"     << std::endl;		

		file << (x[0] + s) * scale << "  " << (x[1] - s) * scale << "  ";
        file << (x[0] - s) * scale << "  " << (x[1] + s) * scale << "  ";
        file << "seg_bold"     << std::endl;    
	}

	// plotting current dart
	if (mps_solver->_current_dart != 0)
	{
		double* x = mps_solver->_current_dart;
		// plot vertex
        file << (x[0] - s) * scale << "  " << (x[1] - s) * scale << "  ";
        file << (x[0] + s) * scale << "  " << (x[1] + s) * scale << "  ";
        file << "seg_blue"     << std::endl;		

		file << (x[0] + s) * scale << "  " << (x[1] - s) * scale << "  ";
        file << (x[0] - s) * scale << "  " << (x[1] + s) * scale << "  ";
        file << "seg_blue"     << std::endl;    
	}

	// plot domain boundaries
	size_t num_boundary_edges = mps_solver->_boundary_faces->size();
	for (size_t iedge = 0; iedge < num_boundary_edges; iedge++)
	{
		size_t ist = (*(mps_solver->_boundary_faces))[iedge][0];
		size_t iend = (*(mps_solver->_boundary_faces))[iedge][1];
		double xst = (*(mps_solver->_boundary_points))[ist][0];
		double yst = (*(mps_solver->_boundary_points))[ist][1];
		double xend = (*(mps_solver->_boundary_points))[iend][0];
		double yend = (*(mps_solver->_boundary_points))[iend][1];

		file << xst * scale << "  " << yst * scale << "  ";
        file << xend * scale << "  " << yend * scale << "  ";
        file << "seg_bold"     << std::endl;    
	}

    file << 0.0 * scale << "  " << 0.0 * scale << "  ";  
    file << 1.0 * scale << "  " << 0.0 * scale << "  ";          
	file << 1.0 * scale << "  " << 1.0 * scale << "  ";          
	file << 0.0 * scale << "  " << 1.0 * scale << "  ";          
    file << "quad_bold"      << std::endl;
	

    file << "showpage" << std::endl;
  
	return 0;
    #pragma endregion
};

