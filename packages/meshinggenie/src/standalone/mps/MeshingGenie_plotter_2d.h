
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

   output: a postscript file that shows various details

 * Last modified: 11/21/2012
********************************************************************************/


#ifndef MESHING_GENIE_PLOTTER_2D_H
#define MESHING_GENIE_PLOTTER_2D_H

#include "MeshingGenie_defs.h"
#include "MeshingGenie_mps_nd.h"

class MeshingGenie_plotter_2d   
{	
    public:
        //! constructor
		MeshingGenie_plotter_2d(){ };
   

        //! Destructor                
        ~MeshingGenie_plotter_2d(){ };

		int plot_active_pool(MeshingGenie_mps_nd* mps_solver);
};
                                
#endif	



