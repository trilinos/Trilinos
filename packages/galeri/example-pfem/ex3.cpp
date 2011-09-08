// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

// ===========================================================================
// This example shows the usage of the grid generator of Galeri. In this
// case, we create a hexahedral grid for a cube. The internal elements are
// defined in the "domain" data structure, while the boundary elements in the
// "boundary" data structure. In both cases, grid elements are
// Galeri::grid::Loadable objects.
//
// The grid generator is quite simple, and works for the following shapes:
// - (0,1) x (0,1) square, divided into squares or triangles
// - (0,1) x (0,1) x (0,1) cube, divided into tets or hexs.
// ===========================================================================

#include "Galeri_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"

#include "Galeri_core_Object.h"
#include "Galeri_core_Workspace.h"
#include "Galeri_grid_Loadable.h"
#include "Galeri_grid_Generator.h"
#include "Galeri_viz_MEDIT.h"

using namespace Galeri;

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  Galeri::core::Workspace::setNumDimensions(3);

  Galeri::grid::Loadable domain, boundary;

  int numGlobalElementsX = 2 * comm.NumProc();
  int numGlobalElementsY = 2;
  int numGlobalElementsZ = 2;

  int mx = comm.NumProc();
  int my = 1;
  int mz = 1;

  Galeri::grid::Generator::
  getCubeWithHexs(comm, numGlobalElementsX, numGlobalElementsY, numGlobalElementsZ,
                  mx, my, mz, domain, boundary);

  Epetra_Vector x(domain.getVertexMap());
  x.Random();

  Galeri::viz::MEDIT::write(domain, "sol", x);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
