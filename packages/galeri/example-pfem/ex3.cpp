// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
