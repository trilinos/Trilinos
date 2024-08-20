// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
