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

#include "Galeri_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"

#include "Galeri_core_Workspace.h"
#include "Galeri_grid_Triangle.h"
#include "Galeri_grid_Quad.h"
#include "Galeri_grid_Loadable.h"
#include "Galeri_viz_MEDIT.h"

using namespace Galeri;

// ============================================================================
// When run on two processors, the code creates the following grid:
//
// [2] -- [5] -- [8]
//  |    / |    / |
//  |   /  |   /  |
//  |  /   |  /   |
//  | /    | /    |
// [1] -- [4] -- [7]
//  |      |      |
//  |      |      |
//  |      |      |
//  |      |      |
// [0] -- [3] -- [6]
//
// which is composed by 2 quadrilateral elements, and 4 triangular elements.
// The quads have size h x h. h is the minimal length of triangle sides as
// well. When run with more processors, additional elements are added
// on the right of the above picture, so that the total number of quads is
// comm.NumProc(), and that of triangles 2 * comm.NumProc().
//
// Then, one each processor we define a vector, living the of grid vertices,
// and we produce a VTK output file. This file can be visualized, for example,
// with MaYaVi.
//
// \author Marzio Sala, ETH
//
// \date Last modified on Aug-06.
// ============================================================================

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  Galeri::core::Workspace::setNumDimensions(2);

  int vertexOffset = 3 * comm.MyPID();
  int elementOffset = 3 * comm.MyPID();
  double h = 1.0;
  double xOffset = h * comm.MyPID();

  // domain1 contains the quads, domain2 the triangles
  grid::Loadable domain1(comm, -1, 1, "Quad", &elementOffset);
  grid::Loadable domain2(comm, -1, 2, "Triangle");

  // domain1, start with connectivity.
  domain1.setGlobalConnectivity(elementOffset, 0, vertexOffset + 0);
  domain1.setGlobalConnectivity(elementOffset, 1, vertexOffset + 3);
  domain1.setGlobalConnectivity(elementOffset, 2, vertexOffset + 4);
  domain1.setGlobalConnectivity(elementOffset, 3, vertexOffset + 1);

  domain1.freezeConnectivity();

  // x-coordinates
  domain1.setGlobalCoordinates(vertexOffset + 0, 0, xOffset);
  domain1.setGlobalCoordinates(vertexOffset + 3, 0, xOffset + h);
  domain1.setGlobalCoordinates(vertexOffset + 4, 0, xOffset + h);
  domain1.setGlobalCoordinates(vertexOffset + 1, 0, xOffset);

  // y-coordinates
  domain1.setGlobalCoordinates(vertexOffset + 0, 1, 0.0);
  domain1.setGlobalCoordinates(vertexOffset + 3, 1, 0.0);
  domain1.setGlobalCoordinates(vertexOffset + 4, 1, h);
  domain1.setGlobalCoordinates(vertexOffset + 1, 1, h);

  domain1.freezeCoordinates();

  cout << domain1;

  // now domain2, start with connectivity
  domain2.setGlobalConnectivity(elementOffset, 0, vertexOffset + 1);
  domain2.setGlobalConnectivity(elementOffset, 1, vertexOffset + 4);
  domain2.setGlobalConnectivity(elementOffset, 2, vertexOffset + 5);

  domain2.setGlobalConnectivity(elementOffset + 1, 0, vertexOffset + 1);
  domain2.setGlobalConnectivity(elementOffset + 1, 1, vertexOffset + 5);
  domain2.setGlobalConnectivity(elementOffset + 1, 2, vertexOffset + 2);

  domain2.freezeConnectivity();

  // x-coordinates
  domain2.setGlobalCoordinates(vertexOffset + 1, 0, xOffset);
  domain2.setGlobalCoordinates(vertexOffset + 4, 0, xOffset + h);
  domain2.setGlobalCoordinates(vertexOffset + 5, 0, xOffset + h);
  domain2.setGlobalCoordinates(vertexOffset + 2, 0, xOffset);

  // y-coordinates
  domain2.setGlobalCoordinates(vertexOffset + 1, 1, h);
  domain2.setGlobalCoordinates(vertexOffset + 4, 1, h);
  domain2.setGlobalCoordinates(vertexOffset + 5, 1, 2 * h);
  domain2.setGlobalCoordinates(vertexOffset + 2, 1, 2 * h);

  domain1.freezeCoordinates();

  // output random values defined on the grid vertices
  
  Epetra_Vector x1(domain1.getVertexMap()), x2(domain2.getVertexMap());
  x1.Random(), x2.Random();

  Galeri::viz::MEDIT::write(domain1, "domain1", x1);
  Galeri::viz::MEDIT::write(domain2, "domain2", x2);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
