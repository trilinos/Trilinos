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
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

#include "Galeri_core_Workspace.h"
#include "Galeri_grid_Segment.h"
#include "Galeri_grid_Loadable.h"
//#include "Galeri_grid_Rebalance.h"

// ============================================================================ 
// The goal of this example is to show the basic usage of class
// Galeri::grid::Loadable. The example builds, step-by-step, a 1D grid,
// associates some values to elements and vertices, then plots a function
// defined on the grid vertices in MEDIT format.
//
// This example can be run with any numbers of processors.
//
// \author Marzio Sala, ETHZ
//
// \date Last modified on Aug-06
// ============================================================================ 
//
using namespace Galeri;

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // ======================================================= //
  // Specifies the dimensionality of the problem: 1, 2, or 3 //
  // Creates a 1D grid on (0, 1) composed by segments. Each  //
  // processor will have 4 elements. We build the grid using //
  // the constructor that accepts a string to specify the    //
  // element type. Custom-defined elements can be specified  //
  // by declaring an empty grid object, then calling method  //
  // initialize() with an instance of the element.           //
  //                                                         //
  // We allocate space to store one additional data per grid //
  // element, and 1 additional data per grid vertex.         //
  // ======================================================= //
  
  core::Workspace::setNumDimensions(1);

  int numMyElements = 4;
  int numGlobalElements = numMyElements * comm.NumProc();
  int numElementData = 1;
  int numVertexData = 1;

  grid::Loadable domain(comm, numGlobalElements, numMyElements, 
                        "Segment", 0, numElementData, numVertexData);

  // ===================================================== //
  // Each processor inserts locally owned elements, then   //
  // call freezeCoordinates(), which computes the set      //
  // of locally owned vertices. Then, the coordinates of   //
  // these vertices is inserted. Finanlly, the grid is     //
  // freezed by calling freezeConnectivity().              //
  // ===================================================== //
  
  for (int LEID = 0; LEID < numMyElements; ++LEID)
  {
    int GEID = domain.getGEID(LEID);

    domain.setGlobalConnectivity(GEID, 0, GEID);
    domain.setGlobalConnectivity(GEID, 1, GEID + 1);

    domain.setElementData(GEID, 0, comm.MyPID()); // element data is proc numbere
  }

  domain.freezeConnectivity();

  double h = 1.0 / numGlobalElements;

  for (int LVID = 0; LVID < domain.getNumMyVertices(); ++LVID)
  {
    int GVID = domain.getGVID(LVID);
    domain.setGlobalCoordinates(GVID, 0, h * GVID);
    domain.setVertexData(GVID, 0, 1.0); // store something on each vertex
  }

  domain.freezeCoordinates();

  // prints out the grid data
  cout << domain;

  // ========================================================= //
  // We now create the set of boundary nodes. For simplicity,  //
  // both nodes (the one on the left and the one on the right) //
  // are defined on processor 0. The nodes have coordinates    //
  // (0.0) and (1.0), and they are of Dirichlet type.          //
  // No data is assigned to elements and vertices.             //
  //                                                           //
  // NOTE: boundary data are just another collection of grid   //
  //       elements! In this case, elements are "points".      //
  //                                                           //
  // NOTE: there is no connection between the data layout used //
  //       for the internal elements, and the one used for the //
  //       boundary elements. In fact, the two data structures //
  //       are completely independent. This makes the code     //
  //       more flexible, and boundary data easier to insert.  //
  //       One can have as many data structures as required,   //
  //       for internal or border elements. Mixed elements can //
  //       be supported by creating a different grid object    //
  //       for each element type.                              //
  //       However, the global vertex ID must be maintained    //
  //       over all grid objects.                              //
  // ========================================================= //

  int numMyBoundaryElements = (comm.MyPID() == 0)?2:0; 

  grid::Loadable boundary(comm, 2, numMyBoundaryElements, "Point");

  if (comm.MyPID() == 0)
  {
    boundary.setGlobalConnectivity(0, 0, 0);
    boundary.setGlobalConnectivity(1, 0, domain.getNumGlobalElements());
  }

  boundary.freezeConnectivity();

  if (comm.MyPID() == 0)
  {
    boundary.setGlobalCoordinates(0, 0, 0.0);
    boundary.setGlobalCoordinates(domain.getNumGlobalElements(), 0, 1.0);
  }

  boundary.freezeCoordinates();

  cout << boundary;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
