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
#include "Galeri_Exception.h"
#include "Galeri_FiniteElements.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

using namespace Galeri;
using namespace Galeri::FiniteElements;

// =========== //
// main driver //
// =========== //

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  try {
    
    // Prepares the computational domain. For simplicity,           
    // the computation domain has (nx * NumProcs, ny, nz) elements, 
    // and it is partitioned into (NumProcs, 1, 1) subdomains. 
    
    int nx = 2 * Comm.NumProc();
    int ny = 2; 
    int nz = 2;
    int mx = Comm.NumProc(), my = 1, mz = 1;
    double lx = 1.0 * Comm.NumProc();
    double ly = 1.0;
    double lz = 1.0;

    // Uncomment one of the following:
    TriangleRectangleGrid Grid(Comm, nx, ny, mx, my, lx, ly);
    //QuadRectangleGrid Grid(Comm, nx, ny, mx, my, lx, ly);
    //TetCubeGrid Grid(Comm, nx, ny, nz, mx, my, mz, lx, ly, lz);
    //HexCubeGrid Grid(Comm, nx, ny, nz, mx, my, mz, lx, ly, lz);

    // just work on processor 0, only to reduce the output
    if (!Comm.MyPID()) 
    {
      // Extracts the information from the Grid using AbstractGrid 
      // methods. First, some general information. 

      cout << "Number of dimensions = " << Grid.NumDimensions() << endl;
      cout << "Number of vertices per element = " << Grid.NumVerticesPerElement() << endl;
      cout << "Number of faces per element = " << Grid.NumFacesPerElement() << endl;
      cout << "Number of vertices per face = " << Grid.NumVerticesPerFace() << endl;
      cout << "Element type = " << Grid.ElementType() << endl;

      cout << "Number of elements: global = " << Grid.NumGlobalElements();
      cout << ", on proc 0 = " << Grid.NumMyElements() << endl;
      cout << "Number of vertices: global = " << Grid.NumGlobalVertices();
      cout << ", on proc 0 = " << Grid.NumMyVertices() << endl;
      cout << "Number of boundary faces: "
        "on proc 0 = " << Grid.NumMyBoundaryFaces() << endl;

      // Now loop over all the local elements, and print  
      // the local ID of each vertex.

      vector<double> coord(3);
      vector<int> vertices(Grid.NumVerticesPerElement());

      for (int i = 0 ; i < Grid.NumMyElements() ; ++i) 
      {
        cout << "Element " << i << ": ";
        Grid.ElementVertices(i, &vertices[0]);
        for (int j = 0 ; j < Grid.NumVerticesPerElement() ; ++j)
          cout << vertices[j] << " ";
        cout << endl;
      }

      // Loop over all local vertices, and print out the coordinates

      for (int i = 0 ; i < Grid.NumMyVertices() ; ++i)
      {
        Grid.VertexCoord(i, &coord[0]);
        cout << "Vertex " << i << ": (" << coord[0];
        cout << ", " << coord[1] << ", " << coord[2] << ")" << endl;
      }

      // Loop over all local boundary faces, and print out the local ID of each
      // vertex belonging to the face

      for (int i = 0 ; i < Grid.NumMyBoundaryFaces() ; ++i)
      {
        cout << "Boundary face " << i << ": ";
        int tag;
        Grid.FaceVertices(i, tag, &vertices[0]);
        for (int j = 0 ; j < Grid.NumVerticesPerFace() ; ++j)
          cout << vertices[j] << " ";
        cout << endl;
      }
    }
    
    // The following instructions can be used to visualize with MEDIT
#if 0
    MEDITInterface MEDIT(Comm);

    // We need to define a vector to plot, in this case constant
    Epetra_Vector Vector(Grid.RowMap());
    MEDIT.Write(Grid, "grid", Vector);
#endif

  }
  catch (Exception& rhs)
  {
    if (Comm.MyPID() == 0)
    {
      cerr << "Caught exception: ";
      rhs.Print();
    }
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
}
