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

#ifndef GALERI_VIZ_VTK_H
#define GALERI_VIZ_VTK_H

#include "Epetra_IntVector.h"

#include "Galeri_grid_Element.h"

#include "Epetra_Import.h"

namespace Galeri {
namespace viz {

class VTK
{
public:

  static
  void write(Galeri::grid::Loadable& patch,
             const string& BaseName,
             const Epetra_MultiVector& vector)
  {
    const Epetra_Comm& comm = patch.getComm();

    string FileName = BaseName + ".vtu";
    std::ofstream vtkFile;

    const Epetra_Map& nonOverlappingVertexMap = patch.getNonOverlappingVertexMap();
    const Epetra_MultiVector& nonOverlappingCoord = patch.getNonOverlappingCoordinates();

    Epetra_IntVector itmp(nonOverlappingVertexMap);
    for (int i = 0; i < itmp.MyLength(); ++i)
      itmp[i] = i + 
      

    int numGlobalVertices = nonOverlappingVertexMap.NumGlobalElements();
    int numGlobalElements = patch.getNumGlobalElements();

    if (comm.MyPID() == 0)
    {
      vtkFile.open(FileName.c_str());

      vtkFile << "<?xml version=\"1.0\"?>" << endl;
      vtkFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">" << endl;
      vtkFile << "<UnstructuredGrid>" << endl;
      vtkFile << "<Piece NumberOfPoints=\"" << numGlobalVertices << "\" ";
      vtkFile << "NumberOfCells=\"" << numGlobalElements << "\">" << endl;
      vtkFile << "<Points>" << endl;
      vtkFile << "<DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
      vtkFile.close();
    }

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    {
      if (iproc == comm.MyPID())
      {
        vtkFile.open(FileName.c_str(),ios::app);

        for (int i = 0; i < nonOverlappingCoord.MyLength(); ++i)
        {
          vtkFile << setw(12) << setiosflags(ios::showpoint) 
            << setw(12) << nonOverlappingCoord[0][i] << " "
            << setw(12) << nonOverlappingCoord[1][i] << " "
            << setw(12) << "0.0" << endl;
        }
        vtkFile.close();
      }
      comm.Barrier();
    }

    if (comm.MyPID() == 0)
    {
      vtkFile.open(FileName.c_str(),ios::app);
      vtkFile << "</DataArray>" << endl;
      vtkFile << "</Points>" << endl;

      vtkFile << "<Cells>" << endl;
      vtkFile << "<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
      vtkFile.close();
    }

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    {
      if (iproc == comm.MyPID())
      {
        vtkFile.open(FileName.c_str(),ios::app);

        for (int i = 0 ; i < patch.getNumMyElements(); ++i) 
        {
          for (int j = 0; j < patch.getNumVerticesPerElement(); ++j)
            vtkFile << patch.getGlobalConnectivity(i, j) + 1 << " ";
          vtkFile << endl;
        }

        vtkFile.close();
      }
      comm.Barrier();
    }

    if (comm.MyPID() == 0)
    {
      vtkFile.open(FileName.c_str(),ios::app);

      vtkFile << "</DataArray>" << endl;
      vtkFile << "<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
      for (int i = 0 ; i < patch.getNumGlobalElements(); ++i) 
        vtkFile << patch.getNumVerticesPerElement() << " " << endl;

      vtkFile << "</DataArray>" << endl;
      vtkFile << "<DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;

      string what;
      if (patch.getElement().getLabel() == "Galeri::grid::Quad")
        what = "9 ";
      else if (patch.getElement().getLabel() == "Galeri::grid::Triangle")
        what = "8 "; // FIXME
      for (int i = 0 ; i < patch.getNumMyElements(); ++i) 
        vtkFile << what << endl;

      vtkFile << "</DataArray>" << endl;
      vtkFile << "</Cells>" << endl;

      vtkFile << "<PointData Vectors=\"Velocity\">" << endl;
      vtkFile << "<DataArray type=\"Float32\" Name=\"Displacement\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;

      vtkFile.close();
    }

    Epetra_Map nonOverlappingMap(vector.GlobalLength(), 0, comm);
    Epetra_Import importer(nonOverlappingMap, vector.Map());
    Epetra_MultiVector nonOverlappingVector(nonOverlappingMap, vector.NumVectors());
    nonOverlappingVector.Import(vector, importer, Insert);

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    {
      vtkFile.open(FileName.c_str(),ios::app);

      for (int i = 0; i < nonOverlappingVector.MyLength(); ++i)
        vtkFile << setiosflags(ios::showpoint) << nonOverlappingVector[0][i] << endl;

      vtkFile.close();
    }

    if (comm.MyPID() == 0)
    {
      vtkFile.open(FileName.c_str(),ios::app);

      vtkFile << "</DataArray>" << endl;
      vtkFile << "</PointData>" << endl;

      vtkFile << "</Piece>" << endl;
      vtkFile << "</UnstructuredGrid>" << endl;
      vtkFile << "</VTKFile>" << endl;

      vtkFile.close();
    }
  }
}; // class VTK

} // namespace viz
} // namespace Galeri

#endif
