// @HEADER
// ************************************************************************
//
//            Galeri: Finite Element and Matrix Generation Package
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

#ifndef GALERI_VIZ_VTK_H
#define GALERI_VIZ_VTK_H

#include "Galeri_grid_Element.h"

#include "Epetra_Import.h"

namespace Galeri {
namespace viz {

class VTK
{
public:

  static
  void Write(const Epetra_Comm& comm,
             Galeri::grid::Loadable& patch,
             const string& BaseName,
             const Epetra_MultiVector& vector)
  {
    string FileName = BaseName + ".vtu";
    std::ofstream vtkFile;

    int numGlobalVertices = patch.getNumGlobalVertices();
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

    const Epetra_MultiVector& linearCoord = patch.getLinearCoordinates();

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    {
      if (iproc == comm.MyPID())
      {
        vtkFile.open(FileName.c_str(),ios::app);

        for (int i = 0; i < linearCoord.MyLength(); ++i)
        {
          vtkFile << setw(12) << setiosflags(ios::showpoint) 
            << setw(12) << linearCoord[0][i] << " "
            << setw(12) << linearCoord[1][i] << " "
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

    Epetra_Map linearMap(vector.GlobalLength(), 0, comm);
    Epetra_Import importer(linearMap, vector.Map());
    Epetra_MultiVector linearVector(linearMap, vector.NumVectors());
    linearVector.Import(vector, importer, Insert);

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    {
      vtkFile.open(FileName.c_str(),ios::app);

      for (int i = 0; i < linearVector.MyLength(); ++i)
        vtkFile << setiosflags(ios::showpoint) << linearVector[0][i] << endl;

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
