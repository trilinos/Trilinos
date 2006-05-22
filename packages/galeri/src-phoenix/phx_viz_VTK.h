#ifndef PHX_VIZ_VTK_H
#define PHX_VIZ_VTK_H

#include "phx_grid_Element.h"

#include "Epetra_Import.h"

namespace phx {
namespace viz {

class VTK
{
public:

  static
  void Write(const Epetra_Comm& comm,
             phx::grid::Loadable& patch,
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
      if (patch.getElement().getLabel() == "phx::grid::Quad")
        what = "9 ";
      else if (patch.getElement().getLabel() == "phx::grid::Triangle")
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
} // namespace phx

#endif
