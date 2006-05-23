#ifndef PHX_VIZ_MEDIT_H
#define PHX_VIZ_MEDIT_H

#include "phx_grid_Element.h"

#include "Epetra_Import.h"

namespace phx {
namespace viz {

class MEDIT
{
public:

  static
  void Write(const Epetra_Comm& comm,
             phx::grid::Loadable& patch,
             const string& BaseName,
             const Epetra_MultiVector& vector)
  {
    std::ofstream medit;
    string FileName = BaseName + ".mesh";

    const Epetra_MultiVector& linearCoord = patch.getLinearCoordinates();

    if (comm.MyPID() == 0)
    {
      medit.open(FileName.c_str());
      medit << "MeshVersionFormatted 1" << endl;
      medit << "Dimension 3" << endl;
      medit << "# mesh from phoenix" << endl << endl;
      medit << "Vertices " << patch.getNumGlobalVertices() << endl;
      medit.close();
    }

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    {
      if (iproc == comm.MyPID())
      {
        medit.open(FileName.c_str(),ios::app);
        if (phx::core::Utils::getNumDimensions() == 2)
        {
          for (int i = 0; i < linearCoord.MyLength(); ++i)
          {
            medit << setw(12) << setiosflags(ios::showpoint) 
              << setw(12) << linearCoord[0][i] << " "
              << setw(12) << linearCoord[1][i] << " "
              << setw(12) << "0.0 1" << endl;
          }
        }
        else if (phx::core::Utils::getNumDimensions() == 3)
        {
          for (int i = 0; i < linearCoord.MyLength(); ++i)
          {
            medit << setw(12) << setiosflags(ios::showpoint) 
              << setw(12) << linearCoord[0][i] << " "
              << setw(12) << linearCoord[1][i] << " "
              << setw(12) << linearCoord[2][i] << " 1" << endl;
          }
        }
          
        medit.close();
      }
      comm.Barrier();
    }

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    {
      if (comm.MyPID() == iproc) 
      {
        medit.open(FileName.c_str(),ios::app);
        if (iproc == 0)
        {
          if (patch.getElementType() == "Triangle")
            medit << "Triangles " << patch.getNumGlobalElements() << endl;
          else if (patch.getElementType() == "Quad")
            medit << "Quadrilaterals " << patch.getNumGlobalElements() << endl;
          else if (patch.getElementType() == "Hex")
            medit << "Hexahedra " << patch.getNumGlobalElements() << endl;
          else
          {
            cout << "NOT SUPPORTED YET" << endl;
            exit(0);
          }
        }

        for (int i = 0 ; i < patch.getNumMyElements(); ++i) 
        {
          for (int j = 0; j < patch.getNumVerticesPerElement(); ++j)
            medit << patch.getMyConnectivity(i, j) + 1 << " ";

          medit << comm.MyPID() << endl;
        }

        if (iproc == comm.NumProc() - 1) 
          medit << endl << "End" << endl;

        medit.close();
      }

      comm.Barrier();
    }

    Epetra_Import importer(linearCoord.Map(), vector.Map());
    Epetra_MultiVector linearVector(linearCoord.Map(), vector.NumVectors());
    linearVector.Import(vector, importer, Insert);

    // ======== //
    // .bb file //
    // ======== //

    FileName = BaseName + ".bb";    

    for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
    {
      if (comm.MyPID() == iproc) 
      {
        if (iproc == 0)
        {
          medit.open(FileName.c_str());
          medit << "3 1 " << patch.getNumGlobalVertices() << " 2" << endl;
        }
        else
          medit.open(FileName.c_str(),ios::app);

        for (int i = 0; i < linearVector.MyLength(); ++i)
          medit << setiosflags(ios::showpoint) << linearVector[0][i] << endl;

        medit.close();
      }
      comm.Barrier();
    }
  }
};

} // namespace viz
} // namespace phx

#endif
