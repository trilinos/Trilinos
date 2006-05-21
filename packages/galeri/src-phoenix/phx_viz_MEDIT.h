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
    if (comm.MyPID() == 0)
    {
      string FileName = BaseName + ".mesh";
      std::ofstream medit;

      medit.open(FileName.c_str());
      medit << "MeshVersionFormatted 1" << endl;
      medit << "Dimension 3" << endl;
      medit << "# mesh from phoenix" << endl << endl;
      medit << "Vertices " << patch.getNumGlobalVertices() << endl;

      for (int i = 0; i < patch.getNumGlobalVertices(); ++i)
      {
        medit << setw(12) << setiosflags(ios::showpoint) 
              << setw(12) << patch.getGlobalCoordinates(i, 0) << " "
              << setw(12) << patch.getGlobalCoordinates(i, 1) << " "
              << setw(12) << "0.0 1" << endl;
      }
      medit.close();
    }

    comm.Barrier();

    for (int ProcID = 0 ; ProcID < comm.NumProc() ; ++ProcID) 
    {
      if (comm.MyPID() == ProcID) {

        string FileName = BaseName + ".mesh";
        std::ofstream medit(FileName.c_str(),ios::app);

        medit << "Quadrilaterals " << patch.getNumGlobalElements() << endl;

        for (int i = 0 ; i < patch.getNumMyElements(); ++i) 
        {
          for (int j = 0; j < 4; ++j)
            medit << patch.getGlobalConnectivity(i, j) + 1 << " ";

          medit << comm.MyPID() << endl;
        }

        if (ProcID == comm.NumProc() - 1) 
          medit << endl << "End" << endl;

        medit.close();
      }

      comm.Barrier();

    } // for Procs, write elements

    Epetra_Map linearMap(vector.GlobalLength(), 0, comm);
    Epetra_Import importer(linearMap, vector.Map());
    Epetra_MultiVector linearVector(linearMap, vector.NumVectors());
    linearVector.Import(vector, importer, Insert);

    // ======== //
    // .bb file //
    // ======== //

    if (comm.MyPID() == 0) {

      string BBName = BaseName + ".bb";    
      std::ofstream bb;

        bb.open(BBName.c_str());
        bb << "3 1 " << patch.getNumGlobalVertices() << " 2" << endl;

      for (int i = 0; i < linearVector.MyLength(); ++i)
        bb << setiosflags(ios::showpoint) << linearVector[0][i] << endl;

      bb.close();

    }
  }
};

} // namespace viz
} // namespace phx

#endif
