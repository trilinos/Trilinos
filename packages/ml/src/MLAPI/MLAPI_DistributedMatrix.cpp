/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "ml_include.h"
//#include "ml_lapack.h"
#include "ml_comm.h"
#include "ml_epetra.h"
#include "MLAPI_Error.h"
#include "MLAPI_BaseObject.h"
#include "MLAPI_Space.h"
#include "MLAPI_DistributedMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include <iomanip>

using namespace std;

ostream& MLAPI::DistributedMatrix::
Print(ostream& os, const bool verbose) const
{
  int Length = MaxNumEntries();
  vector<double> Values(Length);
  vector<int>    Indices(Length);

  if (GetMyPID() == 0) {

    os << endl;
    os << "*** MLAPI::DistributedMatrix ***" << endl;
    os << "Label = " << GetLabel() << endl;
    os << "Number of rows    = " << GetRangeSpace().GetNumGlobalElements() << endl;
    os << "Number of columns = " << GetDomainSpace().GetNumGlobalElements() << endl;
    os << endl;
    os.width(10); os << "row ID";
    os.width(10); os << "col ID";
    os.width(30); os << "value";
    os << endl;
    os << endl;
  }

  for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

    if (GetMyPID() == iproc) {

      if (IsFillCompleted()) {
        for (int i = 0 ; i < NumMyRows() ; ++i) {
          int GRID = RangeMap_->GID(i);
          double* Values;
          int* Indices;
          int NumEntries;
          Matrix_->ExtractMyRowView(i, NumEntries, Values, Indices);
          for (int j = 0 ; j < NumEntries ; ++j) {
            os.width(10); os << GRID;
            os.width(10); os << Matrix_->RowMatrixColMap().GID(Indices[j]);
            os.width(30); os << Values[j];
            os << endl;
          }
        }
      }
      else {
        for (int i = 0 ; i < NumMyRows() ; ++i) {
          int GRID = RangeMap_->GID(i);
          double* Values;
          int* Indices;
          int NumEntries;
          Matrix_->ExtractGlobalRowView(GRID, NumEntries, Values, Indices);
          for (int j = 0 ; j < NumEntries ; ++j) {
            os.width(10); os << GRID;
            os.width(10); os << Indices[j];
            os.width(30); os << Values[j];
            os << endl;
          }
        }
      }
    }
    Barrier();
  }

  if (GetMyPID() == 0)
    os << endl;

  Barrier();

  return(os);
}

#endif
