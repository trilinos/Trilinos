#ifndef PHX_CORE_UTILS_H
#define PHX_CORE_UTILS_H

#include "Epetra_Map.h"

#include "Teuchos_Hashtable.hpp"

namespace phx {
namespace core {

class Utils 
{
public:
  Utils() {}
  ~Utils() {}

  static int getNumDimensions()
  {
    return(numDimensions_);
  }

  static void setNumDimensions(const int numDimensions)
  {
    numDimensions_ = numDimensions;
  }

  static int numDimensions_;

#if 0
  static
  RefCountPtr<Epetra_Map>
  createMatrixMap(const Epetra_Comm& Comm,
                  map<string, RefCountPtr<phx::grid::Loadable> > patches)
  {
    int allocated = 0;
    map<string, RefCountPtr<phx::grid::Loadable> >::iterator iter;

    for (iter = patches.begin(); iter != patches.end(); ++iter)
    {
      allocated += iter->second->getNumGlobalVertices();
    }

    Teuchos::Hashtable<int, short int> hash(allocated * 2);

    vector<int> list(allocated * 2);

    int count = 0;

    for (iter = patches.begin(); iter != patches.end(); ++iter)
    {
      const RefCountPtr<Epetra_Map> vertexMap = iter->second->getVertexMap();

      int* myGlobalElements = vertexMap->MyGlobalElements();

      for (int i = 0; i < vertexMap->NumMyElements(); ++i)
      {
        if (!hash.containsKey(myGlobalElements[i]))
        {
          list[count++] = myGlobalElements[i];
          hash.put(myGlobalElements[i], 1);
        }
      }
    }

    return(rcp(new Epetra_Map(-1, count, &list[0], 0, Comm)));
  }


  // FIXME: DELETE
  static void 
  setDirichletBoundaryConditions(int size, int* oundaryList,
                                 double* boundaryValue, char* boundaryType,
                                 Epetra_FECrsMatrix& matrix, 
                                 Epetra_FEVector& rhs)
  {
  /*
  RefCountPtr<Epetra_Map> boundaryMap = rcp(new Epetra_Map(0, list.size(), &list[0], 0, comm));
  RefCountPtr<Epetra_Vector> boundaryValues = rcp(new Epetra_Vector(*boundaryMap));

  if (comm.MyPID() == 0)
  {
    boundaryValues[0] = values[0];
  }
  else if (comm.MyPID() == comm.NumProc())
  {
    boundaryValues[0] = values[0];
  }
  */

    Teuchos::Hashtable<int, double> dirichletRows;
    dirichletRows.put(0, 0.0);
    dirichletRows.put(matrix.NumGlobalRows() - 1, 0.0);

    for (int i = 0; i < matrix.NumMyRows(); ++i)
    {
      int GID = matrix.RowMatrixRowMap().GID(i);

      bool isDirichlet = false;
      if (dirichletRows.containsKey(GID)) isDirichlet = true;

      int* indices;
      double* values;
      int numEntries;
      matrix.ExtractMyRowView(i, numEntries, values, indices);

      if (isDirichlet)
      {
        for (int j = 0; j < numEntries; ++j)
          if (indices[j] != i) values[j] = 0.0;
          else values[j] = 1.0;
          rhs[0][i] = dirichletRows.get(GID);
      }
      else
      {
        for (int j = 0; j < numEntries; ++j)
        {
          if (indices[j] == i) continue;
          if (dirichletRows.containsKey(matrix.RowMatrixColMap().GID(indices[j]))) values[j] = 0.0;
        }
      }
    }
  }
#endif

private:

}; // class Utils

int Utils::numDimensions_ = 3;

} // namespace core
} // namespace phx

#endif
