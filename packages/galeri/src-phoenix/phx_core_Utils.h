#ifndef PHX_CORE_UTILS_H
#define PHX_CORE_UTILS_H

#include "phx_grid_Loadable.h"

#include "Epetra_Map.h"

#include "Teuchos_Hashtable.hpp"

namespace phx {
namespace core {

class Utils 
{
public:
  Utils() {}
  ~Utils() {}

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
      const Epetra_Map& vertexMap = iter->second->getVertexMap();

      int* myGlobalElements = vertexMap.MyGlobalElements();

      for (int i = 0; i < vertexMap.NumMyElements(); ++i)
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

private:

}; // class Utils

} // namespace core
} // namespace phx

#endif
