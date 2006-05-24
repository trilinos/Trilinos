#ifndef PHX_CORE_UTILS_H
#define PHX_CORE_UTILS_H

#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

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

    static 
    Epetra_MultiVector* createMultiVectorComponent(const Epetra_MultiVector& input)
    {
      const Epetra_Comm& comm = input.Comm();
      const Epetra_BlockMap& inputMap = input.Map();
      const int* myGlobalElements = inputMap.MyGlobalElements();
      const int numMyElements = inputMap.NumMyElements();
      const int indexBase = inputMap.IndexBase();

      Epetra_Map extractMap(-1, numMyElements, myGlobalElements, indexBase, comm);
      Epetra_MultiVector* output = new Epetra_MultiVector(extractMap, input.NumVectors());

      return(output);
    }

    static 
    void extractMultiVectorComponent(const Epetra_MultiVector& input,
                                     const int equation,
                                     Epetra_MultiVector& output)
    {
      const Epetra_BlockMap& inputMap = input.Map();
      const int numMyElements = inputMap.NumMyElements();

      for (int i = 0; i < numMyElements; ++i)
      {
        int j = inputMap.FirstPointInElement(i) + equation;
        for (int k = 0; k < input.NumVectors(); ++k)
          output[k][i] = input[k][j];
      }
    }

  private:

}; // class Utils

} // namespace core
} // namespace phx

#endif
