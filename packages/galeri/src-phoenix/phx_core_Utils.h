// @HEADER
// ************************************************************************
//
//                  Galeri Matrix Generation Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef PHX_CORE_UTILS_H
#define PHX_CORE_UTILS_H

#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

#define PHX_MAX(x,y) (( (x) > (y) ) ? x : y)
#define PHX_MIN(x,y) (( (x) < (y) ) ? x : y) 
#define PHX_SGN(x) (((x) < 0.0) ? -1.0 : 1.0) 

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
