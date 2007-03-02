/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Ifpack_UserPartitioner.h"
#include "Epetra_CrsGraph.h"

//==============================================================================
int Ifpack_UserPartitioner::ComputePartitions()
{
  
  if (Map_ == 0)
    IFPACK_CHK_ERR(-1);

  // simply copy user's vector
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    Partition_[i] = Map_[i];
  }

  // put together all partitions composed by 1 one vertex
  // (if any)
  vector<int> singletons(NumLocalParts());
  for (unsigned int i = 0 ; i < singletons.size() ; ++i) {
    singletons[i] = 0;
  }

#if 0
  // may want to uncomment the following to ensure that no
  // partitions are in fact singletons
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    ++singletons[Partition_[i]];
  }
  
  int count = 0;
  for (unsigned int i = 0 ; i < singletons.size() ; ++i) {
    if (singletons[i] == 1)
      ++count;
  }

  int index = -1;
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    int j = Partition_[i];
    if (singletons[j] == 1) {
      if (index == -1)
        index = j;
      else
        Partition_[i] = index;
    }
  }
#endif

  return(0);
}
