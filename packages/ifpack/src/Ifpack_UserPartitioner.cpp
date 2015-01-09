/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
  std::vector<int> singletons(NumLocalParts());
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
