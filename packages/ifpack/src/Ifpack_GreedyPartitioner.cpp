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
#include "Ifpack_GreedyPartitioner.h"
#include "Ifpack_Graph.h"

#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Teuchos_ParameterList.hpp"

//==============================================================================
int Ifpack_GreedyPartitioner::ComputePartitions()
{
  std::vector<int> ElementsPerPart(NumLocalParts());
  std::vector<int> Count(NumLocalParts());
  for (int i = 0 ; i < NumLocalParts() ; ++i)
    Count[i] = 0;

  // define how many nodes have to be put on each part
  int div = NumMyRows() / NumLocalParts();
  int mod = NumMyRows() % NumLocalParts();

  for (int i = 0 ; i < NumLocalParts() ; ++i) {
    Count[i] = 0;
    ElementsPerPart[i] = div;
    if (i < mod) ElementsPerPart[i]++;
  }

  for( int i=0 ; i<NumMyRows() ; ++i ) {
    Partition_[i] = -1;
  }

  int NumEntries;
  std::vector<int> Indices(MaxNumEntries());
  
  // load root node for partition 0
  int CurrentPartition = 0;
  int TotalCount = 0;

  // filter singletons and empty rows, put all of them in partition 0
  for (int i = 0 ; i < NumMyRows() ; ++i) {
    NumEntries = 0;
    int ierr = Graph_->ExtractMyRowCopy(i, MaxNumEntries(),
                                        NumEntries, &Indices[0]);
    IFPACK_CHK_ERR(ierr);
    if (NumEntries <= 1) {
      Partition_[i] = 0;
      TotalCount++;
    }
  }

  if (TotalCount)
    CurrentPartition = 1;

  std::vector<int> ThisLevel(1);
  ThisLevel[0] = RootNode_;

  // be sure that RootNode is not a singleton or empty row
  if (Partition_[RootNode_] != -1) {
    // look for another RN
    for (int i = 0 ; i < NumMyRows() ; ++i)
      if (Partition_[i] == -1) {
        ThisLevel[0] = i;
        break;
      }
  }
  else {
    Partition_[RootNode_] = CurrentPartition;
  }

  // now aggregate the non-empty and non-singleton rows
  while (ThisLevel.size()) {

    std::vector<int> NextLevel;

    for (unsigned int i = 0 ; i < ThisLevel.size() ; ++i) {

      int CurrentNode = ThisLevel[i];
      int ierr = Graph_->ExtractMyRowCopy(CurrentNode, MaxNumEntries(),
                                          NumEntries, &Indices[0]);
      IFPACK_CHK_ERR(ierr);

      if (NumEntries <= 1)
        continue;

      for (int j = 0 ; j < NumEntries ; ++j) {

        int NextNode = Indices[j];
        if (NextNode >= NumMyRows()) continue;

        if (Partition_[NextNode] == -1) {
          // this is a free node
          NumLocalParts_ = CurrentPartition + 1;
          Partition_[NextNode] = CurrentPartition;
          ++Count[CurrentPartition];
          ++TotalCount;
          NextLevel.push_back(NextNode);
        }
      }
    } // for (i)

    // check whether change partition or not
    if (Count[CurrentPartition] >= ElementsPerPart[CurrentPartition])
      ++CurrentPartition;

    // swap next and this
    ThisLevel.resize(0);
    for (unsigned int i = 0 ; i < NextLevel.size() ; ++i)
      ThisLevel.push_back(NextLevel[i]);

    if (ThisLevel.size() == 0 && (TotalCount != NumMyRows())) {
      // need to look for new RootNode, do this in a simple way
      for (int i = 0 ; i < NumMyRows() ; i++) {
        if (Partition_[i] == -1)
          ThisLevel.push_back(i);
        break;
      }
    }

  } // while (ok)

  return(0);
}

