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
#include "Ifpack_Graph.h"

#include "Epetra_Comm.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"

static const string PrintMsg_ = "(Ifpack_OvPartitioner) ";

//==============================================================================
Ifpack_OverlappingPartitioner::
Ifpack_OverlappingPartitioner(const Ifpack_Graph* Graph) :
  NumLocalParts_(1),
  Graph_(Graph),
  OverlappingLevel_(0),
  IsComputed_(false),
  verbose_(false)
{
}

//==============================================================================
Ifpack_OverlappingPartitioner::~Ifpack_OverlappingPartitioner()
{
}

//==============================================================================
int Ifpack_OverlappingPartitioner::SetParameters(Teuchos::ParameterList& List)
{

  NumLocalParts_ = List.get("partitioner: local parts", NumLocalParts_);
  OverlappingLevel_ = List.get("partitioner: overlap", OverlappingLevel_);
  verbose_ = List.get("partitioner: print level", verbose_);

  if (NumLocalParts_ < 0)
    NumLocalParts_ = Graph_->NumMyRows() / (-NumLocalParts_);
  if (NumLocalParts_ == 0)
    NumLocalParts_ = 1;
  if (NumLocalParts_ < 0)
    IFPACK_CHK_ERR(-1);
  if (NumLocalParts_ > Graph_->NumMyRows())
    IFPACK_CHK_ERR(-1);

  if (OverlappingLevel_ < 0)
    IFPACK_CHK_ERR(-1);

  SetPartitionParameters(List);

  return(0);
}

//==============================================================================
int Ifpack_OverlappingPartitioner::Compute()
{

  if (NumLocalParts_ < 1)
    IFPACK_CHK_ERR(-1); // incorrect value

  if (OverlappingLevel_ < 0)
    IFPACK_CHK_ERR(-1); // incorrect value

  // some output

  if (verbose_ && (Comm().MyPID() == 0)) {
    cout << PrintMsg_ << "Number of local parts  = " << NumLocalParts_ << endl;
    cout << PrintMsg_ << "Number of global parts = " 
         << NumLocalParts_ * Comm().NumProc() << endl;
    cout << PrintMsg_ << "Amount of overlap      = " << OverlappingLevel_ << endl;
  }

  // 1.- allocate memory 

  Partition_.resize(NumMyRows());
  Parts_.resize(NumLocalParts());

  // 2.- sanity checks on input graph
 
  if (Graph_->Filled() == false)
    IFPACK_CHK_ERR(-4); // need FillComplete() called

  if (Graph_->NumGlobalRows64() != Graph_->NumGlobalCols64())
    IFPACK_CHK_ERR(-3); // can partition square matrices only

  if (NumLocalParts_ < 1)
    IFPACK_CHK_ERR(-2); // value not valid
 
  // 3.- perform non-overlapping partition
 
  IFPACK_CHK_ERR(ComputePartitions());

  // 4.- compute the partitions with overlapping
  
  IFPACK_CHK_ERR(ComputeOverlappingPartitions());

  // 5.- return to the user
 
  IsComputed_ = true;

  return(0);
}

// ======================================================================
int Ifpack_OverlappingPartitioner::ComputeOverlappingPartitions()
{

  // FIXME: the first part of this function should be elsewhere
  // start defining the subgraphs for no overlap

  std::vector<int> sizes;
  sizes.resize(NumLocalParts_);

  // 1.- compute how many rows are in each subgraph
  for (int i = 0 ; i < NumLocalParts_ ; ++i)
    sizes[i] = 0;

  for (int i = 0 ; i < NumMyRows() ; ++i) {
    if (Partition_[i] >= NumLocalParts_) {
      cerr << "ERROR: Partition[" << i << "] = "<< Partition_[i] 
	   << ", NumLocalParts = " << NumLocalParts_ << endl;
      cerr << "(file = " << __FILE__ << ", line = "
           << __LINE__ << ")" << endl;
      IFPACK_CHK_ERR(-10);
    }
    // no singletons should be here, as the matrix is
    // supposed to be filtered through Ifpack_SingletonFilter
    if (Partition_[i] == -1)
      IFPACK_CHK_ERR(-1);
    sizes[Partition_[i]]++;
  }

  // 2.- allocate space for each subgraph
  for (int i = 0 ; i < NumLocalParts_ ; ++i)
    Parts_[i].resize(sizes[i]);

  // 3.- cycle over all rows and populate the vectors

  for (int i = 0 ; i < NumLocalParts_ ; ++i)
    sizes[i] = 0;

  for (int i = 0 ; i < NumMyRows() ; ++i) {
    int part = Partition_[i];
    int count = sizes[part];
    Parts_[part][count] = i;
    sizes[part]++;
  }

  if (OverlappingLevel_ == 0)
    return(0);

  // wider overlap requires further computations
  for (int level = 1 ; level <= OverlappingLevel_ ; ++level) {

    std::vector<std::vector<int> > tmp;
    tmp.resize(NumLocalParts_);

    // cycle over all rows in the local graph (that is the overlapping
    // graph). For each row, all columns will belong to the subgraph of
    // row `i'.

    int MaxNumEntries_tmp = Graph_->MaxMyNumEntries();
    std::vector<int> Indices;
    Indices.resize(MaxNumEntries_tmp);

    for (int part = 0 ; part < NumLocalParts_ ; ++part) {

      for (int i = 0; i < (int)Parts_[part].size() ; ++i) {  

	int LRID = Parts_[part][i];
	int NumIndices;
	int ierr = Graph_->ExtractMyRowCopy(LRID, MaxNumEntries_tmp, 
                                            NumIndices, &Indices[0]);
	IFPACK_CHK_ERR(ierr);

	for (int j = 0 ; j < NumIndices ; ++j) {

	  // use *local* indices
	  int col = Indices[j];
          if (col >= NumMyRows())
            continue;

	  // has this column already been inserted?
	  std::vector<int>::iterator
	    where = find(tmp[part].begin(), tmp[part].end(), col);

	  if (where == tmp[part].end()) {
	    tmp[part].push_back(col);
	  }
	}
      }
    }

    // now I convert the STL vectors into Epetra_IntSerialDenseVectors.
    for (int i = 0 ; i < NumLocalParts_ ; ++i) {
      Parts_[i].resize(tmp[i].size());
      for (int j = 0 ; j < (int)tmp[i].size() ; ++j)
	Parts_[i][j] = tmp[i][j];
    }
  }

  return(0);

}

//============================================================================
int Ifpack_OverlappingPartitioner::NumMyRows() const
{
  return(Graph_->NumMyRows());
}

//============================================================================
int Ifpack_OverlappingPartitioner::NumMyNonzeros() const
{
  return(Graph_->NumMyNonzeros());
}

//============================================================================
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
int Ifpack_OverlappingPartitioner::NumGlobalRows() const
{
  return(Graph_->NumGlobalRows());
}
#endif

long long Ifpack_OverlappingPartitioner::NumGlobalRows64() const
{
  return(Graph_->NumGlobalRows64());
}
//============================================================================
int Ifpack_OverlappingPartitioner::MaxNumEntries() const
{
  return(Graph_->MaxMyNumEntries());
}

//============================================================================
const Epetra_Comm& Ifpack_OverlappingPartitioner::Comm() const
{
  return(Graph_->Comm());
}

// ======================================================================
ostream& Ifpack_OverlappingPartitioner::Print(ostream & os) const
{

  if (Comm().MyPID()) 
    return(os);

  os << "================================================================================" << endl;
  os << "Ifpack_OverlappingPartitioner" << endl;
  os << "Number of local rows  = " << Graph_->NumMyRows() << endl;
  os << "Number of global rows = " << Graph_->NumGlobalRows64() << endl;
  os << "Number of local parts = " << NumLocalParts_ << endl;
  os << "Overlapping level     = " << OverlappingLevel_ << endl;
  os << "Is computed           = " << IsComputed_ << endl;
  os << "================================================================================" << endl;

  return(os);
}

