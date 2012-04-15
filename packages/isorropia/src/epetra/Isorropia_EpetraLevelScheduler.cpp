//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER

#include <Isorropia_EpetraLevelScheduler.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_CrsGraph.h>
#endif

namespace Isorropia {

#ifdef HAVE_EPETRA

namespace Epetra {


LevelScheduler::~LevelScheduler() {}


LevelScheduler::LevelScheduler(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			       const Teuchos::ParameterList& paramlist,
		               bool compute_now) 
  : Operator(input_graph, paramlist,0)
{
  if (compute_now)
    schedule(true);
}

void
LevelScheduler::schedule(bool force_scheduling)
{
  if (alreadyComputed() && !force_scheduling)
    return;

  const Epetra_CrsGraph &graph = *input_graph_;

  int nrows = graph.NumMyRows();
  int maxNonZeros = graph.MaxNumIndices();

  properties_.clear();
  properties_.assign(nrows, 0);
  
  if ((nrows < 2) || (maxNonZeros < 1)){
    computeNumberOfProperties();
    return;
  }

  // algorithm from legacy Petra_CRS_Graph.cc 

  if (graph.LowerTriangular()){
    for (int i=0; i < nrows ; i++){
  
      int numIDs = 0;
      int *IDs = NULL;
      graph.ExtractMyRowView(i, numIDs, IDs);
  
      int depth = -1;
  
      for (int j=0; j < numIDs; j++){
        int col = IDs[j];
        if ((col < i) && (properties_[col] > depth))
          depth = properties_[col];
      }
      depth++;
  
      properties_[i] = depth;
    }
  }
  else if (graph.UpperTriangular()){

    for (int i=nrows-1; i >= 0 ; i--){
  
      int numIDs = 0;
      int *IDs = NULL;
      graph.ExtractMyRowView(i, numIDs, IDs);
  
      int depth = -1;
  
      for (int j=0; j < numIDs; j++){
        int col = IDs[j];
        if ((col > i) && (properties_[col] > depth))
          depth = properties_[col];
      }
      depth++;
  
      properties_[i] = depth;
    }
  }
  else{
    // error
    properties_.assign(nrows,-1);
  }

  operation_already_computed_ = true;

  computeNumberOfProperties();
}

} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

