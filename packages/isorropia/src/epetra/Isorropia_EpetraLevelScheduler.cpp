//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
//@HEADER

#include <Isorropia_EpetraLevelScheduler.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_CrsGraph.h>
#endif

namespace Isorropia {

#ifdef HAVE_EPETRA

namespace Epetra {


LevelScheduler::LevelScheduler(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
		 bool compute_now) : Operator(input_graph)
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
  properties_.clear();

  int nrows = graph.NumMyRows();
  int maxNonZeros = graph.MaxNumIndices();
  
  if ((nrows < 2) || (maxNonZeros < 1)){
    properties_.assign(nrows,1);
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
        if (properties_[IDs[j]] > depth)
          depth = properties_[IDs[j]];
      }
      depth++;
  
      properties_.push_back(depth);
    }
  }
  else if (graph.UpperTriangular()){

    properties_.assign(nrows,0);

    for (int i=nrows-1; i >= 0 ; i--){
  
      int numIDs = 0;
      int *IDs = NULL;
      graph.ExtractMyRowView(i, numIDs, IDs);
  
      int depth = -1;
  
      for (int j=0; j < numIDs; j++){
        if (properties_[IDs[j]] > depth)
          depth = properties_[IDs[j]];
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

