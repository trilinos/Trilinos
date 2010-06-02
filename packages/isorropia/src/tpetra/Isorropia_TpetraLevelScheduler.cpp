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

#include <Isorropia_TpetraLevelScheduler.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_ISORROPIA_TPETRA
#include <Tpetra_CrsGraph_decl.hpp>

namespace Isorropia {

namespace Tpetra {


template <class Node>
LevelScheduler<Node>::~LevelScheduler() {}


template <class Node>
LevelScheduler<Node>::LevelScheduler(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph,
			       const Teuchos::ParameterList& paramlist,
		               bool compute_now) 
  : Isorropia::Tpetra::Operator<Node>(input_graph, paramlist,0)
{
  if (compute_now)
    schedule(true);
}

template <class Node>
void LevelScheduler<Node>::schedule(bool force_scheduling)
{
  if (alreadyComputed() && !force_scheduling)
    return;

  const ::Tpetra::CrsGraph<int,int,Node> &graph = *input_graph_;

  int nrows = graph.getNodeNumRows();                  // graph.NumMyRows();
  int maxNonZeros = graph.getNodeMaxNumRowEntries();   // graph.MaxNumIndices();

  this->properties_.clear();              // nondependent template inheritance weirdness
  this->properties_.assign(nrows, 0);     // nondependent template inheritance weirdness
  
  if ((nrows < 2) || (maxNonZeros < 1))
  {
    this->computeNumberOfProperties(); // nondependent template inheritance weirdness
    return;
  }

  // algorithm from legacy Petra_CRS_Graph.cc 

  int numIDs=0;
  Teuchos::ArrayRCP<const int> rowView;

  if (graph.isLowerTriangular())
  {
    for (int i=0; i < nrows ; i++)
    {
      rowView = graph.getLocalRowView(i);
      numIDs = rowView.size();
  
      int depth = -1;
  
      for (int j=0; j < numIDs; j++)
      {
        int col = rowView[j];
        if ((col < i) && (this->properties_[col] > depth))     // nondependent template inheritance weirdness
	{
          depth = this->properties_[col];                      // nondependent template inheritance weirdness
	}
      }
      depth++;
  
      this->properties_[i] = depth;                            // nondependent template inheritance weirdness
    }
  }
  else if (graph.isUpperTriangular())
  {
    for (int i=nrows-1; i >= 0 ; i--)
    {
      rowView = graph.getLocalRowView(i);  
      numIDs = rowView.size();

      int depth = -1;
  
      for (int j=0; j < numIDs; j++)
      {
        int col = rowView[j];
        if ((col > i) && (this->properties_[col] > depth))             // nondependent template inheritance weirdness
	{
          depth = this->properties_[col];
	}
      }
      depth++;
  
      this->properties_[i] = depth;
    }
  }
  else
  {
    // error
    this->properties_.assign(nrows,-1);
  }

  this->operation_already_computed_ = true;

  this->computeNumberOfProperties();
}

} // namespace TPETRA

}//namespace Isorropia

#endif //HAVE_ISORROPIA_TPETRA

