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

