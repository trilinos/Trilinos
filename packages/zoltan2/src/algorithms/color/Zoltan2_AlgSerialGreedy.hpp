// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef _ZOLTAN2_ALGSERIALGREEDY_HPP_
#define _ZOLTAN2_ALGSERIALGREEDY_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_ColoringSolution.hpp>

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgSerialGreedy.hpp
//! \brief Serial greedy first-fit graph coloring (local graph only)

namespace Zoltan2{

template <typename Adapter>
class AlgSerialGreedy
{
  private:
    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;
    typedef typename Adapter::scalar_t scalar_t;
    // Class member variables
    RCP<GraphModel<typename Adapter::base_adapter_t> > model_;
    RCP<Teuchos::Comm<int> > comm_;
  
  public:
  AlgSerialGreedy(
    const RCP<GraphModel<typename Adapter::base_adapter_t> > &model,
    const RCP<Teuchos::Comm<int> > &comm
  ) : model_(model), comm_(comm)
  {
  }

  void color(
    const RCP<ColoringSolution<Adapter> > &solution,
    const RCP<Teuchos::ParameterList> &pl
  )
  {
    HELLO;
  
    // Color local graph. Global coloring is supported in Zoltan (not Zoltan2).
    // Get local graph.
    ArrayView<const lno_t> edgeIds;
    ArrayView<const lno_t> offsets;
    ArrayView<StridedData<lno_t, scalar_t> > wgts; // Not used; needed by getLocalEdgeList
  
    const lno_t nVtx = model_->getLocalNumVertices();
    model_->getLocalEdgeList(edgeIds, offsets, wgts); // Don't need wgts
  
#if 0
    // Debug
    cout << "Debug: Local graph from getLocalEdgeList" << endl;
    cout << "rank " << comm_->getRank() << ": nVtx= " << nVtx << endl;
    cout << "rank " << comm_->getRank() << ": edgeIds: " << edgeIds << endl;
    cout << "rank " << comm_->getRank() << ": offsets: " << offsets << endl;
#endif
  
    // Get color array to fill.
    // TODO: Allow user to input an old coloring.
    ArrayRCP<int> colors = solution->getColorsRCP();
    for (lno_t i=0; i<nVtx; i++){
      colors[i] = 0;
    }

    // Find max degree, since (max degree)+1 is an upper bound.
    lno_t maxDegree = 0; 
    for (lno_t i=0; i<nVtx; i++){
      if (offsets[i+1]-offsets[i] > maxDegree)
        maxDegree = offsets[i+1]-offsets[i];
    }

    // First-fit greedy coloring.
    // Use natural order for now. 
    // TODO: Support better orderings (e.g., Smallest-Last)
    int maxColor = 0;
 
    // array of size #colors: forbidden[i]=v means color[v]=i so i is forbidden
    Teuchos::Array<int> forbidden(maxDegree+2, 0);

    for (lno_t i=0; i<nVtx; i++){
      //std::cout << "Debug: i= " << i << std::endl;
      lno_t v=i; // TODO: Use ordering here.
      for (lno_t j=offsets[v]; j<offsets[v+1]; j++){
        lno_t nbor = edgeIds[j];
        //std::cout << "Debug: nbor= " << nbor << ", color= " << colors[nbor] << std::endl;
        if (colors[nbor] > 0){
          // Neighbors' colors are forbidden
          forbidden[colors[nbor]] = v;
        }
      }
      // Pick first (smallest) available color > 0
      for (int c=1; c < forbidden.length(); c++){
        if (forbidden[c] != v){ 
          colors[v] = c;
          break;
        }
      }
      if (colors[v]==0) colors[v]=1; // Corner case for first vertex
      //std::cout << "Debug: colors[i]= " << colors[v] << std::endl;
      if (colors[v] > maxColor){
        maxColor = colors[v];
      }
    }
  
    // Set numColors in solution
    solution->setNumColors(maxColor); // solution->numColors_ = maxColor;

    return;
  }
  
};
}
#endif
