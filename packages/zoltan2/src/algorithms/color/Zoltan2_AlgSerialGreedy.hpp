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
#ifdef INCLUDE_ZOLTAN2_EXPERIMENTAL
#ifndef _ZOLTAN2_ALGFIRSTFIT_HPP_
#define _ZOLTAN2_ALGFIRSTFIT_HPP_

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
    typedef int color_t; // TODO get from adapter
  
  public:
  AlgSerialGreedy()
  {
  }

  void color(
    const RCP<GraphModel<Adapter> > &model,
    const RCP<ColoringSolution<Adapter> > &solution,
    const RCP<Teuchos::ParameterList> &pl,
    const RCP<Teuchos::Comm<int> > &comm
  )
  {
    HELLO;
  
    // Only color local graph. Global coloring is supported in Zoltan (not Zoltan2).
    // Get local graph.
    ArrayView<const lno_t> edgeIds;
    ArrayView<const lno_t> offsets;
    ArrayView<StridedData<lno_t, scalar_t> > wgts; // Not used; needed by getLocalEdgeList
  
    const lno_t nVtx = model->getLocalNumVertices();
    model->getLocalEdgeList(edgeIds, offsets, wgts); // Don't need wgts
  
#if 0
    // Debug
    cout << "Debug: Local graph from getLocalEdgeList" << endl;
    cout << "rank " << comm->getRank() << ": nVtx= " << nVtx << endl;
    cout << "rank " << comm->getRank() << ": edgeIds: " << edgeIds << endl;
    cout << "rank " << comm->getRank() << ": offsets: " << offsets << endl;
#endif
  
    // Get color array to fill.
    // TODO: Allow user to input an old coloring.
    ArrayRCP<color_t> colors = solution->getColorsRCP();

    // First-fit greedy coloring.
    // Use natural order for now. 
    // TODO: Support better orderings (e.g., Smallest-Last)
    const color_t maxColorGuess = 256; // for array allocation
    color_t maxColor = 0;
 
    // array of size #colors: forbidden[i]=v means color[v]=i so i is forbidden
    Teuchos::Array<color_t> forbidden(maxColorGuess, -1);

    for (lno_t i=0; i<nVtx; i++){
      lno_t v=i; // TODO: Use ordering here.
      for (lno_t j=offsets[v]; j<offsets[v+1]; j++){
        lno_t nbor = edgeIds[j];
        if (colors[nbor] > 0){
          // TODO: Check if we need reallocate forbidden array.
          forbidden[colors[nbor]] = v;
        }
      }
      // Pick first (smallest) available color > 0
      color_t c=1;
      while (forbidden[c]==v) c++;
      colors[v] = c;
      if (c > maxColor){
        maxColor = c;
      }
    }
  
    // Set numColors in solution
    solution->setNumColors(maxColor); // solution->numColors_ = maxColor;

    return;
  }
  
};
}
#endif
#endif //INCLUDE_ZOLTAN2_EXPERIMENTAL
