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
//! \file Zoltan2_AlgFirstFit.hpp
//! \brief Greedy first-fit graph coloring (serial, local graph only)


namespace Zoltan2{

template <typename Adapter>
class AlgFirstFit
{
  private:
    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;
    typedef typename Adapter::scalar_t scalar_t;
  
  public:
  AlgFirstFit()
  {
  }

  void color(
    const RCP<GraphModel<Adapter> > &model,
    const RCP<ColoringSolution<typename Adapter::gid_t,
  			     typename Adapter::lno_t> > &solution,
    const RCP<Teuchos::ParameterList> &pl,
    const RCP<Teuchos::Comm<int> > &comm
  )
  {
    HELLO;
  
    // Check size of communicator: serial only.
    // TODO: Remove this test when RCM works on local graph.
    //if (comm->getSize() > 1){
    //  throw std::runtime_error("RCM currently only works in serial.");
    //}
  
    // Get local graph.
    ArrayView<const lno_t> edgeIds;
    ArrayView<const lno_t> offsets;
    ArrayView<StridedData<lno_t, scalar_t> > wgts;
  
    const size_t nVtx = model->getLocalNumVertices();
    model->getLocalEdgeList(edgeIds, offsets, wgts); 
  
#if 0
    // Debug
    cout << "Debug: Local graph from getLocalEdgeList" << endl;
    cout << "rank " << comm->getRank() << ": nVtx= " << nVtx << endl;
    cout << "rank " << comm->getRank() << ": edgeIds: " << edgeIds << endl;
    cout << "rank " << comm->getRank() << ": offsets: " << offsets << endl;
#endif
  
    // Get color array to fill
    ArrayRCP<lno_t> colors = solution->getColors();

    // TODO: First-fit coloring.
  
    return;
  }
  
};
}
#endif
#endif //INCLUDE_ZOLTAN2_EXPERIMENTAL
