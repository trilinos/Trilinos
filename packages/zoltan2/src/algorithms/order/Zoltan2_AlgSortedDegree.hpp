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
#ifndef _ZOLTAN2_SORTEDDEGREE_HPP_
#define _ZOLTAN2_SORTEDDEGREE_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <algorithm>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgSortedDegree.hpp
//! \brief Order vertices by sorted (increasing) degree.
//! \brief Sorting by decreasing degree is also possible.

namespace Zoltan2{

// Comparison functions for sort. Note sort_inc is defined in some
// versions of the gcc library, so spell it 'zort_inc' here.
bool zort_inc(std::pair<size_t,size_t> a, std::pair<size_t,size_t> b)
{
  return (a.first < b.first);
}
bool zort_dec(std::pair<size_t,size_t> a, std::pair<size_t,size_t> b)
{
  return (a.first > b.first);
}

template <typename Adapter>
class AlgSortedDegree
{
  public:

    AlgSortedDegree()
    {
    }

  int order(
    const RCP<GraphModel<Adapter> > &model, 
    const RCP<OrderingSolution<typename Adapter::gid_t,
                               typename Adapter::lno_t> > &solution,
    const RCP<Teuchos::ParameterList> &pl,
    const RCP<Teuchos::Comm<int> > &comm
  ) 
  {
    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;
    typedef typename Adapter::gid_t gid_t;
    typedef typename Adapter::scalar_t scalar_t;
  
    int ierr= 0;
  
    HELLO;
  
    lno_t *perm;
    //perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());
    perm = (lno_t *) (solution->getPermutation());
    if (perm==0){
      // Throw exception
      cerr << "perm is NULL" << std::endl;
      ierr = -1;
    }
  
    // Get local graph.
    const size_t nVtx = model->getLocalNumVertices();
    ArrayView<const lno_t> edgeIds;
    ArrayView<const lno_t> offsets;
    ArrayView<StridedData<lno_t, scalar_t> > wgts;
    model->getLocalEdgeList(edgeIds, offsets, wgts);
  
    // Store degrees together with index so we can sort.
    std::vector<std::pair<size_t, size_t> >  degrees(nVtx);
    for (lno_t i=0; i<(lno_t)nVtx; i++){
      degrees[i].first  = offsets[i+1] - offsets[i];
      degrees[i].second = i;
    }
  
    // Sort degrees.
    if (1) // TODO: Check parameter for inc/dec order
      std::sort(degrees.begin(), degrees.end(), zort_inc);
    else
      std::sort(degrees.begin(), degrees.end(), zort_dec);
  
    // Copy permuted indices to perm.
    for (lno_t i=0; i<(lno_t)nVtx; i++){
      perm[i] = degrees[i].second;
    }
  
    return ierr;
  }
  
};
}
#endif
