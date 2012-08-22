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
#ifndef _ZOLTAN2_INCDEGREE_HPP_
#define _ZOLTAN2_INCDEGREE_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <algorithm>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgIncDegree.hpp
//! \brief Order vertices by increasing degree

// Comparison function for sort.
bool comp(std::pair<lno_t,lno_t> a, std::pair<lno_t,lno_t> b)
{
  return (a.first < b.first);
}

namespace Zoltan2{

template <typename Adapter>
int AlgIncDegree(
  const RCP<GraphModel<Adapter> > &model, 
  const RCP<OrderingSolution<typename Adapter::gid_t,
                             typename Adapter::lno_t> > &solution,
  const RCP<Teuchos::ParameterList> &pl,
  const RCP<Teuchos::Comm<int> > &comm
) 
{
#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL

  int ierr= 0;
  Z2_THROW_EXPERIMENTAL("Zoltan2 IncreasingDegree ordering is strictly "
                        "experimental software "
                        "while it is being developed and tested.")
  return ierr;

#else //INCLUDE_ZOLTAN2_EXPERIMENTAL

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::scalar_t scalar_t;

  int ierr= 0;

  HELLO;

  lno_t *perm;
  perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());

  // Get local graph.
  ArrayView<const lno_t> edgeIds;
  ArrayView<const lno_t> offsets;
  ArrayView<const scalar_t> wgts;
  model->getLocalEdgeList(edgeIds, offsets, wgts);

  // Store degrees together with index so we can sort.
  std::vector<std::pair<lno_t, lno_t> >  degrees(nVtx);
  for (lno_t i=0; i<nVtx; i++){
    degrees[i].first  = offsets[i+1] - offsets[i];
    degrees[i].second = i;
  }

  // Sort degrees.
  std::sort(degrees.begin(), degrees.end(), comp);

  // Copy permuted indices to perm.
  for (lno_t i=0; i<nVtx; i++){
    perm[i] = degrees[i].second;
  }

  return ierr;
#endif // INCLUDE_ZOLTAN2_EXPERIMENTAL
}

}
#endif
