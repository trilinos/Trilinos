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
#ifndef _ZOLTAN2_ALGSPECTRAL_HPP_
#define _ZOLTAN2_ALGSPECTRAL_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgSpectral.hpp
//! \brief Spectral ordering of a graph (local or global).
//! \brief Sorts the Fiedler vector of the graph Laplacian.


namespace Zoltan2{

template <typename Adapter>
int AlgSpectral(
  const RCP<GraphModel<Adapter> > &model, 
  //const RCP<Adapter> &matrixadapter, // Hack: Use matrix adapter directly

  // TO DO - update this algorithm to have the proper formatting like
  // the others.
  const RCP<LocalOrderingSolution<typename Adapter::lno_t> >
    &solution,
  const RCP<Teuchos::ParameterList> &pl,
  const RCP<const Teuchos::Comm<int> > &comm
) 
{
#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL

  Z2_THROW_EXPERIMENTAL("Zoltan2 Spectral ordering is strictly "
                        "experimental software "
                        "while it is being developed and tested.")

#else //INCLUDE_ZOLTAN2_EXPERIMENTAL

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::scalar_t scalar_t;

  int ierr= 0;

  HELLO;

// TODO: Check params to do local or global ordering.
  bool localOrder = true;

  const size_t nVtx = model->getLocalNumVertices();
  lno_t *perm = solution->getPermutationView();
  for (lno_t i=0; i<nVtx; i++){
    perm[i] = -1;
  }

  // Get local graph.
  ArrayView<const gno_t> edgeIds;
  ArrayView<const offset_t> offsets;
  ArrayView<StridedData<lno_t, scalar_t> > wgts;

  model->getEdgeList(edgeIds, offsets, wgts);

  //cout << "Debug: Local graph from getLocalEdgeList" << endl;
  //cout << "edgeIds: " << edgeIds << endl;
  //cout << "offsets: " << offsets << endl;

  // Form the graph Laplacian: L = D-A
  // Create a new Tpetra matrix, but use views of existing data when possible.
  // TODO

  // TODO: Find smallest eigenvalues using Anasazi

  // TODO: Sort Fiedler vector.

  solution->setHavePerm(true);
  return ierr;
#endif // INCLUDE_ZOLTAN2_EXPERIMENTAL
}

}
#endif
