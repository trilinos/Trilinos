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
#ifndef _ZOLTAN2_ALGRANDOM_HPP_
#define _ZOLTAN2_ALGRANDOM_HPP_

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgRandom.hpp
//! \brief Random ordering using the Knuth shuffle.
//! \brief TODO: Only local permutation, should add global option.


namespace Zoltan2{

template <typename Adapter>
int AlgRandom(
  const RCP<IdentifierModel<Adapter> > &model, 
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

#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL

  Z2_THROW_EXPERIMENTAL("Zoltan2 random ordering is strictly "
                        "experimental software "
                        "while it is being developed and tested.")
  return ierr;

#else //INCLUDE_ZOLTAN2_EXPERIMENTAL

  HELLO;

  // This is the classic Knuth shuffle, also known as Fisher-Yates shuffle.
  // References:
  //   D.E. Knuth, "The Art of Computer Programming", volume 2, 1969.
  //   R. Durstenfeld, "Algorithm 235: Random permutation", CACM, vol. 7, 1964.

  // Start with the identity permutation.
  const size_t n = model->getLocalNumIdentifiers();
  lno_t *perm;
  perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());
  if (perm){
    for (size_t i=0; i<n; i++){
      perm[i] = i;
    }
  }
  else
    // throw exception?
    ierr = -1;

  // Swap random pairs of indices in perm.
  lno_t j, temp;
  for (lno_t i=n-1; i>0; i--){
    // Choose j randomly in [0,i]
    j = rand() % (i+1);
    // Swap (perm[i], perm[j])
    temp = perm[i];
    perm[i] = perm[j];
    perm[j] = temp;
  }

  return ierr;

#endif //INCLUDE_ZOLTAN2_EXPERIMENTAL

}

}
#endif
