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

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>


namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgRandom.hpp
//! \brief Random ordering using the Fisher-Yates (Knuth) shuffle.
//! \brief TODO: Only local permutation, could add global option.

template <typename Adapter>
class AlgRandom : public Algorithm<Adapter>
{
  private:

  const RCP<IdentifierModel<Adapter> > model;
  const RCP<Teuchos::ParameterList> pl;
  const RCP<const Teuchos::Comm<int> > comm;

  public:

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;

  AlgRandom(
    const RCP<IdentifierModel<Adapter> > &model__,
    const RCP<Teuchos::ParameterList> &pl__,
    const RCP<const Teuchos::Comm<int> > &comm__
  ) : model(model__), pl(pl__), comm(comm__)
  {
  }

  int globalOrder(const RCP<GlobalOrderingSolution<gno_t> > &/* solution */)
  {
    throw std::logic_error("AlgRandom does not yet support global ordering.");
  }

  int localOrder(const RCP<LocalOrderingSolution<lno_t> > &solution)
  {
  
    int ierr= 0;
  
    HELLO;
  
    // This is the Fisher-Yates shuffle (aka Knuth shuffle).
    // References:
    //   R.A. Fisher, F. Yates, (1948) [1938], Statistical tables for biological, agricultural and medical research (3rd ed.).
    //   R. Durstenfeld, "Algorithm 235: Random permutation", CACM, vol. 7, 1964.
    //   D.E. Knuth, "The Art of Computer Programming", volume 2, 1969.
  
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
  
    solution->setHavePerm(true);
    return ierr;
  
  }
  
};
}
#endif
