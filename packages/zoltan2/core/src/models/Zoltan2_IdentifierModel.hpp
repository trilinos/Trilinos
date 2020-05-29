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

/*! \file Zoltan2_IdentifierModel.hpp
    \brief Defines the IdentifierModel interface.
*/

#ifndef _ZOLTAN2_IDENTIFIERMODEL_HPP_
#define _ZOLTAN2_IDENTIFIERMODEL_HPP_

#include <Zoltan2_Adapter.hpp>
#include <Zoltan2_Model.hpp>
#include <Zoltan2_StridedData.hpp>

namespace Zoltan2 {

/*!  \brief IdentifierModel defines the interface for all identifier models.

    The constructor of the IdentifierModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    methods should be local methods.

    The template parameter is an InputAdapter, which is an object that
    provides a uniform interface for models to the user's input data.
*/

template <typename Adapter>
class IdentifierModel : public Model<Adapter> 
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
#endif

  /*! \brief Constructor
       \param ia  the input adapter from which to build the model
       \param env   the application environment (including problem parameters)
       \param comm  the problem communicator
       \param modelFlags   bit map of Zoltan2::IdentifierModelFlags
   */
  
  IdentifierModel(const RCP<const Adapter> &ia, 
                  const RCP<const Environment> &env,
                  const RCP<const Comm<int> > &comm, modelFlag_t &modelFlags);

  /*! Returns the number of identifiers on this process.
   */
  inline size_t getLocalNumIdentifiers() const { return gids_.size(); }

  /*! Returns the global number identifiers.
   */
  inline global_size_t getGlobalNumIdentifiers() const {
    return numGlobalIdentifiers_;
  }

  /*! Sets pointers to this process' identifier Ids and their weights.
      \param Ids will on return point to the list of the global Ids for
        each identifier on this process.
      \param wgts will on return point to a list of the weight or weights
         associated with each identifier in the Ids list. Each weight
         is represented as a StridedData object.
       \return The number of ids in the Ids list.
   */
  inline size_t getIdentifierList(ArrayView<const gno_t> &Ids,
                                  ArrayView<input_t> &wgts) const 
  {
    Ids = ArrayView<const gno_t>();
    wgts = weights_.view(0, nUserWeights_);
    size_t n = getLocalNumIdentifiers();
    if (n){
      Ids = ArrayView<const gno_t>(
                      reinterpret_cast<const gno_t*>(gids_.getRawPtr()), n);
    }
    return n;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  inline size_t getLocalNumObjects() const {return getLocalNumIdentifiers();}

  inline size_t getGlobalNumObjects() const {return getGlobalNumIdentifiers();}

private:
  gno_t numGlobalIdentifiers_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  ArrayRCP<const gno_t> gids_;
  int nUserWeights_;
  ArrayRCP<input_t> weights_;
};

////////////////////////////////////////////////////
template <typename Adapter>
  IdentifierModel<Adapter>::IdentifierModel( 
    const RCP<const Adapter> &ia,
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    modelFlag_t &/* modelFlags */):
      numGlobalIdentifiers_(), env_(env), comm_(comm),
      gids_(), nUserWeights_(0), weights_()
{
  // Get the local and global problem size
  size_t nLocalIds = ia->getLocalNumIDs();
  gno_t lsum = nLocalIds;
  reduceAll<int, gno_t>(*comm_, Teuchos::REDUCE_SUM, 1, &lsum,
    &numGlobalIdentifiers_);

  // Get the number of weights
  int tmp = ia->getNumWeightsPerID();
  // Use max number of weights over all processes as nUserWeights_
  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_MAX, 1,
      &tmp, &nUserWeights_);

  // Prepare to store views from input adapter
  // TODO:  Do we have to store these views, or can we get them on an 
  // TODO:  as-needed basis?
  Array<const scalar_t *> wgts(nUserWeights_, (const scalar_t *)NULL);
  Array<int> wgtStrides(nUserWeights_, 0);

  if (nUserWeights_ > 0){
    input_t *w = new input_t [nUserWeights_];
    weights_ = arcp<input_t>(w, 0, nUserWeights_);
  }

  const gno_t *gids=NULL;
  
  // Get the input adapter's views
  try{
    ia->getIDsView(gids);
    for (int idx=0; idx < nUserWeights_; idx++)
      ia->getWeightsView(wgts[idx], wgtStrides[idx], idx);
  }
  Z2_FORWARD_EXCEPTIONS;

  if (nLocalIds){
    gids_ = arcp(gids, 0, nLocalIds, false);

    if (nUserWeights_ > 0){
      for (int idx=0; idx < nUserWeights_; idx++){
        ArrayRCP<const scalar_t> wgtArray(wgts[idx], 0,
                                          nLocalIds*wgtStrides[idx], false);
        weights_[idx] = input_t(wgtArray, wgtStrides[idx]);
      }
    }
  }

  env_->memory("After construction of identifier model");
}

}  // namespace Zoltan2

#endif
