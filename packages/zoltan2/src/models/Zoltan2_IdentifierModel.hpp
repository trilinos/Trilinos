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

#include <Zoltan2_Model.hpp>
#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_VectorInput.hpp>
#include <Zoltan2_GraphInput.hpp>
#include <Zoltan2_IdentifierInput.hpp>
#include <Zoltan2_CoordinateInput.hpp>
#include <Zoltan2_StridedData.hpp>

namespace Zoltan2 {

/*!  \brief IdentifierModel defines the interface for all identifier models.

    The constructor of the IdentifierModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    methods should be local methods.

    The template parameter is an InputAdapter, which is an object that
    provides a uniform interface for models to the user's input data.

    Explicit instantiations exist for:
      \li MatrixAdapter
      \li IdentifierAdapter
      \li VectorAdapter
      \li CoordinateAdapter

    \todo Add instantiations for GraphAdapter, MeshInput
*/

template <typename Adapter>
class IdentifierModel : public Model<Adapter> 
{
private:

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
#endif
  
  IdentifierModel(){
    throw std::logic_error("a specific instantiation should be used");
  }

  ////////////////////////////////////////////////////
  // The IdentifierModel interface.
  ////////////////////////////////////////////////////

  /*! \brief Returns the number identifiers on this process.
   */
  size_t getLocalNumIdentifiers() const { return 0; }

  /*! \brief Returns the global number identifiers.
   */
  global_size_t getGlobalNumIdentifiers() const { return 0; }

  /*! \brief Returns the dimension (0 or greater) of identifier weights.
   */
  int getIdentifierWeightDim() const { return 0; }

  /*! \brief Sets pointers to this process' identifier Ids and their weights.

      \param Ids will on return point to the list of the global Ids for
        each identifier on this process.
      \param wgts will on return point to a list of the weight or weights
         associated with each identifier in the Ids list. Each weight
         is represented as a StridedData object.

       \return The number of ids in the Ids list.
   */

  size_t getIdentifierList(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &wgts) const {return 0;}

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumIdentifiers();
  }

  size_t getGlobalNumObjects() const
  {
    return getGlobalNumIdentifiers();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const { return ; }
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS

////////////////////////////////////////////////////////////////
// Identifier model derived from IdentifierAdapter.
////////////////////////////////////////////////////////////////

template <typename User>
class IdentifierModel<IdentifierAdapter<User> > : public Model<IdentifierAdapter<User> >
{
public:

  typedef typename IdentifierAdapter<User>::scalar_t  scalar_t;
  typedef typename IdentifierAdapter<User>::gno_t     gno_t;
  typedef typename IdentifierAdapter<User>::lno_t     lno_t;
  typedef typename IdentifierAdapter<User>::gid_t     gid_t;
  typedef IdentifierMap<User> idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  /*! \brief Constructor

       \param ia  the input adapter from which to build the model
       \param env   the application environment (including problem parameters)
       \param comm  the problem communicator
       \param modelFlags   bit map of Zoltan2::IdentifierModelFlags
   */
  
  IdentifierModel( const IdentifierAdapter<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &modelFlags);

  /*! Returns the number identifiers on this process.
   */
  size_t getLocalNumIdentifiers() const { return gids_.size(); }

  /*! Returns the global number identifiers.
   */
  global_size_t getGlobalNumIdentifiers() const {return numGlobalIdentifiers_;}

  /*! Returns the dimension (0 or greater) of identifier weights.
   */
  int getIdentifierWeightDim() const { return this->getNumWeights(); }

  /*! Sets pointers to this process' identifier Ids and their weights.
      \param Ids will on return point to the list of the global Ids for
        each identifier on this process.
      \param wgts will on return point to a list of the weight or weights
         associated with each identifier in the Ids list. Each weight
         is represented as a StridedData object.
         
       \return The number of ids in the Ids list.
   */

  size_t getIdentifierList(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &wgts) const 
  {
    Ids = ArrayView<const gno_t>();
    wgts = weights_.view(0, userWeightDim_);

    size_t n = getLocalNumIdentifiers();
    if (n){
      if (gnosAreGids_)
        Ids = gids_(0, n);
      else
        Ids = gnosConst_(0, n);
    }

    return n;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumIdentifiers();
  }

  size_t getGlobalNumObjects() const
  {
    return getGlobalNumIdentifiers();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const 
  { 
    ArrayView<input_t> weights;
    getIdentifierList(gnos, weights);
  }

private:

  bool gnosAreGids_;
  gno_t numGlobalIdentifiers_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  ArrayRCP<const gid_t> gids_;
  int userWeightDim_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;
};

template <typename User>
  IdentifierModel<IdentifierAdapter<User> >::IdentifierModel( 
    const IdentifierAdapter<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags):
      gnosAreGids_(false), numGlobalIdentifiers_(), env_(env), comm_(comm),
      gids_(), userWeightDim_(0), weights_(), gnos_(), gnosConst_()
{
  userWeightDim_ = ia->getNumberOfWeights();
  size_t nLocalIds = ia->getLocalNumberOfIdentifiers();

  Model<IdentifierAdapter<User> >::maxCount(*comm, userWeightDim_);

  Array<const scalar_t *> wgts(userWeightDim_, (const scalar_t *)NULL);
  Array<int> wgtStrides(userWeightDim_, 0);
  Array<lno_t> weightArrayLengths(userWeightDim_, 0);

  if (userWeightDim_ > 0){
    input_t *w = new input_t [userWeightDim_];
    weights_ = arcp<input_t>(w, 0, userWeightDim_);
  }

  const gid_t *gids=NULL;

  try{
    ia->getIdentifierList(gids);
    for (int dim=0; dim < userWeightDim_; dim++)
      ia->getIdentifierWeights(dim, wgts[dim], wgtStrides[dim]);
  }
  Z2_FORWARD_EXCEPTIONS;

  if (nLocalIds){
    gids_ = arcp(gids, 0, nLocalIds, false);

    if (userWeightDim_ > 0){
      for (int i=0; i < userWeightDim_; i++){
        if (wgts[i] != NULL){
          ArrayRCP<const scalar_t> wgtArray(
            wgts[i], 0, nLocalIds*wgtStrides[i], false);
          weights_[i] = input_t(wgtArray, wgtStrides[i]);
          weightArrayLengths[i] = nLocalIds;
        }
      }
    }
  }

  this->setWeightArrayLengths(weightArrayLengths, *comm_);

  RCP<const idmap_t> idMap;

  try{
    if (modelFlags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE))
      idMap = rcp(new idmap_t(env_, comm_, gids_, true));
    else
      idMap = rcp(new idmap_t(env_, comm_, gids_, false));
  }
  Z2_FORWARD_EXCEPTIONS;

  gnosAreGids_ = idMap->gnosAreGids();

  this->setIdentifierMap(idMap);

  gno_t lsum = nLocalIds;
  reduceAll<int, gno_t>(*comm_, Teuchos::REDUCE_SUM, 1, &lsum,
    &numGlobalIdentifiers_);

  if (!gnosAreGids_ && nLocalIds>0){
    gno_t *tmpGno = new gno_t [nLocalIds];
    env_->localMemoryAssertion(__FILE__, __LINE__, nLocalIds, tmpGno);
    gnos_ = arcp(tmpGno, 0, nLocalIds);

    try{
      ArrayRCP<gid_t> gidsNonConst = arcp_const_cast<gid_t>(gids_);
      idMap->gidTranslate( gidsNonConst(0,nLocalIds),  gnos_(0,nLocalIds),
        TRANSLATE_APP_TO_LIB);
    }
    Z2_FORWARD_EXCEPTIONS;
  }

  gnosConst_ = arcp_const_cast<const gno_t>(gnos_);

  env_->memory("After construction of identifier model");
}

////////////////////////////////////////////////////////////////
// Identifier model derived from CoordinateAdapter.
////////////////////////////////////////////////////////////////

template <typename User>
  class IdentifierModel<CoordinateAdapter<User> > : 
    public Model<CoordinateAdapter<User> >
{
public:

  typedef typename CoordinateAdapter<User>::scalar_t  scalar_t;
  typedef typename CoordinateAdapter<User>::gno_t     gno_t;
  typedef typename CoordinateAdapter<User>::lno_t     lno_t;
  typedef typename CoordinateAdapter<User>::gid_t     gid_t;
  typedef IdentifierMap<User> idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  
  IdentifierModel( const CoordinateAdapter<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &modelFlags);

  /*! Returns the number identifiers on this process.
   */
  size_t getLocalNumIdentifiers() const { return gids_.size(); }

  /*! Returns the global number identifiers.
   */
  global_size_t getGlobalNumIdentifiers() const {return numGlobalIdentifiers_;}

  /*! Returns the dimension (0 or greater) of identifier weights.
   */
  int getIdentifierWeightDim() const { return this->getNumWeights();}

  /*! Sets pointers to this process' identifier Ids and their weights.
      \param Ids will on return point to the list of the global Ids for
        each identifier on this process.
      \param wgts will on return point to a list of the weight or weights
         associated with each identifier in the Ids list. Each weight
         is represented as a StridedData object.

       \return The number of ids in the Ids list.
   */

  size_t getIdentifierList(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &wgts) const            
  {
    Ids = ArrayView<const gno_t>();
    wgts = weights_.view(0, userWeightDim_);

    size_t n = getLocalNumIdentifiers();
    if (n){
      if (gnosAreGids_)
        Ids = gids_(0, n);
      else
        Ids = gnosConst_(0, n);
    }

    return n;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumIdentifiers();
  }

  size_t getGlobalNumObjects() const
  {
    return getGlobalNumIdentifiers();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const 
  { 
    ArrayView<input_t> weights;
    getIdentifierList(gnos, weights);
  }

private:

  bool gnosAreGids_;
  gno_t numGlobalIdentifiers_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  ArrayRCP<const gid_t> gids_;
  int userWeightDim_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;
};

template <typename User>
  IdentifierModel<CoordinateAdapter<User> >::IdentifierModel( 
    const CoordinateAdapter<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &modelFlags):
      gnosAreGids_(false), numGlobalIdentifiers_(), env_(env), comm_(comm),
      gids_(), userWeightDim_(0), weights_(), gnos_(), gnosConst_()
{
  /////////////////////////////////////////
  // Get global IDs.

  size_t nLocalIds = ia->getLocalNumberOfCoordinates();
  const gid_t *gids=NULL;

  if (nLocalIds > 0){
    try{
      size_t coordListSize; 
      const scalar_t *coords;
      int stride;
      coordListSize = ia->getCoordinates(0, gids, coords, stride);
    }
    Z2_FORWARD_EXCEPTIONS;

    gids_ = arcp(gids, 0, nLocalIds, false);
  }

  RCP<const idmap_t> idMap;

  try{
    if (modelFlags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE) )
      idMap = rcp(new idmap_t(env_, comm_, gids_, true));
    else
      idMap = rcp(new idmap_t(env_, comm_, gids_, false));
  }
  Z2_FORWARD_EXCEPTIONS;

  numGlobalIdentifiers_ = idMap->getGlobalNumberOfIds();
  gnosAreGids_ = idMap->gnosAreGids();
  this->setIdentifierMap(idMap);   // Base Model method

  /////////////////////////////////////////
  // Get weights.

  userWeightDim_ = ia->getNumberOfWeights();
  Array<lno_t> weightListSizes(userWeightDim_, 0);

  Model<CoordinateAdapter<User> >::maxCount(*comm, userWeightDim_);

  if (userWeightDim_ > 0){
    input_t *weightObj = new input_t [userWeightDim_];
    weights_ = arcp(weightObj, 0, userWeightDim_);

    if (nLocalIds > 0){
      const scalar_t *wgts=NULL;
      int stride = 0;
      size_t wgtListSize;

      for (int wdim=0; wdim < userWeightDim_; wdim++){
        wgtListSize = ia->getCoordinateWeights(wdim, wgts, stride);

        if (wgtListSize > 0){  // non-uniform weights
          ArrayRCP<const scalar_t> wgtArray(wgts, 0, wgtListSize, false);
          weightObj[wdim] = StridedData<lno_t, scalar_t>(wgtArray, stride);
          weightListSizes[wdim] = wgtListSize;
        }
      }
    }
  }

  this->setWeightArrayLengths(weightListSizes, *comm_);

  /////////////////////////////////////////////
  // Get internal global numbers if necessary.

  if (!gnosAreGids_ && nLocalIds>0){
    gno_t *tmpGno = new gno_t [nLocalIds];
    env_->localMemoryAssertion(__FILE__, __LINE__, nLocalIds, tmpGno);
    gnos_ = arcp(tmpGno, 0, gids_.size());

    try{
     ArrayRCP<gid_t> gidsNonConst = arcp_const_cast<gid_t>(gids_);
     idMap->gidTranslate(gidsNonConst(0, nLocalIds), gnos_(0, nLocalIds), 
       TRANSLATE_APP_TO_LIB);
    }
    Z2_FORWARD_EXCEPTIONS;
  }

  gnosConst_ = arcp_const_cast<const gno_t>(gnos_);
  env_->memory("After construction of identifier model");
}


////////////////////////////////////////////////////////////////
// Identifier model derived from MatrixAdapter.
////////////////////////////////////////////////////////////////

template <typename User>
class IdentifierModel<MatrixAdapter<User> > : public Model<MatrixAdapter<User> >
{
public:

  typedef typename MatrixAdapter<User>::scalar_t  scalar_t;
  typedef typename MatrixAdapter<User>::gno_t     gno_t;
  typedef typename MatrixAdapter<User>::lno_t     lno_t;
  typedef typename MatrixAdapter<User>::gid_t     gid_t;
  typedef IdentifierMap<User> idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  
  IdentifierModel( const MatrixAdapter<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &modelFlags);

  /*! Returns the number identifiers on this process.
   */
  size_t getLocalNumIdentifiers() const { return gids_.size(); }

  /*! Returns the global number identifiers.
   */
  global_size_t getGlobalNumIdentifiers() const {return numGlobalIdentifiers_;}

  /*! Returns the dimension (0 or greater) of identifier weights.
   *    Weights are not yet implemented in MatrixAdapter.
   */
  int getIdentifierWeightDim() const { return 0; }

  /*! Sets pointers to this process' identifier Ids and their weights.
      \param Ids will on return point to the list of the global Ids for
        each identifier on this process.
      \param wgts will on return point to a list of the weight or weights
         associated with each identifier in the Ids list. Each weight
         is represented as a StridedData object.

       \return The number of ids in the Ids list.
   */

  size_t getIdentifierList(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &wgts) const            
  {
    Ids = ArrayView<const gno_t>();
    wgts = weights_.view(0, userWeightDim_);

    size_t n = getLocalNumIdentifiers();
    if (n){
      if (gnosAreGids_)
        Ids = gids_(0, n);
      else
        Ids = gnosConst_(0, n);
    }

    return n;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumIdentifiers();
  }

  size_t getGlobalNumObjects() const
  {
    return getGlobalNumIdentifiers();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const 
  { 
    ArrayView<input_t> weights;
    getIdentifierList(gnos, weights);
  }

private:

  bool gnosAreGids_;
  gno_t numGlobalIdentifiers_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  ArrayRCP<const gid_t> gids_;
  int userWeightDim_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;
};

  
template <typename User>
  IdentifierModel<MatrixAdapter<User> >::IdentifierModel( 
    const MatrixAdapter<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &modelFlags):
      gnosAreGids_(false), numGlobalIdentifiers_(), env_(env), comm_(comm),
      gids_(), userWeightDim_(0), weights_(), gnos_(), gnosConst_()
{
  size_t nLocalIds;
  const gid_t *gids;
  const gid_t *colIds;
  const lno_t *offsets;

  try{
    nLocalIds = ia->getRowListView(gids, offsets, colIds);
  }
  Z2_FORWARD_EXCEPTIONS;

  if (nLocalIds){
    gids_ = arcp(gids, 0, nLocalIds, false);
  }

  RCP<const idmap_t> idMap;

  try{
    if (modelFlags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE) )
      idMap = rcp(new idmap_t(env_, comm_, gids_, true));
    else
      idMap = rcp(new idmap_t(env_, comm_, gids_, false));
  }
  Z2_FORWARD_EXCEPTIONS;

  numGlobalIdentifiers_ = idMap->getGlobalNumberOfIds();
  gnosAreGids_ = idMap->gnosAreGids();

  this->setIdentifierMap(idMap);   // Base Model methods
  Array<lno_t> weightListSizes;
  this->setWeightArrayLengths(weightListSizes, *comm_);

  if (!gnosAreGids_ && nLocalIds>0){
    gno_t *tmpGno = new gno_t [nLocalIds];
    env_->localMemoryAssertion(__FILE__, __LINE__, nLocalIds, tmpGno);
    gnos_ = arcp(tmpGno, 0, gids_.size());

    try{
     ArrayRCP<gid_t> gidsNonConst = arcp_const_cast<gid_t>(gids_);
      idMap->gidTranslate(gidsNonConst(0, nLocalIds),
        gnos_(0, nLocalIds), TRANSLATE_APP_TO_LIB);
    }
    Z2_FORWARD_EXCEPTIONS;
  }

  gnosConst_ = arcp_const_cast<const gno_t>(gnos_);
  env_->memory("After construction of identifier model");
}

////////////////////////////////////////////////////////////////
// Identifier model derived from VectorAdapter.
////////////////////////////////////////////////////////////////

template <typename User>
class IdentifierModel<VectorAdapter<User> > : public Model<VectorAdapter<User> >
{
public:

  typedef typename VectorAdapter<User>::scalar_t  scalar_t;
  typedef typename VectorAdapter<User>::gno_t     gno_t;
  typedef typename VectorAdapter<User>::lno_t     lno_t;
  typedef typename VectorAdapter<User>::gid_t     gid_t;
  typedef IdentifierMap<User> idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  
  IdentifierModel( const VectorAdapter<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &modelFlags);

  /*! Returns the number identifiers on this process.
   */
  size_t getLocalNumIdentifiers() const { return gids_.size(); }

  /*! Returns the global number identifiers.
   */
  global_size_t getGlobalNumIdentifiers() const {return numGlobalIdentifiers_;}

  /*! Returns the dimension (0 or greater) of identifier weights.
   *    Weights are not yet implemented in VectorAdapter.
   */
  int getIdentifierWeightDim() const { return 0; }

  /*! Sets pointers to this process' identifier Ids and their weights.
      \param Ids will on return point to the list of the global Ids for
        each identifier on this process.
      \param wgts will on return point to a list of the weight or weights
         associated with each identifier in the Ids list. Each weight
         is represented as a StridedData object.

       \return The number of ids in the Ids list.
   */

  size_t getIdentifierList(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &wgts) const            
  {
    Ids = ArrayView<const gno_t>();
    wgts = weights_.view(0, userWeightDim_);

    size_t n = getLocalNumIdentifiers();
    if (n){
      if (gnosAreGids_)
        Ids = gids_(0, n);
      else
        Ids = gnosConst_(0, n);
    }

    return n;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumIdentifiers();
  }

  size_t getGlobalNumObjects() const
  {
    return getGlobalNumIdentifiers();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const 
  { 
    ArrayView<input_t> weights;
    getIdentifierList(gnos, weights);
  }

private:

  bool gnosAreGids_;
  gno_t numGlobalIdentifiers_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  ArrayRCP<const gid_t> gids_;
  int userWeightDim_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;
};

  
template <typename User>
  IdentifierModel<VectorAdapter<User> >::IdentifierModel( 
    const VectorAdapter<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &modelFlags):
      gnosAreGids_(false), numGlobalIdentifiers_(), env_(env), comm_(comm),
      gids_(), userWeightDim_(0), weights_(), gnos_(), gnosConst_()
{
  size_t nLocalIds;
  const gid_t *gids;
  const scalar_t *elements;
  int stride;

  try{
    nLocalIds = ia->getVector(gids, elements, stride);
  }
  Z2_FORWARD_EXCEPTIONS;

  if (nLocalIds){
    gids_ = arcp(gids, 0, nLocalIds, false);
  }

  RCP<const idmap_t> idMap;

  try{
    if (modelFlags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE) )
      idMap = rcp(new idmap_t(env_, comm_, gids_, true));
    else
      idMap = rcp(new idmap_t(env_, comm_, gids_, false));
  }
  Z2_FORWARD_EXCEPTIONS;

  numGlobalIdentifiers_ = idMap->getGlobalNumberOfIds();
  gnosAreGids_ = idMap->gnosAreGids();

  this->setIdentifierMap(idMap);   // Base Model methods
  Array<lno_t> weightListSizes;
  this->setWeightArrayLengths(weightListSizes, *comm_);

  if (!gnosAreGids_ && nLocalIds>0){
    gno_t *tmpGno = new gno_t [nLocalIds];
    env_->localMemoryAssertion(__FILE__, __LINE__, nLocalIds, tmpGno);
    gnos_ = arcp(tmpGno, 0, gids_.size());

    try{
     ArrayRCP<gid_t> gidsNonConst = arcp_const_cast<gid_t>(gids_);
      idMap->gidTranslate(gidsNonConst(0, nLocalIds),
        gnos_(0, nLocalIds), TRANSLATE_APP_TO_LIB);
    }
    Z2_FORWARD_EXCEPTIONS;
  }

  gnosConst_ = arcp_const_cast<const gno_t>(gnos_);
  env_->memory("After construction of identifier model");
}

////////////////////////////////////////////////////////////////
// Identifier model derived from GraphAdapter.
////////////////////////////////////////////////////////////////

template <typename User>
class IdentifierModel<GraphAdapter<User> > : public Model<GraphAdapter<User> >
{
public:

  typedef typename GraphAdapter<User>::scalar_t  scalar_t;
  typedef typename GraphAdapter<User>::gno_t     gno_t;
  typedef typename GraphAdapter<User>::lno_t     lno_t;
  typedef typename GraphAdapter<User>::gid_t     gid_t;
  typedef IdentifierMap<User> idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  /*! \brief Constructor

       \param ia  the input adapter from which to build the model
       \param env   the application environment (including problem parameters)
       \param comm  the problem communicator
       \param modelFlags   bit map of Zoltan2::IdentifierModelFlags
   */
  
  IdentifierModel( const GraphAdapter<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &modelFlags);

  /*! Returns the number identifiers on this process.
   */
  size_t getLocalNumIdentifiers() const { return gids_.size(); }

  /*! Returns the global number identifiers.
   */
  global_size_t getGlobalNumIdentifiers() const {return numGlobalIdentifiers_;}

  /*! Returns the dimension (0 or greater) of identifier weights.
   */
  int getIdentifierWeightDim() const { return this->getNumWeights(); }

  /*! Sets pointers to this process' identifier Ids and their weights.
      \param Ids will on return point to the list of the global Ids for
        each identifier on this process.
      \param wgts will on return point to a list of the weight or weights
         associated with each identifier in the Ids list. Each weight
         is represented as a StridedData object.
         
       \return The number of ids in the Ids list.
   */

  size_t getIdentifierList(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &wgts) const 
  {
    Ids = ArrayView<const gno_t>();
    wgts = weights_.view(0, userWeightDim_);

    size_t n = getLocalNumIdentifiers();
    if (n){
      if (gnosAreGids_)
        Ids = gids_(0, n);
      else
        Ids = gnosConst_(0, n);
    }

    return n;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumIdentifiers();
  }

  size_t getGlobalNumObjects() const
  {
    return getGlobalNumIdentifiers();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const 
  { 
    ArrayView<input_t> weights;
    getIdentifierList(gnos, weights);
  }

private:

  bool gnosAreGids_;
  gno_t numGlobalIdentifiers_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  ArrayRCP<const gid_t> gids_;
  int userWeightDim_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;
};

template <typename User>
  IdentifierModel<GraphAdapter<User> >::IdentifierModel( 
    const GraphAdapter<User> *ia,
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm,
    modelFlag_t &modelFlags):
      gnosAreGids_(false), numGlobalIdentifiers_(), env_(env), comm_(comm),
      gids_(), userWeightDim_(0), weights_(), gnos_(), gnosConst_()
{
  userWeightDim_ = ia->getVertexWeightDimension();
  size_t nLocalIds = ia->getLocalNumberOfVertices();

  Model<GraphAdapter<User> >::maxCount(*comm, userWeightDim_);

  Array<const scalar_t *> wgts(userWeightDim_, NULL);
  Array<int> wgtStrides(userWeightDim_, 0);
  Array<lno_t> weightArrayLengths(userWeightDim_, 0);

  if (userWeightDim_ > 0){
    input_t *w = new input_t [userWeightDim_];
    weights_ = arcp<input_t>(w, 0, userWeightDim_);
  }

  const gid_t *gids=NULL;

  try{
    const lno_t *offsets;
    const gid_t *nbors;
    ia->getVertexListView(gids, offsets, nbors);
    for (int dim=0; dim < userWeightDim_; dim++)
      ia->getVertexWeights(dim, wgts[dim], wgtStrides[dim]);
  }
  Z2_FORWARD_EXCEPTIONS;

  if (nLocalIds){
    gids_ = arcp(gids, 0, nLocalIds, false);

    if (userWeightDim_ > 0){
      for (int i=0; i < userWeightDim_; i++){
        if (wgts[i] != NULL){
          ArrayRCP<const scalar_t> wgtArray(
            wgts[i], 0, nLocalIds*wgtStrides[i], false);
          weights_[i] = input_t(wgtArray, wgtStrides[i]);
          weightArrayLengths[i] = nLocalIds;
        }
      }
    }
  }

  this->setWeightArrayLengths(weightArrayLengths, *comm_);

  RCP<const idmap_t> idMap;

  try{
    if (modelFlags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE))
      idMap = rcp(new idmap_t(env_, comm_, gids_, true));
    else
      idMap = rcp(new idmap_t(env_, comm_, gids_, false));
  }
  Z2_FORWARD_EXCEPTIONS;

  gnosAreGids_ = idMap->gnosAreGids();

  this->setIdentifierMap(idMap);

  gno_t lsum = nLocalIds;
  reduceAll<int, gno_t>(*comm_, Teuchos::REDUCE_SUM, 1, &lsum,
    &numGlobalIdentifiers_);

  if (!gnosAreGids_ && nLocalIds>0){
    gno_t *tmpGno = new gno_t [nLocalIds];
    env_->localMemoryAssertion(__FILE__, __LINE__, nLocalIds, tmpGno);
    gnos_ = arcp(tmpGno, 0, nLocalIds);

    try{
      ArrayRCP<gid_t> gidsNonConst = arcp_const_cast<gid_t>(gids_);
      idMap->gidTranslate( gidsNonConst(0,nLocalIds),  gnos_(0,nLocalIds),
        TRANSLATE_APP_TO_LIB);
    }
    Z2_FORWARD_EXCEPTIONS;
  }

  gnosConst_ = arcp_const_cast<const gno_t>(gnos_);

  env_->memory("After construction of identifier model");
}

#endif // DOXYGEN_SHOULD_SKIP_THIS

}  // namespace Zoltan2

#endif
