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

/*! \file Zoltan2_CoordinateModel.hpp
    \brief Defines the CoordinateModel classes.
*/


#ifndef _ZOLTAN2_COORDINATEMODEL_HPP_
#define _ZOLTAN2_COORDINATEMODEL_HPP_

#include <Zoltan2_Model.hpp>
#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_GraphInput.hpp>
#include <Zoltan2_IdentifierInput.hpp>
#include <Zoltan2_CoordinateInput.hpp>
#include <Zoltan2_VectorInput.hpp>
#include <Zoltan2_StridedData.hpp>

namespace Zoltan2 {

/*!  \brief This class provides geometric coordinates with optional weights 
           to the Zoltan2 algorithm.

    The template parameter is an Input Adapter.  Input adapters are
    templated on the basic user input type.
*/
template <typename Adapter>
class CoordinateModel : public Model<Adapter> 
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
#endif
  
  CoordinateModel(){
    throw std::logic_error("a specific instantiation should be used");
  }

  ////////////////////////////////////////////////////
  // The CoordinateModel interface.
  ////////////////////////////////////////////////////

  /*! \brief Returns the dimension of the coordinates.
   */
  int getCoordinateDim() const { return 0; }

  /*! \brief Returns the number of coordinates on this process.
   */
  size_t getLocalNumCoordinates() const { return 0; }

  /*! \brief Returns the global number coordinates.
   */
  global_size_t getGlobalNumCoordinates() const { return 0; }

  /*! \brief Returns the dimension (0 or greater) of coordinate weights.
   */
  int getCoordinateWeightDim() const { return 0; }

  /*! \brief Returns the coordinate ids, values and optional weights.

      \param Ids will on return point to the list of the global Ids for
        each coordinate on this process.

      \param xyz on return is a list of getCoordinateDim() 
          StridedData objects, each containing the coordinates for 
          one dimension. If the coordinate dimension is three, then
          the coordinates for <tt>Ids[k]</tt> are 
          <tt>xyz[0][k], xyz[1][k], xyz[2][k]</tt>.

      \param wgts on return is a list of getCoordinateWeightDim() 
          StridedData objects, each containing the weights for 
          one weight dimension.  For the dimension \d,
          the weight for <tt>Ids[k]</tt> is <tt>wgts[d][k]</tt>.

       \return The number of ids in the Ids list.

      Memory for this data is allocated either by the user or the Model.
      The caller gets a view of the data.
   */

  size_t getCoordinates(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const {return 0;}

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumCoordinates();
  }

  size_t getGlobalNumObjects() const
  {
    return getGlobalNumCoordinates();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const { return ; }
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS

////////////////////////////////////////////////////////////////
// Coordinate model derived from CoordinateInput.
////////////////////////////////////////////////////////////////

template <typename User>
class CoordinateModel<CoordinateInput<User> > : 
  public Model<CoordinateInput<User> >
{

public:

  typedef typename CoordinateInput<User>::scalar_t  scalar_t;
  typedef typename CoordinateInput<User>::gno_t     gno_t;
  typedef typename CoordinateInput<User>::lno_t     lno_t;
  typedef typename CoordinateInput<User>::gid_t     gid_t;
  typedef IdentifierMap<User> idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  /*! \brief Constructor

       \param ia  the input adapter from which to build the model
       \param env   the application environment (including problem parameters)
       \param comm  the problem communicator
       \param  modelFlags    a bit map of ModelFlags
   */
  
  CoordinateModel( const CoordinateInput<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags);

  int getCoordinateDim() const { return coordinateDim_; }

  size_t getLocalNumCoordinates() const { return gids_.size(); }

  global_size_t getGlobalNumCoordinates() const {return numGlobalCoordinates_;}

  int getCoordinateWeightDim() const { return userNumWeights_;}

  size_t getCoordinates(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const 
  {
    xyz = xyz_.view(0, coordinateDim_);
    wgts = weights_.view(0, userNumWeights_);

    size_t nCoord = getLocalNumCoordinates();
    Ids =  ArrayView<const gno_t>();

    if (nCoord){
      if (gnosAreGids_)
        Ids = Teuchos::arrayView<const gno_t>(
          reinterpret_cast<const gno_t *>(gids_.getRawPtr()), nCoord);
      else
        Ids = gnosConst_.view(0, nCoord);
    }
    
    return nCoord;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumCoordinates();
  }

  size_t getGlobalNumObjects() const
  {
    return getGlobalNumCoordinates();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const 
  { 
    ArrayView<input_t> xyz;
    ArrayView<input_t> weights;
    getCoordinates(gnos, xyz, weights);
  }

private:

  bool gnosAreGids_;
  gno_t numGlobalCoordinates_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  int coordinateDim_;
  ArrayRCP<const gid_t> gids_;
  ArrayRCP<input_t> xyz_;
  int userNumWeights_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;

};

template <typename User>
  CoordinateModel<CoordinateInput<User> >::CoordinateModel( 
    const CoordinateInput<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags):
      gnosAreGids_(false), numGlobalCoordinates_(), env_(env), comm_(comm),
      coordinateDim_(), gids_(), xyz_(), userNumWeights_(0), weights_(), 
      gnos_(), gnosConst_()
{
  bool consecutiveIds = flags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

  size_t nLocalIds = ia->getLocalNumberOfCoordinates();

  // Get coordinates and weights (if any)

  coordinateDim_ = ia->getCoordinateDimension();
  userNumWeights_ = ia->getNumberOfWeights();

  Model<CoordinateInput<User> >::maxCount(*comm, coordinateDim_, 
    userNumWeights_);

  env_->localBugAssertion(__FILE__, __LINE__, "coordinate dimension",
    coordinateDim_ > 0, COMPLEX_ASSERTION);

  input_t *coordArray = new input_t [coordinateDim_];
  input_t *weightArray = NULL;
  if (userNumWeights_)
    weightArray = new input_t [userNumWeights_];

  env_->localMemoryAssertion(__FILE__, __LINE__, userNumWeights_+coordinateDim_,
    coordArray && (!userNumWeights_ || weightArray));

  Array<lno_t> arrayLengths(userNumWeights_, 0);

  if (nLocalIds){
    for (int dim=0; dim < coordinateDim_; dim++){
      int stride;
      const gid_t *gids=NULL;
      const scalar_t *coords=NULL;
      try{
        ia->getCoordinates(dim, gids, coords, stride);
      }
      Z2_FORWARD_EXCEPTIONS;

      ArrayRCP<const scalar_t> cArray(coords, 0, nLocalIds*stride, false);
      coordArray[dim] = input_t(cArray, stride);

      if (dim==0)
        gids_ = arcp(gids, 0, nLocalIds, false);
    }

    for (int dim=0; dim < userNumWeights_; dim++){
      int stride;
      const scalar_t *weights;
      try{
        ia->getCoordinateWeights(dim, weights, stride);
      }
      Z2_FORWARD_EXCEPTIONS;

      if (weights){
        ArrayRCP<const scalar_t> wArray(weights, 0, nLocalIds*stride, false);
        weightArray[dim] = input_t(wArray, stride);
        arrayLengths[dim] = nLocalIds;
      }
    }
  }

  xyz_ = arcp(coordArray, 0, coordinateDim_);

  if (userNumWeights_)
    weights_ = arcp(weightArray, 0, userNumWeights_);

  // Sets this->weightDim_, this->uniform_
  this->setWeightArrayLengths(arrayLengths, *comm_);

  // Create identifier map.

  RCP<const idmap_t> idMap;

  try{
    idMap = rcp(new idmap_t(env_, comm_, gids_, consecutiveIds));
  }
  Z2_FORWARD_EXCEPTIONS;

  numGlobalCoordinates_ = idMap->getGlobalNumberOfIds();
  gnosAreGids_ = idMap->gnosAreGids();

  this->setIdentifierMap(idMap);

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

  env_->memory("After construction of coordinate model");
}

////////////////////////////////////////////////////////////////
// Coordinate model derived from MatrixInput.
////////////////////////////////////////////////////////////////

template <typename User>
class CoordinateModel<MatrixInput<User> > : 
  public Model<MatrixInput<User> >
{
public:

  typedef typename MatrixInput<User>::scalar_t  scalar_t;
  typedef typename MatrixInput<User>::gno_t     gno_t;
  typedef typename MatrixInput<User>::lno_t     lno_t;
  typedef typename MatrixInput<User>::gid_t     gid_t;
  typedef IdentifierMap<User> idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  /*! \brief Constructor

       \param ia  the input adapter from which to build the model
       \param env   the application environment (including problem parameters)
       \param comm  the problem communicator
       \param  modelFlags    a bit map of ModelFlags
   */
  
  CoordinateModel( const MatrixInput<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags);

  int getCoordinateDim() const { return coordinateDim_; }

  size_t getLocalNumCoordinates() const { return gids_.size(); }

  global_size_t getGlobalNumCoordinates() const {return numGlobalCoordinates_;}

  /* \brief Weights are not implemented in MatrixInput.
   *
   *   Whereas a computational model may create weights to use with
   *   a matrix problem, they are not inherent in the input.
   */
  int getCoordinateWeightDim() const { return 0; }

  size_t getCoordinates(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const 
  {
    xyz = xyz_.view(0, coordinateDim_);
    wgts = weights_.view(0, userNumWeights_);

    size_t nCoord = getLocalNumCoordinates();
    Ids =  ArrayView<const gno_t>();

    if (nCoord){
      if (gnosAreGids_)
        Ids = Teuchos::arrayView<const gno_t>(
          reinterpret_cast<const gno_t *>(gids_.getRawPtr()), nCoord);
      else
        Ids = gnosConst_.view(0, nCoord);
    }

    return nCoord;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumCoordinates();
  }

  size_t getGlobalNumObjects() const
  {
    return getGlobalNumCoordinates();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const 
  { 
    ArrayView<input_t> xyz;
    ArrayView<input_t> weights;
    getCoordinates(gnos, xyz, weights);
  }

private:

  bool gnosAreGids_;
  gno_t numGlobalCoordinates_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  int coordinateDim_;
  ArrayRCP<const gid_t> gids_;
  ArrayRCP<input_t> xyz_;
  int userNumWeights_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;
};

template <typename User>
  CoordinateModel<MatrixInput<User> >::CoordinateModel( 
    const MatrixInput<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags):
      gnosAreGids_(false), numGlobalCoordinates_(), env_(env), comm_(comm),
      coordinateDim_(), gids_(), xyz_(), userNumWeights_(0), weights_(), 
      gnos_(), gnosConst_()
{
  bool consecutiveIds = flags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

  userNumWeights_= 0;  // matrix input does not have weights

  coordinateDim_ = ia->getCoordinateDimension();

  Model<MatrixInput<User> >::maxCount(*comm, coordinateDim_, userNumWeights_);

  env_->localBugAssertion(__FILE__, __LINE__, "coordinate dimension",
    coordinateDim_ > 0, COMPLEX_ASSERTION);

  size_t nLocalIds = (coordinateDim_ ? ia->getLocalNumRows() : 0);

  // Get coordinates

  input_t *coordArray = new input_t [coordinateDim_];
  env_->localMemoryAssertion(__FILE__, __LINE__, coordinateDim_, coordArray);

  if (nLocalIds){
    for (int dim=0; dim < coordinateDim_; dim++){
      int stride;
      const scalar_t *coords=NULL;
      try{
        ia->getRowCoordinates(dim, coords, stride);
      }
      Z2_FORWARD_EXCEPTIONS;

      ArrayRCP<const scalar_t> cArray(coords, 0, nLocalIds*stride, false);
      coordArray[dim] = input_t(cArray, stride);
    }
  }

  xyz_ = arcp(coordArray, 0, coordinateDim_);

  // Create identifier map.

  const gid_t *rowGids;
  const lno_t *offsets;
  const gid_t *colGids;

  ia->getRowListView(rowGids, offsets, colGids);

  gids_ = arcp<const gid_t>(rowGids, 0, nLocalIds, false);

  RCP<const idmap_t> idMap;

  try{
    idMap = rcp(new idmap_t(env_, comm_, gids_, consecutiveIds));
  }
  Z2_FORWARD_EXCEPTIONS;

  numGlobalCoordinates_ = idMap->getGlobalNumberOfIds();
  gnosAreGids_ = idMap->gnosAreGids();

  this->setIdentifierMap(idMap);

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

  env_->memory("After construction of coordinate model");
}

////////////////////////////////////////////////////////////////
// Coordinate model derived from VectorInput.
////////////////////////////////////////////////////////////////

template <typename User>
class CoordinateModel<VectorInput<User> > : 
  public Model<VectorInput<User> >
{

public:

  typedef typename VectorInput<User>::scalar_t  scalar_t;
  typedef typename VectorInput<User>::gno_t     gno_t;
  typedef typename VectorInput<User>::gid_t     gid_t;
  typedef typename VectorInput<User>::lno_t     lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  typedef IdentifierMap<User> idmap_t;

  CoordinateModel( const VectorInput<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags);

  int getCoordinateDim() const { return coordinateDim_;}

  size_t getLocalNumCoordinates() const { return gids_.size();}

  global_size_t getGlobalNumCoordinates() const {return numGlobalCoordinates_;}

  int getCoordinateWeightDim() const { return userNumWeights_;}

  size_t getCoordinates(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const
  {
    xyz = xyz_.view(0, coordinateDim_);
    wgts = weights_.view(0, userNumWeights_);

    size_t nCoord = getLocalNumCoordinates();
    Ids =  ArrayView<const gno_t>();

    if (nCoord){
      if (gnosAreGids_)
        Ids = Teuchos::arrayView<const gno_t>(
          reinterpret_cast<const gno_t *>(gids_.getRawPtr()), nCoord);
      else
        Ids = gnosConst_.view(0, nCoord);
    }

    return nCoord;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const {return gids_.size();}

  size_t getGlobalNumObjects() const {return numGlobalCoordinates_;}

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const {
    ArrayView<input_t> a, b;
    getCoordinates(gnos, a, b);
    return;
  }

private:

  bool gnosAreGids_;
  gno_t numGlobalCoordinates_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  int coordinateDim_;
  ArrayRCP<const gid_t> gids_;
  ArrayRCP<input_t> xyz_;
  int userNumWeights_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;
};

template <typename User>
CoordinateModel<VectorInput<User> >::CoordinateModel( 
    const VectorInput<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags):
      gnosAreGids_(false), numGlobalCoordinates_(), env_(env), comm_(comm),
      coordinateDim_(), gids_(), xyz_(), userNumWeights_(0), weights_(), 
      gnos_(), gnosConst_()
{
  bool consecutiveIds = flags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

  size_t nLocalIds = ia->getLocalLength();

  // Get coordinates and weights (if any)

  coordinateDim_ = ia->getNumberOfVectors();
  userNumWeights_ = ia->getNumberOfWeights();

  Model<VectorInput<User> >::maxCount(*comm, coordinateDim_, userNumWeights_);

  env_->localBugAssertion(__FILE__, __LINE__, "coordinate dimension",
    coordinateDim_ > 0, COMPLEX_ASSERTION);

  input_t *coordArray = new input_t [coordinateDim_];
  input_t *weightArray = NULL;
  if (userNumWeights_)
    weightArray = new input_t [userNumWeights_];

  env_->localMemoryAssertion(__FILE__, __LINE__, userNumWeights_+coordinateDim_,
    coordArray && (!userNumWeights_|| weightArray));

  Array<lno_t> arrayLengths(userNumWeights_, 0);

  if (nLocalIds){
    for (int dim=0; dim < coordinateDim_; dim++){
      int stride;
      const gid_t *gids=NULL;
      const scalar_t *coords=NULL;
      try{
        ia->getVector(dim, gids, coords, stride);
      }
      Z2_FORWARD_EXCEPTIONS;

      ArrayRCP<const scalar_t> cArray(coords, 0, nLocalIds*stride, false);
      coordArray[dim] = input_t(cArray, stride);

      if (dim==0)
        gids_ = arcp(gids, 0, nLocalIds, false);
    }

    for (int dim=0; dim < userNumWeights_; dim++){
      int stride;
      const scalar_t *weights;
      try{
        ia->getVectorWeights(dim, weights, stride);
      }
      Z2_FORWARD_EXCEPTIONS;

      if (weights){
        ArrayRCP<const scalar_t> wArray(weights, 0, nLocalIds*stride, false);
        weightArray[dim] = input_t(wArray, stride);
        arrayLengths[dim] = nLocalIds;
      }
    }
  }

  xyz_ = arcp(coordArray, 0, coordinateDim_);

  if (userNumWeights_)
    weights_ = arcp(weightArray, 0, userNumWeights_);

  // Sets this->weightDim_ and this->uniform_
  this->setWeightArrayLengths(arrayLengths, *comm_);

  // Create identifier map.

  RCP<const idmap_t> idMap;

  try{
    idMap = rcp(new idmap_t(env_, comm_, gids_, consecutiveIds));
  }
  Z2_FORWARD_EXCEPTIONS;

  numGlobalCoordinates_ = idMap->getGlobalNumberOfIds();
  gnosAreGids_ = idMap->gnosAreGids();

  this->setIdentifierMap(idMap);

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

  env_->memory("After construction of coordinate model");
}

////////////////////////////////////////////////////////////////
// Coordinate model derived from IdentifierInput.
// A coordinate model can not be built from IdentifierInput.
// This specialization exists so that other code can compile.
////////////////////////////////////////////////////////////////

template <typename User>
class CoordinateModel<IdentifierInput<User> > : 
  public Model<IdentifierInput<User> >
{
public:

  typedef typename IdentifierInput<User>::scalar_t  scalar_t;
  typedef typename IdentifierInput<User>::gno_t     gno_t;
  typedef typename IdentifierInput<User>::lno_t     lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  
  CoordinateModel( const IdentifierInput<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags)
  {
    throw logic_error(
      "a coordinate model can not be build from identifier input");
  }

  int getCoordinateDim() const { return 0;}
  size_t getLocalNumCoordinates() const { return 0;}
  global_size_t getGlobalNumCoordinates() const {return 0;}
  int getCoordinateWeightDim() const { return 0;}
  size_t getCoordinates(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const { return 0;}

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const {return 0;}
  size_t getGlobalNumObjects() const {return 0;}
  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const {return;}
};


////////////////////////////////////////////////////////////////
// Coordinate model derived from GraphInput.
////////////////////////////////////////////////////////////////

template <typename User>
class CoordinateModel<GraphInput<User> > : 
  public Model<GraphInput<User> >
{

public:

  typedef typename GraphInput<User>::scalar_t  scalar_t;
  typedef typename GraphInput<User>::gno_t     gno_t;
  typedef typename GraphInput<User>::gid_t     gid_t;
  typedef typename GraphInput<User>::lno_t     lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  typedef IdentifierMap<User> idmap_t;

  CoordinateModel( const GraphInput<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags);

  int getCoordinateDim() const { return coordinateDim_;}

  size_t getLocalNumCoordinates() const { return gids_.size();}

  global_size_t getGlobalNumCoordinates() const {return numGlobalCoordinates_;}

  int getCoordinateWeightDim() const { return userNumWeights_;}

  size_t getCoordinates(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const
  {
    xyz = xyz_.view(0, coordinateDim_);
    wgts = weights_.view(0, userNumWeights_);

    size_t nCoord = getLocalNumCoordinates();
    Ids =  ArrayView<const gno_t>();

    if (nCoord){
      if (gnosAreGids_)
        Ids = Teuchos::arrayView<const gno_t>(
          reinterpret_cast<const gno_t *>(gids_.getRawPtr()), nCoord);
      else
        Ids = gnosConst_.view(0, nCoord);
    }

    return nCoord;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const {return gids_.size();}

  size_t getGlobalNumObjects() const {return numGlobalCoordinates_;}

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const {
    ArrayView<input_t> a, b;
    getCoordinates(gnos, a, b);
    return;
  }

private:

  bool gnosAreGids_;
  gno_t numGlobalCoordinates_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  int coordinateDim_;
  ArrayRCP<const gid_t> gids_;
  ArrayRCP<input_t> xyz_;
  int userNumWeights_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;
};

template <typename User>
CoordinateModel<GraphInput<User> >::CoordinateModel( 
    const GraphInput<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags):
      gnosAreGids_(false), numGlobalCoordinates_(), env_(env), comm_(comm),
      coordinateDim_(), gids_(), xyz_(), userNumWeights_(0), weights_(), 
      gnos_(), gnosConst_()
{
  coordinateDim_ = ia->getCoordinateDimension();

  env->localInputAssertion(__FILE__, __LINE__, 
   "graph input does not have vertex coordinates",
   coordinateDim_>0, BASIC_ASSERTION);

  // Get coordinates and weights (if any)

  userNumWeights_ = ia->getVertexWeightDimension();

  Model<GraphInput<User> >::maxCount(*comm, coordinateDim_, userNumWeights_);

  input_t *coordArray = new input_t [coordinateDim_];
  input_t *weightArray = NULL;
  if (userNumWeights_)
    weightArray = new input_t [userNumWeights_];

  env_->localMemoryAssertion(__FILE__, __LINE__, userNumWeights_+coordinateDim_,
    coordArray && (!userNumWeights_|| weightArray));

  Array<lno_t> arrayLengths(userNumWeights_, 0);

  size_t nLocalIds = ia->getLocalNumberOfVertices();

  if (nLocalIds){
    const gid_t *globalIds=NULL;
    const lno_t *offsets=NULL;
    const gid_t *edgeIds=NULL;

    size_t numIds = ia->getVertexListView(globalIds, offsets, edgeIds);

    gids_ = arcp(globalIds, 0, nLocalIds, false);

    for (int dim=0; dim < coordinateDim_; dim++){
      int stride;
      const scalar_t *coords=NULL;
      try{
        ia->getVertexCoordinates(dim, coords, stride);
      }
      Z2_FORWARD_EXCEPTIONS;

      ArrayRCP<const scalar_t> cArray(coords, 0, nLocalIds*stride, false);
      coordArray[dim] = input_t(cArray, stride);
    }

    for (int dim=0; dim < userNumWeights_; dim++){
      int stride;
      const scalar_t *weights;
      try{
        ia->getVertexWeights(dim, weights, stride);
      }
      Z2_FORWARD_EXCEPTIONS;

      if (weights){
        ArrayRCP<const scalar_t> wArray(weights, 0, nLocalIds*stride, false);
        weightArray[dim] = input_t(wArray, stride);
        arrayLengths[dim] = nLocalIds;
      }
    }
  }

  xyz_ = arcp(coordArray, 0, coordinateDim_);

  if (userNumWeights_)
    weights_ = arcp(weightArray, 0, userNumWeights_);

  // Sets this->weightDim_ and this->uniform_
  this->setWeightArrayLengths(arrayLengths, *comm_);

  // Create identifier map.

  bool consecutiveIds = flags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

  RCP<const idmap_t> idMap;

  try{
    idMap = rcp(new idmap_t(env_, comm_, gids_, consecutiveIds));
  }
  Z2_FORWARD_EXCEPTIONS;

  numGlobalCoordinates_ = idMap->getGlobalNumberOfIds();
  gnosAreGids_ = idMap->gnosAreGids();

  this->setIdentifierMap(idMap);

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

  env_->memory("After construction of coordinate model");
}
#endif   // DOXYGEN_SHOULD_SKIP_THIS

}   // namespace Zoltan2

#endif
