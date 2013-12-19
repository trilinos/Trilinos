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

// disable clang warnings
#ifdef __clang__
#pragma clang system_header
#endif

#include <Zoltan2_Model.hpp>
#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_IdentifierAdapter.hpp>
#include <Zoltan2_VectorAdapter.hpp>
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
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS

////////////////////////////////////////////////////////////////
// Coordinate model derived from MatrixAdapter.
////////////////////////////////////////////////////////////////

template <typename User>
class CoordinateModel<MatrixAdapter<User> > : 
  public Model<MatrixAdapter<User> >
{
public:

  typedef typename MatrixAdapter<User>::scalar_t  scalar_t;
  typedef typename MatrixAdapter<User>::gno_t     gno_t;
  typedef typename MatrixAdapter<User>::lno_t     lno_t;
  typedef typename MatrixAdapter<User>::gid_t     gid_t;
  typedef IdentifierMap<User> idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  /*! \brief Constructor

       \param ia  the input adapter from which to build the model
       \param env   the application environment (including problem parameters)
       \param comm  the problem communicator
       \param  modelFlags    a bit map of ModelFlags
   */
  
  CoordinateModel( const MatrixAdapter<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags);

  int getCoordinateDim() const { return coordinateDim_; }

  size_t getLocalNumCoordinates() const { return gids_.size(); }

  global_size_t getGlobalNumCoordinates() const {return numGlobalCoordinates_;}

  /* \brief Weights are not implemented in MatrixAdapter.
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
  CoordinateModel<MatrixAdapter<User> >::CoordinateModel( 
    const MatrixAdapter<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags):
      gnosAreGids_(false), numGlobalCoordinates_(), env_(env), comm_(comm),
      coordinateDim_(), gids_(), xyz_(), userNumWeights_(0), weights_(), 
      gnos_(), gnosConst_()
{
  this->env_->debug(DETAILED_STATUS, "CoordinateModel<MatrixAdapter>");

  bool consecutiveIds = flags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

  userNumWeights_= 0;  // matrix input does not have weights

  coordinateDim_ = ia->getCoordinateDimension();

  Model<MatrixAdapter<User> >::maxCount(*comm, coordinateDim_, userNumWeights_);

  env_->localBugAssertion(__FILE__, __LINE__, "coordinate dimension",
    coordinateDim_ > 0, COMPLEX_ASSERTION);

  size_t nLocalIds = (coordinateDim_ ? ia->getLocalNumIDs() : 0);

  // Get coordinates

  this->env_->debug(DETAILED_STATUS, "    getting coordinates");
  input_t *coordArray = new input_t [coordinateDim_];
  env_->localMemoryAssertion(__FILE__, __LINE__, coordinateDim_, coordArray);

  if (nLocalIds){
    for (int dim=0; dim < coordinateDim_; dim++){
      int stride;
      const scalar_t *coords=NULL;
      try{
        ia->getCoordinatesView(coords, stride, dim);
      }
      Z2_FORWARD_EXCEPTIONS;

      ArrayRCP<const scalar_t> cArray(coords, 0, nLocalIds*stride, false);
      coordArray[dim] = input_t(cArray, stride);
    }
  }

  xyz_ = arcp(coordArray, 0, coordinateDim_);

  // Create identifier map.

  this->env_->debug(DETAILED_STATUS, "    getting identifiers");
  const gid_t *gids;
  ia->getIDsView(gids);
  gids_ = arcp<const gid_t>(gids, 0, nLocalIds, false);

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
  this->env_->debug(DETAILED_STATUS, "    finished");
}

////////////////////////////////////////////////////////////////
// Coordinate model derived from VectorAdapter.
////////////////////////////////////////////////////////////////

template <typename User>
class CoordinateModel<VectorAdapter<User> > : 
  public Model<VectorAdapter<User> >
{

public:

  typedef typename VectorAdapter<User>::scalar_t  scalar_t;
  typedef typename VectorAdapter<User>::gno_t     gno_t;
  typedef typename VectorAdapter<User>::gid_t     gid_t;
  typedef typename VectorAdapter<User>::lno_t     lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  typedef IdentifierMap<User> idmap_t;

  CoordinateModel( const VectorAdapter<User> *ia, 
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
CoordinateModel<VectorAdapter<User> >::CoordinateModel( 
    const VectorAdapter<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t flags):
      gnosAreGids_(false), numGlobalCoordinates_(), env_(env), comm_(comm),
      coordinateDim_(), gids_(), xyz_(), userNumWeights_(0), weights_(), 
      gnos_(), gnosConst_()
{
  bool consecutiveIds = flags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

  size_t nLocalIds = ia->getLocalNumIDs();

  // Get coordinates and weights (if any)

  coordinateDim_ = ia->getNumEntriesPerID();
  userNumWeights_ = ia->getNumWeightsPerID();

  Model<VectorAdapter<User> >::maxCount(*comm, coordinateDim_, userNumWeights_);

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
    const gid_t *gids=NULL;
    ia->getIDsView(gids);
    gids_ = arcp(gids, 0, nLocalIds, false);

    for (int dim=0; dim < coordinateDim_; dim++){
      int stride;
      const scalar_t *coords=NULL;
      try{
        ia->getEntriesView(coords, stride, dim);
      }
      Z2_FORWARD_EXCEPTIONS;

      ArrayRCP<const scalar_t> cArray(coords, 0, nLocalIds*stride, false);
      coordArray[dim] = input_t(cArray, stride);
    }

    for (int idx=0; idx < userNumWeights_; idx++){
      int stride;
      const scalar_t *weights;
      try{
        ia->getWeightsView(weights, stride, idx);
      }
      Z2_FORWARD_EXCEPTIONS;

      if (weights){
        ArrayRCP<const scalar_t> wArray(weights, 0, nLocalIds*stride, false);
        weightArray[idx] = input_t(wArray, stride);
        arrayLengths[idx] = nLocalIds;
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
// Coordinate model derived from IdentifierAdapter.
// A coordinate model can not be built from IdentifierAdapter.
// This specialization exists so that other code can compile.
////////////////////////////////////////////////////////////////

template <typename User>
class CoordinateModel<IdentifierAdapter<User> > : 
  public Model<IdentifierAdapter<User> >
{
public:

  typedef typename IdentifierAdapter<User>::scalar_t  scalar_t;
  typedef typename IdentifierAdapter<User>::gno_t     gno_t;
  typedef typename IdentifierAdapter<User>::lno_t     lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  
  CoordinateModel( const IdentifierAdapter<User> *ia, 
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
};


////////////////////////////////////////////////////////////////
// Coordinate model derived from GraphAdapter.
////////////////////////////////////////////////////////////////

template <typename User>
class CoordinateModel<GraphAdapter<User> > : 
  public Model<GraphAdapter<User> >
{

public:

  typedef typename GraphAdapter<User>::scalar_t  scalar_t;
  typedef typename GraphAdapter<User>::gno_t     gno_t;
  typedef typename GraphAdapter<User>::gid_t     gid_t;
  typedef typename GraphAdapter<User>::lno_t     lno_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  typedef IdentifierMap<User> idmap_t;

  CoordinateModel( const GraphAdapter<User> *ia, 
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
CoordinateModel<GraphAdapter<User> >::CoordinateModel( 
    const GraphAdapter<User> *ia, 
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

  // CoordinateModel is built with points == GRAPH_VERTEX from GraphAdapter.
  // It is not ready to use points == GRAPH_EDGE from GraphAdapter.
  env->localInputAssertion(__FILE__, __LINE__, 
   "CoordinateModel from GraphAdapter is implemented only for "
   "Graph Vertices as primary object, not for Graph Edges", 
   ia->getPrimaryEntityType() == Zoltan2::GRAPH_VERTEX, BASIC_ASSERTION);

  // Get coordinates and weights (if any)

  userNumWeights_ = ia->getNumWeightsPerVertex();

  Model<GraphAdapter<User> >::maxCount(*comm, coordinateDim_, userNumWeights_);

  input_t *coordArray = new input_t [coordinateDim_];
  input_t *weightArray = NULL;
  if (userNumWeights_)
    weightArray = new input_t [userNumWeights_];

  env_->localMemoryAssertion(__FILE__, __LINE__, userNumWeights_+coordinateDim_,
    coordArray && (!userNumWeights_|| weightArray));

  Array<lno_t> arrayLengths(userNumWeights_, 0);

  size_t nLocalIds = ia->getLocalNumVertices();

  if (nLocalIds){
    const gid_t *globalIds=NULL;
    const lno_t *offsets=NULL;
    const gid_t *edgeIds=NULL;

    ia->getVertexIDsView(globalIds, offsets, edgeIds);
    gids_ = arcp(globalIds, 0, nLocalIds, false);

    for (int dim=0; dim < coordinateDim_; dim++){
      int stride;
      const scalar_t *coords=NULL;
      try{
        ia->getVertexCoordinatesView(coords, stride, dim);
      }
      Z2_FORWARD_EXCEPTIONS;

      ArrayRCP<const scalar_t> cArray(coords, 0, nLocalIds*stride, false);
      coordArray[dim] = input_t(cArray, stride);
    }

    for (int idx=0; idx < userNumWeights_; idx++){
      int stride;
      const scalar_t *weights;
      try{
        ia->getVertexWeightsView(weights, stride, idx);
      }
      Z2_FORWARD_EXCEPTIONS;

      if (weights){
        ArrayRCP<const scalar_t> wArray(weights, 0, nLocalIds*stride, false);
        weightArray[idx] = input_t(wArray, stride);
        arrayLengths[idx] = nLocalIds;
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
