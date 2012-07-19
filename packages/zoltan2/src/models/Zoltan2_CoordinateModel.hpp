// @HEADER
// ***********************************************************************
//                Copyright message goes here.   
// ***********************************************************************
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

  int getCoordinateWeightDim() const { return this->getNumWeights(); }

  size_t getCoordinates(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const 
  {
    size_t n = getLocalNumCoordinates();

    Ids =  ArrayView<const gno_t>();
    xyz = ArrayView<input_t>();
    wgts = ArrayView<input_t>();

    if (n){
      if (gnosAreGids_)
        Ids = Teuchos::arrayView<const gno_t>(
          reinterpret_cast<const gno_t *>(gids_.getRawPtr()), n);
      else
        Ids = gnosConst_.view(0, n);

      xyz =  xyz_.view(0, coordinateDim_);

      int wdim = this->getNumWeights();
      if (wdim > 0)
	wgts = weights_.view(0, wdim);
    }
    
    return n;
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
      coordinateDim_(), gids_(), xyz_(), weights_(), 
      gnos_(), gnosConst_()
{
  bool consecutiveIds = flags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

  size_t nLocalIds = ia->getLocalNumberOfCoordinates();

  // Get coordinates and weights (if any)

  coordinateDim_ = ia->getCoordinateDimension();
  int wdim = ia->getNumberOfWeights();

  Model<CoordinateInput<User> >::maxCount(*comm, coordinateDim_, wdim);

  env_->localBugAssertion(__FILE__, __LINE__, "coordinate dimension",
    coordinateDim_ > 0, COMPLEX_ASSERTION);

  input_t *coordArray = new input_t [coordinateDim_];
  input_t *weightArray = NULL;
  if (wdim)
    weightArray = new input_t [wdim];

  env_->localMemoryAssertion(__FILE__, __LINE__, wdim+coordinateDim_,
    coordArray && (!wdim || weightArray));

  Array<lno_t> arrayLengths(wdim, 0);

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

    for (int dim=0; dim < wdim; dim++){
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

  if (wdim)
    weights_ = arcp(weightArray, 0, wdim);

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
    size_t n = getLocalNumCoordinates();

    Ids =  ArrayView<const gno_t>();
    xyz = ArrayView<input_t>();
    wgts = ArrayView<input_t>();

    if (n){
      if (gnosAreGids_)
        Ids = Teuchos::arrayView<const gno_t>(
          reinterpret_cast<const gno_t *>(gids_.getRawPtr()), n);
      else
        Ids = gnosConst_.view(0, n);

      xyz =  xyz_.view(0, coordinateDim_);

      int wdim = this->getNumWeights();
      if (wdim)
	wgts = weights_.view(0, wdim);
    }
    
    return n;
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
      coordinateDim_(), gids_(), xyz_(), weights_(), 
      gnos_(), gnosConst_()
{
  bool consecutiveIds = flags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

  int wdim = 0;  // at the present time, matrix input does not have weights

  coordinateDim_ = ia->getCoordinateDimension();

  Model<MatrixInput<User> >::maxCount(*comm, coordinateDim_, wdim);

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

  int getCoordinateWeightDim() const { return this->getNumWeights();}

  size_t getCoordinates(ArrayView<const gno_t>  &Ids,
    ArrayView<input_t> &xyz,
    ArrayView<input_t> &wgts) const
  {
    size_t n = getLocalNumCoordinates();

    Ids =  ArrayView<const gno_t>();
    xyz = ArrayView<input_t>();
    wgts = ArrayView<input_t>();

    if (n){
      if (gnosAreGids_)
        Ids = Teuchos::arrayView<const gno_t>(
          reinterpret_cast<const gno_t *>(gids_.getRawPtr()), n);
      else
        Ids = gnosConst_.view(0, n);

      xyz =  xyz_.view(0, coordinateDim_);

      int wdim = this->getNumWeights();
      if (wdim > 0)
        wgts = weights_.view(0, wdim);
    }

    return n;
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
      coordinateDim_(), gids_(), xyz_(), weights_(), 
      gnos_(), gnosConst_()
{
  bool consecutiveIds = flags.test(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

  size_t nLocalIds = ia->getLocalLength();

  // Get coordinates and weights (if any)

  coordinateDim_ = ia->getNumberOfVectors();
  int wdim = ia->getNumberOfWeights();

  Model<VectorInput<User> >::maxCount(*comm, coordinateDim_, wdim);

  env_->localBugAssertion(__FILE__, __LINE__, "coordinate dimension",
    coordinateDim_ > 0, COMPLEX_ASSERTION);

  input_t *coordArray = new input_t [coordinateDim_];
  input_t *weightArray = NULL;
  if (wdim)
    weightArray = new input_t [wdim];

  env_->localMemoryAssertion(__FILE__, __LINE__, wdim+coordinateDim_,
    coordArray && (!wdim|| weightArray));

  Array<lno_t> arrayLengths(wdim, 0);

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

    for (int dim=0; dim < wdim; dim++){
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

  if (wdim)
    weights_ = arcp(weightArray, 0, wdim);

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


#endif   // DOXYGEN_SHOULD_SKIP_THIS

}   // namespace Zoltan2

#endif
