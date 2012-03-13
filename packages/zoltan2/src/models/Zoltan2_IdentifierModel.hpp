// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_IdentifierModel.hpp
    \brief Defines the IdentifierModel interface.
*/

#ifndef _ZOLTAN2_IDENTIFIERMODEL_HPP_
#define _ZOLTAN2_IDENTIFIERMODEL_HPP_

#include <Zoltan2_Model.hpp>
#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_IdentifierInput.hpp>
#include <Zoltan2_StridedData.hpp>

namespace Zoltan2 {

/*!  \brief IdentifierModel defines the interface for all identifier models.

    The constructor of the IdentifierModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    methods should be local methods.

    The template parameter is an InputAdapter, which is an object that
    provides a uniform interface for models to the user's input data.

    Explicit instantiations exist for:
      \li MatrixInput
      \li IdentifierInput

    \todo Add instantiations for CoordinateInput, GraphInput, VectorInput
               and MeshInput
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

  int getNumWeights() const { return 0; }
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS

template <typename User>
class IdentifierModel<IdentifierInput<User> > : public Model<IdentifierInput<User> >
{
public:

  typedef typename IdentifierInput<User>::scalar_t  scalar_t;
  typedef typename IdentifierInput<User>::gno_t     gno_t;
  typedef typename IdentifierInput<User>::lno_t     lno_t;
  typedef typename IdentifierInput<User>::gid_t     gid_t;
  typedef IdentifierMap<User> idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  /*! \brief Constructor

       \param ia  the input adapter from which to build the model
       \param env   the application environment (including problem parameters)
       \param comm  the problem communicator
       \param modelFlags   bit map of Zoltan2::IdentifierModelFlags
   */
  
  IdentifierModel( const IdentifierInput<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &modelFlags):
      gnosAreGids_(false), numGlobalIdentifiers_(), env_(env), comm_(comm),
      gids_(), weights_(), gnos_(), gnosConst_()
  {
    int weightDim = ia->getNumberOfWeights();
    size_t nLocalIds = ia->getLocalNumberOfIdentifiers();

    const scalar_t **wgts=NULL;
    int *wgtStrides = NULL;

    if (nLocalIds && weightDim){
      wgts = new const scalar_t * [weightDim];
      wgtStrides = new lno_t [weightDim];
      env_->localMemoryAssertion(__FILE__, __LINE__, nLocalIds, 
        wgts && wgtStrides);
    }

    const gid_t *gids=NULL;

    try{
      ia->getIdentifierList(gids);
      for (int dim=0; dim < weightDim; dim++)
        ia->getIdentifierWeights(dim, wgts[dim], wgtStrides[dim]);
    }
    Z2_FORWARD_EXCEPTIONS;

    if (nLocalIds){
      gids_ = arcp(gids, 0, nLocalIds, false);
  
      if (weightDim > 0){
        input_t *w = new input_t [weightDim];
        for (int i=0; i < weightDim; i++){
          ArrayView<const scalar_t> wgtArray(wgts[i], nLocalIds*wgtStrides[i]);
          w[i] = input_t(wgtArray, wgtStrides[i]);
        }
        weights_ = arcp<input_t>(w, 0, weightDim);
      }
    }

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
  }

  /*! Returns the number identifiers on this process.
   */
  size_t getLocalNumIdentifiers() const { return gids_.size(); }

  /*! Returns the global number identifiers.
   */
  global_size_t getGlobalNumIdentifiers() const {return numGlobalIdentifiers_;}

  /*! Returns the dimension (0 or greater) of identifier weights.
   */
  int getIdentifierWeightDim() const { return weights_.size(); }

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
    size_t n = getLocalNumIdentifiers();
    size_t nweights = n * weights_.size();

    Ids =  ArrayView<const gno_t>();
    wgts = ArrayView<input_t>();

    if (n){
      if (gnosAreGids_)
        Ids = gids_(0, n);
      else
        Ids = gnosConst_(0, n);

      if (nweights)
        wgts = weights_.view(0, weights_.size());
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

  int getNumWeights() const { return getIdentifierWeightDim(); }

private:

  bool gnosAreGids_;
  gno_t numGlobalIdentifiers_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  ArrayRCP<const gid_t> gids_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;
};

template <typename User>
class IdentifierModel<MatrixInput<User> > : public Model<MatrixInput<User> >
{
public:

  typedef typename MatrixInput<User>::scalar_t  scalar_t;
  typedef typename MatrixInput<User>::gno_t     gno_t;
  typedef typename MatrixInput<User>::lno_t     lno_t;
  typedef typename MatrixInput<User>::gid_t     gid_t;
  typedef IdentifierMap<User> idmap_t;
  typedef StridedData<lno_t, scalar_t> input_t;
  
  IdentifierModel( const MatrixInput<User> *ia, 
    const RCP<const Environment> &env, const RCP<const Comm<int> > &comm, 
    modelFlag_t &modelFlags):
      gnosAreGids_(false), numGlobalIdentifiers_(), env_(env), 
      comm_(comm), gids_(), weights_(), gnos_(), gnosConst_()
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

    gnosAreGids_ = idMap->gnosAreGids();

    this->setIdentifierMap(idMap);   // Base Model method

    gno_t lsum = nLocalIds;
    reduceAll<int, gno_t>(*comm_, Teuchos::REDUCE_SUM, 1, &lsum,
      &numGlobalIdentifiers_);

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
  }

  /*! Returns the number identifierson this process.
   */
  size_t getLocalNumIdentifiers() const { return gids_.size(); }

  /*! Returns the global number identifiers.
   */
  global_size_t getGlobalNumIdentifiers() const {return numGlobalIdentifiers_;}

  /*! Returns the dimension (0 or greater) of identifier weights.
   *    Weights are not yet implemented in MatrixInput.
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
    size_t n = getLocalNumIdentifiers();
    size_t nweights = 0;

    Ids = ArrayView<const gno_t>(Teuchos::null);
    wgts = ArrayView<input_t>(Teuchos::null);

    if (n){
      if (gnosAreGids_)
        Ids = gids_(0, n);
      else
        Ids = gnosConst_(0, n);

      if (nweights)
        wgts = weights_(0, nweights);
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

  int getNumWeights() const { return getIdentifierWeightDim(); }

private:

  bool gnosAreGids_;
  gno_t numGlobalIdentifiers_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  ArrayRCP<const gid_t> gids_;
  ArrayRCP<input_t> weights_;
  ArrayRCP<gno_t> gnos_;
  ArrayRCP<const gno_t> gnosConst_;
};

#endif // DOXYGEN_SHOULD_SKIP_THIS

}  // namespace Zoltan2

#endif
