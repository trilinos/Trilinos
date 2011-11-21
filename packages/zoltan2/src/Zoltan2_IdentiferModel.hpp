// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_IdentifierModel.hpp

    \brief The interface and implementations of a simple identifier model.
*/


#ifndef _ZOLTAN2_IDENTIFIERMODEL_HPP_
#define _ZOLTAN2_IDENTIFIERMODEL_HPP_

#include <Zoltan2_Model.hpp>
#include <Zoltan2_IdentifierMap.hpp>

namespace Zoltan2 {

/*! Zoltan2::IdentifierModel
    \brief This provides simple IDs and weights to the Zoltan2 algorithm.

    The template parameter is an Input Adapter.  Input adapters are
    templated on the basic user input type.
*/
template <typename Adapter>
class IdentifierModel : public Model<Adapter> {
private:

public:

  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  typedef typename Adapter::gid_t     gid_t;
  typedef typename Adapter::lid_t     lid_t;
  typedef typename Adapter::node_t    node_t;
  typedef typename Adapter::user_t    user_t;
  
  IdentifierModel(Environment &env, Adapter &ia){
    throw logic_error("a specific instantiation should be used");
  }

  /*! Returns the number identifierson this process.
   */
  size_t getLocalNumIdentifiers() const { return 0; }

  /*! Returns the global number identifiers.
   */
  global_size_t getGlobalNumIdentifiers() const { return 0; }

  /*! Returns the dimension (0 or greater) of identifier weights.
   */
  int getIdentifierWeightDim() const { return 0; }

  /*! Returns the base ID, typically 0 or 1.
   */
  gno_t getIndexBase() const {return 0;}

  /*! Sets pointers to this process' identifier Ids and their weights.
      \param Ids will on return point to the list of the global Ids for
        each identifier on this process.
      \param wgts will on return point to a list of the weight or weights
         associated with each identifier in the Ids list.  Weights are listed by
         identifier by weight component.
       \return The number of ids in the Ids list.
   */

  size_t getIdentifierList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &wgts) const { return 0; }

  /*! Returns the IdentiferMap - the Problem may need it to
   translate Zoltan2 gnos to application gids and lids.
   */
  int getIdentifierMap() const { return ; }
};

template <typename User>
class IdentifierModel<IdentifierInput<User> >
{
private:

  bool gnosAreGids;
  int weightDim_;
  gno_t numGlobalIdentifiers_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  ArrayRCP<const gid_t> gids_;
  ArrayRCP<const lid_t> lids_;
  ArrayRCP<const scalar_t> weights_;
  RCP<IdentifierMap> idMap_;
  ArrayRCP<const gno_t> gnos_;

public:

  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  typedef typename Adapter::gid_t     gid_t;
  typedef typename Adapter::lid_t     lid_t;
  typedef typename Adapter::node_t    node_t;
  typedef typename Adapter::user_t    user_t;
  
  IdentifierModel(const RCP<const Environment> &env, 
    IdentifierModel<User> &ia, bool gnosMustBeConsecutive=false): 
      gnosAreGids(false), weightDim_(0), numGlobalIdentifiers_(), env_(env), 
      comm_(env.comm_), gids_(), lids_(), weights_(), idMap_(), gnos_()
  {
    size_t nLocalIds;
    const gid_t *gids;
    const lid_t *lids;
    const scalar_t *wgts;

    weightDim_ = ia.getIdentifierWeightDim();

    try{
      nLocalIds = ia.getIdentifierView(gids, lids, wgts);
    }
    catch (const std::exception &e){
      Z2_THROW_ZOLTAN2_ERROR(env_, e);
    }

    if (nlocalIds){
      gids_ = arcp(const_cast<gid_t *>(gids), 0, nlocalIds, false);
  
      if (lids){
        lids_ = arcp(const_cast<lid_t *>(lids), 0, nlocalIds, false);
      }
  
      if (wgts){
        weights_ = arcp(const_cast<scalar_t *>(wgts), 0, nlocalIds*weightDim_, 
          false);
      }
    }

    try{
      idMap_ = rcp(new IdentifierMap(comm_, env_, gids_, lids_,
                 gnosMustBeConsecutive));
    }
    catch (const std::exception &e){
      Z2_THROW_ZOLTAN2_ERROR(env_, e);
    }

    gnosAreGids = idMap_->gnosAreGids();

    gno_t lsum = nLocalIds;
    Teuchos::reduceAll<int, gno_t>(*comm_, Teuchos::REDUCE_SUM, 1, &lsum
      &numGlobalIdentifiers_);
  }

  /*! Returns the number identifierson this process.
   */
  size_t getLocalNumIdentifiers() const { return gids_.size(); }

  /*! Returns the global number identifiers.
   */
  global_size_t getGlobalNumIdentifiers() const {return numGlobalIdentifiers_;}

  /*! Returns the dimension (0 or greater) of identifier weights.
   */
  int getIdentifierWeightDim() const { return weightDim_; }

  /*! Returns the base ID, typically 0 or 1.
   *    TODO: base has to be added to the IdentifierMap.
   */
  gno_t getIndexBase() const {return 0;}

  /*! Sets pointers to this process' identifier Ids and their weights.
      \param Ids will on return point to the list of the global Ids for
        each identifier on this process.
      \param wgts will on return point to a list of the weight or weights
         associated with each identifier in the Ids list.  Weights are listed by
         identifier by weight component.
       \return The number of ids in the Ids list.
   */

  size_t getIdentifierList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &wgts) const 
  {
    wgts = weights_
    if (gnosAreGids_){
      Ids = gids_;
    }
    else{
      if (gnos_.size() != gids_.size()){
        // TODO check for failures
        gnos_ = arcp(new gno_t [gids_.size()], 0, gids_.size(), true);
        idMap_->gidTranslate( gids_.view(0, gids.size()), 
          gnos_.view(0, gnos.size()), TRANSLATE_APP_TO_LIB);
      }
      Ids = gnos_;
    }
    return gids.size();
  }

  /*! Returns the IdentiferMap - the Problem may need it to
   *  translate Zoltan2 gnos to application gids and lids.
   */
  RCP<IdentifierMap> getIdentifierMap() const {return idMap_;}

};

template <typename User>
class IdentifierModel<XpetraCrsMatrix<User> >
{
private:

  bool gnosAreGids;
  int weightDim_;
  gno_t numGlobalIdentifiers_;
  const RCP<const Environment> env_;
  const RCP<const Comm<int> > comm_;
  ArrayRCP<const gid_t> gids_;
  ArrayRCP<const lid_t> lids_;
  ArrayRCP<const scalar_t> weights_;
  RCP<IdentifierMap> idMap_;
  ArrayRCP<const gno_t> gnos_;

public:

  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  typedef typename Adapter::gid_t     gid_t;
  typedef typename Adapter::lid_t     lid_t;
  typedef typename Adapter::node_t    node_t;
  typedef typename Adapter::user_t    user_t;
  
  IdentifierModel(const RCP<const Environment> &env, 
    XpetraCrsMatrix<User> &ia, bool gnosMustBeConsecutive=false): 
      gnosAreGids(false), weightDim_(0), numGlobalIdentifiers_(), env_(env), 
      comm_(env.comm_), gids_(), lids_(), weights_(), idMap_(), gnos_()
  {
    size_t nLocalIds;
    const gid_t *gids;
    const lid_t *lids;
    const gid_t *colIds;
    const lno_t *offsets;

    weightDim_ = 0;    // TODO not implemented yet

    try{
      nLocalIds = ia.getRowListView(gids, lids, offsets, colIds);
    }
    catch (const std::exception &e){
      Z2_THROW_ZOLTAN2_ERROR(env_, e);
    }

    if (nlocalIds){
      gids_ = arcp(const_cast<gid_t *>(gids), 0, nlocalIds, false);
  
      if (lids){
        lids_ = arcp(const_cast<lid_t *>(lids), 0, nlocalIds, false);
      }
  
      if (wgts){
        weights_ = arcp(const_cast<scalar_t *>(wgts), 0, nlocalIds*weightDim_, 
          false);
      }
    }

    try{
      idMap_ = rcp(new IdentifierMap(comm_, env_, gids_, lids_,
                 gnosMustBeConsecutive));
    }
    catch (const std::exception &e){
      Z2_THROW_ZOLTAN2_ERROR(env_, e);
    }

    gnosAreGids = idMap_->gnosAreGids();

    gno_t lsum = nLocalIds;
    Teuchos::reduceAll<int, gno_t>(*comm_, Teuchos::REDUCE_SUM, 1, &lsum
      &numGlobalIdentifiers_);
  }

  /*! Returns the number identifierson this process.
   */
  size_t getLocalNumIdentifiers() const { return gids_.size(); }

  /*! Returns the global number identifiers.
   */
  global_size_t getGlobalNumIdentifiers() const {return numGlobalIdentifiers_;}

  /*! Returns the dimension (0 or greater) of identifier weights.
   */
  int getIdentifierWeightDim() const { return weightDim_; }

  /*! Returns the base ID, typically 0 or 1.
   *    TODO: base has to be added to the IdentifierMap.
   */
  gno_t getIndexBase() const {return 0;}

  /*! Sets pointers to this process' identifier Ids and their weights.
      \param Ids will on return point to the list of the global Ids for
        each identifier on this process.
      \param wgts will on return point to a list of the weight or weights
         associated with each identifier in the Ids list.  Weights are listed by
         identifier by weight component.
       \return The number of ids in the Ids list.
   */

  size_t getIdentifierList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &wgts) const 
  {
    wgts = weights_
    if (gnosAreGids_){
      Ids = gids_;
    }
    else{
      if (gnos_.size() != gids_.size()){
        // TODO check for failures
        gnos_ = arcp(new gno_t [gids_.size()], 0, gids_.size(), true);
        idMap_->gidTranslate( gids_.view(0, gids.size()), 
          gnos_.view(0, gnos.size()), TRANSLATE_APP_TO_LIB);
      }
      Ids = gnos_;
    }
    return gids.size();
  }

  /*! Returns the IdentiferMap - the Problem may need it to
   *  translate Zoltan2 gnos to application gids and lids.
   */
  RCP<IdentifierMap> getIdentifierMap() const {return idMap_;}
};

}  // namespace Zoltan2

#endif
