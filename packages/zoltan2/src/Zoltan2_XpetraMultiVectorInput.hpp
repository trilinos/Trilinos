// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraMultiVectorInput.hpp

    \brief An input adapter for Xpetra::MultiVector.
*/

#ifndef _ZOLTAN2_XPETRAMULTIVECTORINPUT_HPP_
#define _ZOLTAN2_XPETRAMULTIVECTORINPUT_HPP_

#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_MultiVectorInput.hpp>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////////
/*! Zoltan2::XpetraMultiVectorInput
    \brief Provides access for Zoltan2 to Xpetra::MultiVector.

    The template parameter is the user's input object: 
     Epetra_MultiVector
     Tpetra::MultiVector
     Xpetra::MultiVector
*/

template <typename User>
class XpetraMultiVectorInput : public MultiVectorInput<User> {
public:

  typedef typename InputAdapter<User>::scalar_t scalar_t;
  typedef typename InputAdapter<User>::lno_t    lno_t;
  typedef typename InputAdapter<User>::gno_t    gno_t;
  typedef typename InputAdapter<User>::lid_t    lid_t;
  typedef typename InputAdapter<User>::gid_t    gid_t;
  typedef typename InputAdapter<User>::node_t   node_t;

  typedef Xpetra::MultiVector<
    scalar_t, lno_t, gno_t, node_t> x_mvector_t;
  typedef Xpetra::TpetraMultiVector<
    scalar_t, lno_t, gno_t, node_t> xt_mvector_t;
  typedef Xpetra::EpetraMultiVector xe_mvector_t;

  /*! Destructor
   */
  ~XpetraMultiVectorInput() { }

  /*! Constructor   
   */
  // Constructor 
  XpetraMultiVectorInput(const RCP<const User> &invector):
    invector_(invector), vector_(), map_(), env_(), base_()
  {
    vector_ = XpetraTraits<User>::convertToXpetra(invector);
    map_ = vector_->getMap();
    base_ = map_->getIndexBase();
  };

  ////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  std::string inputAdapterName()const {
    return std::string("XpetraMultiVector");}

  bool haveLocalIds() const { return true;}

  bool haveConsecutiveLocalIds(size_t &base)const {
    base = base_;
    return true;
  }

  ////////////////////////////////////////////////////
  // The MultiVectorInput interface.
  ////////////////////////////////////////////////////

  size_t getNumVectors() const {return vector_->getNumVectors();}
  
  size_t getLocalLength() const {return vector_->getLocalLength();}
  
  size_t getGlobalLength() const {return vector_->getGlobalLength();}

  lno_t getMultiVectorView(int i, const gid_t *&Ids, 
    const lid_t *&localIds, const scalar_t *&elements, 
    const scalar_t *&wgts) const
  {
    if (map_->lib() == Xpetra::UseTpetra){
      RCP<xt_mvector_t> tvector = 
        rcp_implicit_cast<xt_mvector_t>(vector_);
      ArrayRCP<const scalar_t> data = tvector->getData(i);
      elements = data->get();
    }
    else if (map_->lib() == Xpetra::UseEpetra){
      RCP<xe_mvector_t> evector = 
        rcp_implicit_cast<xe_mvector_t>(vector_);
      ArrayRCP<const scalar_t> data = evector->getData(i);
      elements = data->get();
    }
    else{
      throw std::logic_error("invalid underlying lib");
    }

    ArrayView<const gid_t> gids = map_->getNodeElementList();
    Ids = gids->getRawPtr();
    localIds = NULL;  // Implies 0 through numElements-1
    wgts = NULL; // Not implemented
  }

  ////////////////////////////////////////////////////
  // End of MatrixInput interface.
  ////////////////////////////////////////////////////

  /*! Access to xpetra vector
   */

  const RCP<const x_mvector_t> &getVector() const
  {
    return vector_;
  }

  /*! Apply a partitioning solution to the vector.
   *   Every gid that was belongs to this process must
   *   be on the list, or the Import will fail.
   */
  void applyPartitioningSolution(const User &in, User *&out,
    lno_t numIds, lno_t numParts,
    const gid_t *gid, const lid_t *lid, const lno_t *partition)
  {
    // Get an import list

    ArrayView<const gid_t> gidList(gid, numIds);
    ArrayView<const lno_t> partList(partition, numIds);
    ArrayView<const lno_t> dummyIn;
    ArrayRCP<gid_t> importList;
    ArrayRCP<int> dummyOut;
    size_t numNewRows;

    const RCP<const Comm<int> > comm = map_->getComm();

    try{
      numNewRows = convertPartitionListToImportList<gid_t, lno_t, lno_t>(
        comm, partList, gidList, dummyIn, importList, dummyOut);
    }
    catch (std::exception &e){
      Z2_THROW_ZOLTAN2_ERROR(env, e);
    }

    RCP<const User> inPtr = rcp(&in, false);
    gno_t globalNumElts = invector_->getGlobalLength();
    lno_t localNumElts = numNewRows;

    RCP<User> outPtr = XpetraTraits<User>::doImport(
     inPtr, localNumElts, importList->get(), base_);

    out = outPtr.get();
    outPtr.release();
   }

private:

  RCP<const User> invector_;
  RCP<const x_mvector_t> vector_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > map_;
  RCP<Environment> env_;    // for error messages, etc.
  lno_t base_;
};
  
}  //namespace Zoltan2
  
#endif
