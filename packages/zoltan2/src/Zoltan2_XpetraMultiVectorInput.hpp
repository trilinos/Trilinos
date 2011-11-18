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
#include <Zoltan2_Util.hpp>

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

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::lid_t    lid_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;

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

  size_t getMultiVectorView(int i, const gid_t *&Ids, 
    const lid_t *&localIds, const scalar_t *&elements, 
    const scalar_t *&wgts) const
  {
    elements = NULL;
    if (map_->lib() == Xpetra::UseTpetra){
      const xt_mvector_t *tvector = 
        dynamic_cast<const xt_mvector_t *>(vector_.get());
       
      if (tvector->getLocalLength() > 0){
        ArrayRCP<const scalar_t> data = tvector->getData(i);
        elements = data.get();
      }
    }
    else if (map_->lib() == Xpetra::UseEpetra){
      const xe_mvector_t *evector = 
        dynamic_cast<const xe_mvector_t *>(vector_.get());
        
      if (evector->getLocalLength() > 0){
        ArrayRCP<const double> data = evector->getData(i);
        elements = data.get();
      }
    }
    else{
      throw std::logic_error("invalid underlying lib");
    }

    ArrayView<const gid_t> gids = map_->getNodeElementList();
    Ids = gids.getRawPtr();
    localIds = NULL;  // Implies 0 through numElements-1
    wgts = NULL; // Not implemented
    return gids.size();
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
  size_t applyPartitioningSolution(const User &in, User *&out,
    size_t numParts, size_t numIds,
    const gid_t *gid, const lid_t *lid, const size_t *partition)
  {
    ArrayView<const gid_t> gidList(gid, numIds);
    ArrayView<const size_t> partList(partition, numIds);
    ArrayView<const lno_t> dummyIn;
    ArrayRCP<gid_t> importList;
    ArrayRCP<lno_t> dummyOut;
    size_t numNewRows;

    const RCP<const Comm<int> > comm = map_->getComm();

    try{
      // Get an import list
      numNewRows = convertPartitionListToImportList<gid_t, lno_t, lno_t>(
        *comm, partList, gidList, dummyIn, importList, dummyOut);
    }
    catch (std::exception &e){
      Z2_THROW_ZOLTAN2_ERROR(env, e);
    }

    RCP<const User> inPtr = rcp(&in, false);
    lno_t localNumElts = numNewRows;

    RCP<const User> outPtr = XpetraTraits<User>::doMigration(
     inPtr, localNumElts, importList.get(), base_);

    out = const_cast<User *>(outPtr.get());
    outPtr.release();
    return numNewRows;
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
