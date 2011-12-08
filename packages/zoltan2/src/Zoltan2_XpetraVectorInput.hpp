// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraVectorInput.hpp

    \brief An input adapter for a Xpetra::Vector, Epetra_Vector and 
            Tpetra::Vector.
*/

#ifndef _ZOLTAN2_XPETRAVECTORINPUT_HPP_
#define _ZOLTAN2_XPETRAVECTORINPUT_HPP_

#include <Zoltan2_VectorInput.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_Util.hpp>

#include <Xpetra_EpetraVector.hpp>
#include <Xpetra_TpetraVector.hpp>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////////
/*! Zoltan2::XpetraVectorInput
    \brief Provides access for Zoltan2 to an Xpetra::Vector .

    The template parameter is the user's input object: 
     Epetra_Vector
     Tpetra::Vector
     Xpetra::Vector
*/

template <typename User>
class XpetraVectorInput : public VectorInput<User> {
public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef VectorInput<User>       base_adapter_t;

  typedef Xpetra::Vector<
    scalar_t, lno_t, gno_t, node_t> x_vector_t;
  typedef Xpetra::TpetraVector<
    scalar_t, lno_t, gno_t, node_t> xt_vector_t;
  typedef Xpetra::EpetraVector xe_vector_t;

  /*! Destructor
   */
  ~XpetraVectorInput() { }

  /*! Constructor   
   */
  // Constructor 
  XpetraVectorInput(const RCP<const User> &invector):
    invector_(invector), 
    vector_(),
    map_(),
    env_(),
    base_()
  {
    vector_ = XpetraTraits<User>::convertToXpetra(invector);
    map_ = vector_->getMap();
    base_ = map_->getIndexBase();
  };

  /*! Access to xpetra vector
   */

  const RCP<const x_vector_t> &getVector() const
  {
    return vector_;
  }

  ////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  std::string inputAdapterName()const {
    return std::string("XpetraVector");}

  ////////////////////////////////////////////////////
  // The VectorInput interface.
  ////////////////////////////////////////////////////

  size_t getLocalLength() const {return vector_->getLocalLength();}
  
  size_t getGlobalLength() const {return vector_->getGlobalLength();}

  size_t getVectorView(const gid_t *&Ids, const scalar_t *&elements, 
    const scalar_t *&wgts) const
  {
    elements = NULL;
    const x_vector_t *vec =  vector_.get();

    if (map_->lib() == Xpetra::UseTpetra){
      const xt_vector_t *tvector = dynamic_cast<const xt_vector_t *>(vec);

      if (tvector->getLocalLength() > 0){
        // getData hangs if vector length is 0
        ArrayRCP<const scalar_t> data = tvector->getData(0);
        elements = data.get();
      }
    }
    else if (map_->lib() == Xpetra::UseEpetra){
      const xe_vector_t *evector = dynamic_cast<const xe_vector_t *>(vec);
        
      if (evector->getLocalLength() > 0){
        // getData hangs if vector length is 0
        ArrayRCP<const double> data = evector->getData(0);
        elements = data.get();
      }
    }
    else{
      throw std::logic_error("invalid underlying lib");
    }

    ArrayView<const gid_t> gids = map_->getNodeElementList();
    Ids = gids.getRawPtr();
    wgts = NULL; // Not implemented

    return getLocalLength();
  }

  /*! Apply a partitioning solution to the vector.
   *   Every gid that was belongs to this process must
   *   be on the list, or the Import will fail.
   */
  size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<gid_t, lno_t> &solution)
  { 
    // Get an import list

    ArrayRCP<gid_t> &gidList = solution.getGidsRCPConst();
    ArrayRCP<size_t> &partList = solution.getPartsRCPConst();
    ArrayRCP<lno_t> dummyIn;
    ArrayRCP<gid_t> importList;
    ArrayRCP<lno_t> dummyOut;
    size_t numNewRows;

    const RCP<const Comm<int> > comm = map_->getComm(); 

    try{
      numNewRows = convertPartListToImportList<gid_t, lno_t, lno_t>(
        *comm, partList, gidList, dummyIn, importList, dummyOut);
    }
    Z2_FORWARD_EXCEPTIONS;

    RCP<const User> inPtr = rcp(&in, false);
    lno_t localNumElts = numNewRows;

    RCP<const User> outPtr = XpetraTraits<User>::doMigration(
     inPtr, localNumElts, importList.get());

    out = const_cast<User *>(outPtr.get());
    outPtr.release();
    return numNewRows;
  }

private:

  RCP<const User> invector_;
  RCP<const x_vector_t> vector_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > map_;
  RCP<Environment> env_;
  lno_t base_;

};
  
}  //namespace Zoltan2
  
#endif
