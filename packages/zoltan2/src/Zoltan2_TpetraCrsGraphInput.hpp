// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_EpetraCrsGraphInput.hpp

    \brief An input adapter for a Epetra::CrsGraph.
*/

#ifndef _ZOLTAN2_TPETRACRSGRAPHINPUT_HPP_
#define _ZOLTAN2_TPETRACRSGRAPHINPUT_HPP_

#include <Tpetra_CrsGraph.hpp>
#include <Zoltan2_XpetraCrsGraphInput.hpp>
#include <Zoltan2_Standards.hpp>

namespace Zoltan2 {

/*! Zoltan2::TpetraCrsGraphInput
    \brief Provides access for Zoltan2 to Tpetra::CrsGraph data plus weights.

    TpetraCrsGraphInput provides access for Zoltan2 to a
    Tpetra::CrsGraph using an
    interface that is the same for all graph input.  Most of the methods used by
    Zoltan2 are in the XpetraGraphInput parent class.  The methods used
    by the Zoltan2 caller are in this class.

    The template parameter is the weight type.
*/

template<typename LNO, typename GNO, typename LID=LNO, typename GID=GNO, 
  typename Node=Kokkos::DefaultNode::DefaultNodeType>
class TpetraCrsGraphInput : 
  public XpetraCrsGraphInput<LNO, GNO, LID, GID, Node>{
private:
  typedef Tpetra::CrsGraph<LNO, GNO, Node> crsGraph;
  RCP<const crsGraph > ingraph_;

public:
  /*! Name of input adapter type.
   */
  std::string inputAdapterName()const {return std::string("TpetraCrsGraph");}

  /*! Destructor
   */
  ~TpetraCrsGraphInput() { }

  /*! Constructor
   */
  TpetraCrsGraphInput(const RCP<const crsGraph> graph): 
      XpetraCrsGraphInput<LNO, GNO, LID, GID, Node>(
        Teuchos::rcp(new Xpetra::TpetraCrsGraph<LNO, GNO, Node>(
          Teuchos::rcp_const_cast<crsGraph>(graph)))) 
  {
    ingraph_ = graph;
  }


  /*! Access to the graph that instantiated the adapter
   */
  RCP<const crsGraph> getGraph()
  { 
    return ingraph_;
  }
};

} // namespace

#endif
