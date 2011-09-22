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

#include <Xpetra_TpetraCrsGraph.hpp>
#include <Zoltan2_GraphInput.hpp>
#include <Teuchos_RCP.hpp>

namespace Zoltan2 {

/*! Zoltan2::TpetraCrsGraphInput
    \brief Provides access for Zoltan2 to Tpetra::CrsGraph data plus weights.

    TpetraCrsGraphInput provides access for Zoltan2 to a
    Tpetra::CrsGraph using an
    interface that is the same for all graph input.  Most of the methods used by
    Zoltan2 are in the XpetraGraphInput parent class.  The methods used
    by the Zoltan2 caller are in this class.

    The template parameter is the weight type.

    TODO: include NODE
*/

CONSISTENT_CLASS_TEMPLATE_LINE
class TpetraCrsGraphInput : 
  public XpetraCrsGraphInput<CONSISTENT_TEMPLATE_PARAMS>{
private:

public:
  std::string inputAdapterName()const {return std::string("TpetraCrsGraph");}

  //~TpetraCrsGraphInput() { }

  /*! Default constructor 
   */
  TpetraCrsGraphInput():XpetraCrsGraphInput() {}

  /*! Constructor with graph only
   */
  TpetraCrsGraphInput(
    Teuchos::RCP<Tpetra::CrsGraph<CONSISTENT_TEMPLATE_PARAMS> > graph):
      XpetraCrsGraphInput(graph) {}

  /*! Constructor with weights
   */
  TpetraCrsGraphInput(
    Teuchos::RCP<Tpetra::CrsGraph<CONSISTENT_TEMPLATE_PARAMS> > graph,
    Teuchos::ArrayRCP<Scalar> vertexWeights,
    Teuchos::ArrayRCP<Scalar> edgeWeights):
      XpetraCrsGraphInput(graph, vertexWeights, edgeWeights) {}
};

} // namespace

#endif
