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
  TpetraCrsGraphInput(): XpetraCrsGraphInput<CONSISTENT_TEMPLATE_PARAMS>() {}

  /*! Constructor with graph only
   */
  TpetraCrsGraphInput(
    RCP<Tpetra::CrsGraph<LNO, GNO, Node> > graph):
      XpetraCrsGraphInput<CONSISTENT_TEMPLATE_PARAMS>(graph) {}

  /*! Constructor with weights
   */
  TpetraCrsGraphInput(
    RCP<Tpetra::CrsGraph<LNO, GNO, Node> > graph,
    ArrayRCP<Scalar> vertexWeights,
    ArrayRCP<Scalar> edgeWeights):
      XpetraCrsGraphInput<CONSISTENT_TEMPLATE_PARAMS>(
        graph, vertexWeights, edgeWeights) {}
};

} // namespace

#endif
