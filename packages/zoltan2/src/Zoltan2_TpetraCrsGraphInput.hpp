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

#include <Zoltan2_XpetraGraphInput.hpp>
#include <Tpetra_CrsGraph.h>

namespace Zoltan2 {

/*! Zoltan2::TpetraCrsGraphInput
    \brief TpetraCrsGraphInput provides access for Zoltan2 to Tpetra::CrsGraph data.

    TpetraCrsGraphInput provides access for Zoltan2 to an Tpetra::CrsGraph using an
    interface that is the same for all graph input.  Most of the methods used by
    Zoltan2 objects are in the XpetraGraphInput parent class.  The methods used
    by the Zoltan2 caller are in this class.

    The template parameter is the weight type.

    TODO: include NODE
*/

template <typename Scalar, typename LNO, typename GNO>
class TpetraCrsGraphInput : public XpetraGraphInput<Scalar, LNO, GNO>{
private:

public:

  /*! Default constructor
   */
  TpetraCrsGraphInput(){
    _inputComplete = false;
    _adapterProcessingComplete= false;
  }

  /*! Constructor with graph only
   */
  TpetraCrsGraphInput(Teuchos::RCP<const Tpetra::CrsGraph> graph){ this->setGraph(graph);}

  /*! Constructor with weights
   */
  TpetraCrsGraphInput(Teuchos::RCP<Tpetra::CrsGraph> graph,
                      Teuchos::ArrayRCP<Scalar> vertexWeights,
                      Teuchos::ArrayRCP<Scalar> edgeWeights){

     this->setGraph(graph,vertexWeights,edgeWeights);
  }

  /*! Post construction update
   */

  void setGraph(Teuchos::RCP<Tpetra::CrsGraph> graph){
    _inputComplete = false;
    _adapterProcessingComplete = false;
    Teuchos::ENull null;
    Teuchos::ArrayRCP<Scalar> vwgt(null);
    Teuchos::ArrayRCP<Scalar> ewgt(null);

    Xpetra::TpetraCrsGraph<Scalar, LNO, GNO> xgraph(graph);

    this->setXpetraCrsGraph(xgraph, vwgt, ewgt);

    _inputComplete = true;
    _adapterProcessingComplete = true;
  }

  /*! Post construction update
   */
  void setGraph(Teuchos::RCP<Tpetra::CrsGraph> graph);
                Teuchos::ArrayRCP<Scalar> vertexWeights,
                Teuchos::ArrayRCP<Scalar> edgeWeights){
    _inputComplete = false;
    _adapterProcessingComplete = false;
    Xpetra::TpetraCrsGraph<Scalar, LNO, GNO> xgraph(graph);

    this->setXpetraCrsGraph(xgraph, vertexWeights, edgeWeights);

    _inputComplete = true;
    _adapterProcessingComplete = true;
  }
                
};
  
  
}  //namespace Zoltan2
  
#endif
