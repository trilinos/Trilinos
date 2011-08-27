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

#ifndef _ZOLTAN2_EPETRACRSGRAPHINPUT_HPP_
#define _ZOLTAN2_EPETRACRSGRAPHINPUT_HPP_

#include <Zoltan2_XpetraGraphInput.hpp>
#include <Epetra_CrsGraph.h>

namespace Zoltan2 {

/*! Zoltan2::EpetraCrsGraphInput
    \brief EpetraCrsGraphInput provides access for Zoltan2 to Epetra::CrsGraph data.

    EpetraCrsGraphInput provides access for Zoltan2 to an Epetra::CrsGraph using an
    interface that is the same for all graph input.  Most of the methods used by
    Zoltan2 objects are in the XpetraGraphInput parent class.  The methods used
    by the Zoltan2 caller are in this class.

    The template parameter is the weight type.  Epetra local and global IDs are ints.
*/

template <Scalar>
class EpetraCrsGraphInput : public XpetraGraphInput<Scalar, int, int>{
private:

public:

  /*! Default constructor
   */
  EpetraCrsGraphInput(){
    _inputComplete = false;
    _adapterProcessingComplete= false;
  }

  /*! Constructor with graph only
   */
  EpetraCrsGraphInput(Teuchos::RCP<const Epetra::CrsGraph> graph){ this->setGraph(graph);}

  /*! Constructor with weights
   */
  EpetraCrsGraphInput(Teuchos::RCP<Epetra::CrsGraph> graph,
                      Teuchos::ArrayRCP<Scalar> vertexWeights,
                      Teuchos::ArrayRCP<Scalar> edgeWeights){

     this->setGraph(graph,vertexWeights,edgeWeights);
  }

  /*! Post construction update
   */

  void setGraph(Teuchos::RCP<Epetra::CrsGraph> graph){
    _inputComplete = false;
    _adapterProcessingComplete = false;
    Teuchos::ENull null;
    Teuchos::ArrayRCP<Scalar> vwgt(null);
    Teuchos::ArrayRCP<Scalar> ewgt(null);

    Xpetra::EpetraCrsGraph xgraph(graph);

    this->setXpetraCrsGraph(xgraph, vwgt, ewgt);

    _inputComplete = true;
    _adapterProcessingComplete = true;
  }

  /*! Post construction update
   */
  void setGraph(Teuchos::RCP<Epetra::CrsGraph> graph);
                Teuchos::ArrayRCP<Scalar> vertexWeights,
                Teuchos::ArrayRCP<Scalar> edgeWeights){
    _inputComplete = false;
    _adapterProcessingComplete = false;
    Xpetra::EpetraCrsGraph xgraph(graph);

    this->setXpetraCrsGraph(xgraph, vertexWeights, edgeWeights);

    _inputComplete = true;
    _adapterProcessingComplete = true;
  }
                
};
  
  
}  //namespace Zoltan2
  
#endif
