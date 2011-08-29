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
    \brief This objects provides access for Zoltan2 to Tpetra::CrsGraph data plus weights.

    TpetraCrsGraphInput provides access for Zoltan2 to an Tpetra::CrsGraph using an
    interface that is the same for all graph input.  Most of the methods used by
    Zoltan2 are in the XpetraGraphInput parent class.  The methods used
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

     try{
       this->setGraph(graph,vertexWeights,edgeWeights);
     } catch(const std::exception &e)
         Z2_THROW_ZOLTAN2_ERROR(*_env, e)
  }

  /*! Post construction update
   */

  void setGraph(Teuchos::RCP<Tpetra::CrsGraph> graph){
    _inputComplete = false;
    _adapterProcessingComplete = false;
    Teuchos::ENull null;
    Teuchos::ArrayRCP<Scalar> vwgt(null);
    Teuchos::ArrayRCP<Scalar> ewgt(null);

    try{
      Xpetra::TpetraCrsGraph<Scalar, LNO, GNO> xgraph(graph);
    } catch(const std::exception &e)
        Z2_THROW_ZOLTAN2_ERROR(*_env, e)

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

    try{
      this->setXpetraCrsGraph(xgraph, vertexWeights, edgeWeights);
    } catch(const std::exception &e)
        Z2_THROW_ZOLTAN2_ERROR(*_env, e)

    _inputComplete = true;
    _adapterProcessingComplete = true;
  }
                
};
  
  
}  //namespace Zoltan2
  
#endif
