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
    \brief This provides access for Zoltan2 to Epetra::CrsGraph data and weights.

    EpetraCrsGraphInput provides access for Zoltan2 to an Epetra::CrsGraph and
    vertex and edges weights.  Most of the methods used by Zoltan2 are in the 
    XpetraGraphInput parent class.  The methods used by the Zoltan2 caller 
    to initialize the objects are in this class.

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

      TODO: Perhaps the caller would like to to provide a matrix instead of
            graph, the nonzeros of which are the edge and vertex weights.
      
      \param vertexWeights are given in vertex local number order
      \param edgeWeights when queried, the Epetra::CrsGraph gives edges in a 
               certain order.  The edgeWeights follow this order, omitting
               any self edges.
   */
  EpetraCrsGraphInput(Teuchos::RCP<Epetra::CrsGraph> graph,
                      Teuchos::ArrayRCP<Scalar> vertexWeights,
                      Teuchos::ArrayRCP<Scalar> edgeWeights){

     try{
       this->setGraph(graph,vertexWeights,edgeWeights);
     } catch(const std::exception &e)
         Z2_THROW_ZOLTAN2_ERROR(*_env, e)
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

    try{
      this->setXpetraCrsGraph(xgraph, vwgt, ewgt);
    } catch(const std::exception &e)
         Z2_THROW_ZOLTAN2_ERROR(*_env, e)


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
