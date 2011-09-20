// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_EpetraCrsMatrixInput.hpp

    \brief An input adapter for a Epetra::CrsMatrix.

    \author Siva Rajamanickam
*/

#ifndef _ZOLTAN2_TPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_TPETRACRSMATRIXINPUT_HPP_

#include <Tpetra_CrsMatrix.h>
#include <Xpetra_CrsMatrix.hpp>

namespace Zoltan2 {

/*! Zoltan2::TpetraCrsMatrixInput
    \brief This objects provides access for Zoltan2 to Tpetra::CrsMatrix data.

*/

CONSISTENT_CLASS_TEMPLATE_LINE
class TpetraCrsMatrixInput : public MatrixInput {
private:
      //Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LNO, GNO>> _xmatrix;
      Teuchos::RCP<Xpetra::XpetraCrsMatrix<Scalar, LNO, GNO>> _xmatrix;

public:

  /*! Default constructor
   */
  TpetraCrsMatrixInput(){
    //_inputComplete = false;
    //_adapterProcessingComplete= false;
  }

  /*! Constructor with matrix only
   */
  TpetraCrsMatrixInput(Teuchos::RCP<const Tpetra::CrsMatrix> matrix){
      this->setMatrix(matrix);
  }

  /*! Constructor with weights
    TODO: Do the weight have to be arrays?
   */
  /*TpetraCrsGraphInput(Teuchos::RCP<Tpetra::CrsGraph> graph,
                      Teuchos::ArrayRCP<Scalar> vertexWeights,
                      Teuchos::ArrayRCP<Scalar> edgeWeights){

     try{
       this->setGraph(graph,vertexWeights,edgeWeights);
     } catch(const std::exception &e)
         Z2_THROW_ZOLTAN2_ERROR(*_env, e)
  }*/

  /*! Post construction update
    // TODO: Should be private
   */

  void setMatrix(Teuchos::RCP<Tpetra::CrsMatrix> matrix){
    _inputComplete = false;
    _adapterProcessingComplete = false;
    //Teuchos::ENull null;
    //Teuchos::ArrayRCP<Scalar> vwgt(null);
    //Teuchos::ArrayRCP<Scalar> ewgt(null);

    try{
      //_xmatrix = Teuchos::RCP(new Xpetra::TpetraCrsMatrix<Scalar, LNO, GNO>
                                    //(matrix));
      _xmatrix = Teuchos::RCP(new Xpetra::XpetraCrsMatrix<Scalar, LNO, GNO>
                                    (matrix));
    } catch(const std::exception &e)
        Z2_THROW_ZOLTAN2_ERROR(*_env, e)

    _inputComplete = true;
    _adapterProcessingComplete = true;
  }

  /*! Post construction update
    TODO: Do the weight have to be arrays?
   */
  /*void setGraph(Teuchos::RCP<Tpetra::CrsGraph> graph);
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
  }*/
                
  /*! Returns the global number of rows in the matrix.
   */
  GNO getGlobalNumRows(){
      _xmatrix->getGlobalNumRows();
  }

  /*! Returns the global number columns in the matrix.
   */
  virtual GNO getGlobalNumCols(){
      _xmatrix->getGlobalNumCols();
  }

  /*! Returns the number rows on this process.
   */
  virtual LNO getNodeNumRows(){
      _xmatrix->getNodeNumRows();
  }

  /*! Returns the number edges on this process.
   */
  virtual LNO getNodeNumCols(){
      _xmatrix->getNodeNumCols();
  }

  //TODO: Add just the required functions.
};
  
  
}  //namespace Zoltan2
  
#endif
