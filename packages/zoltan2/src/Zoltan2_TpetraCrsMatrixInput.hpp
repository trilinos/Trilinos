// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_TpetraCrsMatrixInput.hpp

    \brief An input adapter for a Epetra::CrsMatrix.

    \author Siva Rajamanickam
*/

#ifndef _ZOLTAN2_TPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_TPETRACRSMATRIXINPUT_HPP_

#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_TemplateMacros.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Teuchos_RCP.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace Zoltan2 {

/*! Zoltan2::TpetraCrsMatrixInput
    \brief This objects provides access for Zoltan2 to Tpetra::CrsMatrix data.

*/

CONSISTENT_CLASS_TEMPLATE_LINE
class TpetraCrsMatrixInput : public MatrixInput<CONSISTENT_TEMPLATE_PARAMS>
{
private:
      Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LNO, GNO> > _xmatrix;

public:

  /*! Default constructor
   */
  TpetraCrsMatrixInput()
  {
    //_inputComplete = false;
    //_adapterProcessingComplete= false;
  }

  /*! Constructor with matrix only
   */
  // TODO TpetraCrsMatrixInput(RCP<const Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> >
  TpetraCrsMatrixInput(RCP<Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> >
                            matrix)
  {
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

  void setMatrix(RCP<Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> > matrix)
  {
    //_inputComplete = false;
    //_adapterProcessingComplete = false;
    //Teuchos::ENull null;
    //Teuchos::ArrayRCP<Scalar> vwgt(null);
    //Teuchos::ArrayRCP<Scalar> ewgt(null);

    try{
      _xmatrix = Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LNO, GNO> >(new
                     Xpetra::TpetraCrsMatrix<Scalar, LNO, GNO> (matrix));
    } catch(const std::exception &e)
    {
        throw(e);
    }

    //_inputComplete = true;
    //_adapterProcessingComplete = true;
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
  GNO getGlobalNumRows() const{
      return(_xmatrix->getGlobalNumRows());
  }

  /*! Returns the global number columns in the matrix.
   */
  GNO getGlobalNumCols() const{
      return(_xmatrix->getGlobalNumCols());
  }

  /*! Returns the number rows on this process.
   */
  LNO getNodeNumRows() const{
      return(_xmatrix->getNodeNumRows());
  }

  /*! Returns the number edges on this process.
   */
  LNO getNodeNumCols() const{
      return (_xmatrix->getNodeNumCols());
  }

  const Teuchos::Comm<int> &getTeuchosComm(){
      //TODO: Xpetra has to add getComm
      //return(*(_xmatrix->getComm()));
  }

  const MPI_Comm &getMpiComm(){
#ifdef HAVE_MPI
      //TODO: Xpetra has to add getComm
      //return (*((_xmatrix->getComm()).getRawMpiComm()));
#else
      // TODO: What is the valid value here ?
      return 0;
#endif
  }
  //TODO: Add just the required functions.
};
  
  
}  //namespace Zoltan2
  
#endif
