// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsMatrixInput.hpp

    \brief An input adapter for a Xpetra::CrsMatrix.

    \author Siva Rajamanickam
*/

#ifndef _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_

#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_TemplateMacros.hpp>

#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Teuchos_RCP.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace Zoltan2 {

/*! Zoltan2::XpetraCrsMatrixInput
    \brief This objects provides access for Zoltan2 to Xpetra::CrsMatrix data.

*/

CONSISTENT_TRILINOS_CLASS_TEMPLATE_LINE
class XpetraCrsMatrixInput : public
             MatrixInput<CONSISTENT_TRILINOS_TEMPLATE_PARAMS>
{
private:
      Teuchos::RCP<Xpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > _xmatrix;

public:

  /*! Default constructor
   */
  XpetraCrsMatrixInput()
  {
  }

  /*! Constructor
   */
  //TODO i/p param should be RCP<const Xpetra::CrsMatrix... >
  XpetraCrsMatrixInput(RCP<Xpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> >
                            matrix) :
                            _xmatrix(matrix)
  {
      HELLO;
  }

  //TODO i/p param should be RCP<const Xpetra::TpetraCrsMatrix... >
  XpetraCrsMatrixInput(RCP<Xpetra::TpetraCrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> >
                            matrix)
  {
      HELLO;
      _xmatrix = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix
                             <CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > (matrix);
  }

  //TODO i/p param should be RCP<const Xpetra::TpetraCrsMatrix... >
  XpetraCrsMatrixInput(RCP<Xpetra::EpetraCrsMatrix> matrix)
  {
      HELLO;
      _xmatrix = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix
                             <CONSISTENT_TRILINOS_TEMPLATE_PARAMS> > (matrix);
  }

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
  LNO getLocalNumRows() const{
      return(_xmatrix->getNodeNumRows());
  }

  /*! Returns the number edges on this process.
   */
  LNO getLocalNumCols() const{
      return (_xmatrix->getNodeNumCols());
  }

  const Teuchos::Comm<int> &getTeuchosComm(){
      //TODO: Need Xpetra support for getComm
      //return(*(_xmatrix->getComm()));
  }

  const MPI_Comm &getMpiComm(){
#ifdef HAVE_MPI
      //TODO: Need Xpetra support for getComm
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
