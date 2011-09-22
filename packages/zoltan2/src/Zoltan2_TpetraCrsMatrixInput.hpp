// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_TpetraCrsMatrixInput.hpp

    \brief An input adapter for a Tpetra::CrsMatrix.

    \author Siva Rajamanickam
*/

#ifndef _ZOLTAN2_TPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_TPETRACRSMATRIXINPUT_HPP_

#include <Zoltan2_XpetraCrsMatrixInput.hpp>
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

CONSISTENT_TRILINOS_CLASS_TEMPLATE_LINE
class TpetraCrsMatrixInput : public
         XpetraCrsMatrixInput<CONSISTENT_TRILINOS_TEMPLATE_PARAMS>
{
private:

public:

  /*! Default constructor
   */
  TpetraCrsMatrixInput()
  {
  }

  /*! Constructor
   */
  //TODO TpetraCrsMatrixInput(RCP<const Tpetra::CrsMatrix<Scalar,LNO,GNO,Node> >
  TpetraCrsMatrixInput(RCP<Tpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> >
        matrix): XpetraCrsMatrixInput<CONSISTENT_TRILINOS_TEMPLATE_PARAMS>(
        rcp(new Xpetra::TpetraCrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS>
         (matrix)))
  {
      HELLO;
  }

};
  
  
}  //namespace Zoltan2
  
#endif
