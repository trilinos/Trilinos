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

#ifndef _ZOLTAN2_EPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_EPETRACRSMATRIXINPUT_HPP_

#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_TemplateMacros.hpp>

#include <Epetra_CrsMatrix.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Teuchos_RCP.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

namespace Zoltan2 {

/*! Zoltan2::EpetraCrsMatrixInput
    \brief This objects provides access for Zoltan2 to Epetra_CrsMatrix data.

*/

class EpetraCrsMatrixInput : public XpetraCrsMatrixInput<double, int, int>
{
private:

public:

  /*! Default constructor
   */
  EpetraCrsMatrixInput()
  {
  }

  /*! Constructor
   */
  //TODO EpetraCrsMatrixInput(RCP<const Epetra_CrsMatrix >
  EpetraCrsMatrixInput(RCP<Epetra_CrsMatrix matrix) :
      XpetraCrsMatrixInput<double, int, int>(rcp(new Xpetra::EpetraCrsMatrix
                                                     (matrix)))
  {
      HELLO;
  }

};
  
  
}  //namespace Zoltan2
  
#endif
