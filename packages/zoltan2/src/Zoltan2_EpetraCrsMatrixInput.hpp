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


#include <Epetra_CrsMatrix.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>

namespace Zoltan2 {

/*! Zoltan2::EpetraCrsMatrixInput
    \brief This objects provides access for Zoltan2 to Epetra_CrsMatrix data.

*/

class EpetraCrsMatrixInput : public XpetraCrsMatrixInput<double, int, int>
{
private:

public:

  /*! Default constructor   TODO - remove?
   */
  EpetraCrsMatrixInput() { }

  /*! Constructor
   */

  EpetraCrsMatrixInput(RCP<Epetra_CrsMatrix> matrix) :
      XpetraCrsMatrixInput<double, int, int>(matrix){}
};
  
  
}  //namespace Zoltan2
  
#endif
