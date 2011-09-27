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
#include <Zoltan2_Standards.hpp>

namespace Zoltan2 {

/*! Zoltan2::TpetraCrsMatrixInput
    \brief This objects provides access for Zoltan2 to Tpetra::CrsMatrix data.

*/

CONSISTENT_CLASS_TEMPLATE_LINE
class TpetraCrsMatrixInput : public
         XpetraCrsMatrixInput<CONSISTENT_TEMPLATE_PARAMS>
{
private:

public:

  /*! Default constructor
   */
  TpetraCrsMatrixInput() { }

  /*! Constructor
   */
  TpetraCrsMatrixInput(RCP<Tpetra::CrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> >
        matrix): XpetraCrsMatrixInput<CONSISTENT_TEMPLATE_PARAMS>(
        rcp(new Xpetra::TpetraCrsMatrix<CONSISTENT_TRILINOS_TEMPLATE_PARAMS> (matrix)))
// KDDKDD matrix)
// TODO:  Siva, please review this change.  I was following a Zoltan2 email 
// TODO:  from Erik on August 5.
// TODO:  Since the XpetraCrsMatrixInput constructors need 
// TODO:  Xpetra::TpetraCrsMatrix, I am doing the conversion here.
  {
  }
};
  
  
}  //namespace Zoltan2
  
#endif
