// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_TpetraCrsMatrixInput.hpp

    \brief An input adapter for a Tpetra::CrsMatrix.
*/

#ifndef _ZOLTAN2_TPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_TPETRACRSMATRIXINPUT_HPP_

#include <Tpetra_CrsMatrix.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>

namespace Zoltan2 {

/*! Zoltan2::TpetraCrsMatrixInput
    \brief Provides access for Zoltan2 to Tpetra::CrsMatrix data. 
*/

CONSISTENT_CLASS_TEMPLATE_LINE
class TpetraCrsMatrixInput : 
  public XpetraCrsMatrixInput<CONSISTENT_TEMPLATE_PARAMS>{
private:

public:
  std::string inputAdapterName()const {return std::string("TpetraCrsMatrix");}

  ~TpetraCrsMatrixInput() { }

  /*! Default constructor   TODO - remove?
   */
  TpetraCrsMatrixInput(): XpetraCrsMatrixInput<CONSISTENT_TEMPLATE_PARAMS>() {}

  /*! Constructor with matrix only
   */
  TpetraCrsMatrixInput(
    RCP<Tpetra::CrsMatrix<Scalar, LNO, GNO, Node> > matrix):
      XpetraCrsMatrixInput<CONSISTENT_TEMPLATE_PARAMS>(matrix) {}

};
} // namespace

#endif
