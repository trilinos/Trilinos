// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_EpetraCrsMatrixInput.hpp

    \brief An input adapter for a Epetra_CrsMatrix.
*/

#ifndef _ZOLTAN2_EPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_EPETRACRSMATRIXINPUT_HPP_

#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Epetra_CrsMatrix.h>

namespace Zoltan2 {

/*! Zoltan2::EpetraCrsMatrixInput
    \brief Provides access for Zoltan2 to Epetra_CrsMatrix data.
*/

class EpetraCrsMatrixInput : public XpetraCrsMatrixInput<double, int, int>{
private:

public:

  std::string inputAdapterName()const {return std::string("EpetraCrsMatrix");}

  ~EpetraCrsMatrixInput() { }

  /*! Constructor with matrix only
   */
  EpetraCrsMatrixInput(RCP<Epetra_CrsMatrix> matrix):
    XpetraCrsMatrixInput<double, int, int> (matrix) {}
};

} //namespace Zoltan2

#endif
