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

    RCP<const Epetra_CrsMatrix> inmatrix_;
public:

  // KDDKDD This may not be right.  May need to specify double, int, int.
  typedef typename double scalar_t;
  typedef typename int gno_t;
  typedef typename int lno_t;
  typedef typename gno_t gid_t;
  typedef typename lno_t lid_t;
  typedef typename node_t;

  std::string inputAdapterName()const {return std::string("EpetraCrsMatrix");}

  ~EpetraCrsMatrixInput() { }

  /*! Constructor
   */
  EpetraCrsMatrixInput(const RCP<const Epetra_CrsMatrix> matrix):
    XpetraCrsMatrixInput<double, int, int> (
      Teuchos::rcp(new Xpetra::EpetraCrsMatrix(
        Teuchos::rcp_const_cast<Epetra_CrsMatrix>(matrix)))) 
  {
    inmatrix_ = matrix;
  }

  /*! Access to matrix that instantiated adapter
   */
  RCP<const Epetra_CrsMatrix> getMatrix() const 
  { 
    return inmatrix_;
  }

};

} //namespace Zoltan2

#endif
