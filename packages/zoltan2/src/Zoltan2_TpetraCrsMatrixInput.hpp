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

template <typename User>
class TpetraCrsMatrixInput : public XpetraCrsMatrixInput<User>{
private:

  typedef typename XpetraCrsMatrixInput<User>::scalar_t  scalar_t;
  typedef typename XpetraCrsMatrixInput<User>::gno_t     gno_t;
  typedef typename XpetraCrsMatrixInput<User>::lno_t     lno_t;
  typedef typename XpetraCrsMatrixInput<User>::gid_t     gid_t;
  typedef typename XpetraCrsMatrixInput<User>::lid_t     lid_t;
  typedef typename XpetraCrsMatrixInput<User>::node_t    node_t;
  typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> crsMatrix;

  RCP<const crsMatrix > _inmatrix;

public:

  std::string inputAdapterName()const {return std::string("TpetraCrsMatrix");}

  ~TpetraCrsMatrixInput() { }

  /*! Constructor 
   */
  TpetraCrsMatrixInput(const RCP<const crsMatrix> matrix):
    XpetraCrsMatrixInput<User>(
      Teuchos::rcp(new Xpetra::TpetraCrsMatrix<scalar_t, lno_t, gno_t, node_t>(
        Teuchos::rcp_const_cast<crsMatrix>(matrix))))
    
  {
    _inmatrix = matrix;
  }

  /*! Access to matrix that instantiated adapter
   */

  RCP<const crsMatrix> getMatrix() const
  { 
    return _inmatrix;
  }
};
} // namespace

#endif
