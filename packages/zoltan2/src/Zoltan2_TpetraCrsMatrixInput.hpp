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

public:

  typedef typename InputAdapter<User>::scalar_t scalar_t;
  typedef typename InputAdapter<User>::lno_t    lno_t;
  typedef typename InputAdapter<User>::gno_t    gno_t;
  typedef typename InputAdapter<User>::lid_t    lid_t;
  typedef typename InputAdapter<User>::gid_t    gid_t;
  typedef typename InputAdapter<User>::node_t   node_t;
  typedef Xpetra::TpetraCrsMatrix<scalar_t, lno_t, gno_t, node_t> xtmatrix_t;
  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> xmatrix_t;

  /*! Name of input adapter type.
   */
  std::string inputAdapterName()const {return std::string("TpetraCrsMatrix");}

  /*! Destructor
   */
  ~TpetraCrsMatrixInput() { }

  /*! Constructor 
   */
  TpetraCrsMatrixInput(const RCP<const User> &matrix):
    XpetraCrsMatrixInput<User>(rcp_implicit_cast<const xmatrix_t>(
      rcp(new xtmatrix_t(rcp_const_cast<User>(matrix)))))

  {
    inmatrix_ = matrix;
  }

private:

  RCP<const User> inmatrix_;
};
} // namespace

#endif
