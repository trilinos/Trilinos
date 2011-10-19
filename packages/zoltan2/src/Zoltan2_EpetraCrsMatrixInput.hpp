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

template <typename User>
class EpetraCrsMatrixInput : public XpetraCrsMatrixInput<User>{
private:

    RCP<const User> inmatrix_;
public:

  typedef typename InputAdapter<User>::scalar_t scalar_t;
  typedef typename InputAdapter<User>::lno_t    lno_t;
  typedef typename InputAdapter<User>::gno_t    gno_t;
  typedef typename InputAdapter<User>::lid_t    lid_t;
  typedef typename InputAdapter<User>::gid_t    gid_t;
  typedef typename InputAdapter<User>::node_t   node_t;

  std::string inputAdapterName()const {return std::string("EpetraCrsMatrix");}

  ~EpetraCrsMatrixInput() { }

  /*! Constructor
   */
  EpetraCrsMatrixInput(const RCP<const User> &matrix):
    XpetraCrsMatrixInput<User> (
      rcp(new Xpetra::EpetraCrsMatrix(rcp_const_cast<User>(matrix)))) 
  {
    inmatrix_ = matrix;
  }
};

} //namespace Zoltan2

#endif
