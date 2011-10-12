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

template < >
struct InputTraits<Epetra_CrsMatrix>
{
  typedef double scalar_t;
  typedef int lno_t;
  typedef int gno_t;
  typedef int lid_t;
  typedef int gid_t;
  typedef Kokkos::DefaultNode::DefaultNodeType node_t;
};

/*! Zoltan2::EpetraCrsMatrixInput
    \brief Provides access for Zoltan2 to Epetra_CrsMatrix data.
*/

template <typename User>
class EpetraCrsMatrixInput : public XpetraCrsMatrixInput<User>{
private:

    RCP<const Epetra_CrsMatrix> inmatrix_;
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
  EpetraCrsMatrixInput(const RCP<const Epetra_CrsMatrix> matrix):
    XpetraCrsMatrixInput<User> (
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
