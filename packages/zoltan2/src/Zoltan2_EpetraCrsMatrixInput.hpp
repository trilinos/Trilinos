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
  static inline std::string name() {return "Epetra_CrsMatrix";}
  static inline RCP<const Xpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> >
    convertToXpetra(const RCP<const Epetra_CrsMatrix> &a)
    {
      return rcp(new Xpetra::EpetraCrsMatrix(
                             rcp_const_cast<Epetra_CrsMatrix>(a)));
    }
};

/*! Zoltan2::EpetraCrsMatrixInput
    \brief Provides access for Zoltan2 to Epetra_CrsMatrix data.
*/

template <typename User>
class EpetraCrsMatrixInput : public XpetraCrsMatrixInput<User>{
private:

    //KDDRCP<const User> inmatrix_;
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
    XpetraCrsMatrixInput<User> (matrix)
  {
    //KDDinmatrix_ = matrix;
  }
};

} //namespace Zoltan2

#endif
