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

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef LocalOrdinal  lid_t;
  typedef GlobalOrdinal gid_t;
  typedef Node          node_t;
};

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

  typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> crsMatrix;

  /*! Name of input adapter type.
   */
  std::string inputAdapterName()const {return std::string("TpetraCrsMatrix");}

  /*! Destructor
   */
  ~TpetraCrsMatrixInput() { }

  /*! Constructor 
   */
  TpetraCrsMatrixInput(const RCP<const crsMatrix> &matrix):
    XpetraCrsMatrixInput<User>(
      Teuchos::rcp(new Xpetra::TpetraCrsMatrix<scalar_t, lno_t, gno_t, node_t>(
        Teuchos::rcp_const_cast<crsMatrix>(matrix))))
    
  {
    inmatrix_ = matrix;
  }

  /*! Access to matrix that instantiated adapter
   */
  RCP<const crsMatrix> getMatrix() const
  { 
    return inmatrix_;
  }

private:

  RCP<const crsMatrix > inmatrix_;
};
} // namespace

#endif
