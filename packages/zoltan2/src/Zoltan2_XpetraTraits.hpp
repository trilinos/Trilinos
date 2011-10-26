// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsMatrixInput.hpp

    \brief An input adapter for a Xpetra::CrsMatrix.
*/

#ifndef _ZOLTAN2_XPETRATRAITS_HPP_
#define _ZOLTAN2_XPETRATRAITS_HPP_

#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////////
// Extra traits needed only for Xpetra matrices and graphs
template <typename User>
struct XpetraTraits
{
  static inline RCP<const User> convertToXpetra(const RCP<const User> &a) 
  {
    return a;  // Default implementation
  }
};

//////////////////////////////////////////////////////////////////////////////
// Xpetra::CrsMatrix
// KDDKDD:  Do we need specializations for Xpetra::EpetraCrsMatrix and
// KDDKDD:  Xpetra::TpetraCrsMatrix
template <typename scalar_t,
          typename lno_t,
          typename gno_t,
          typename node_t>
struct XpetraTraits<Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> >
{
  static inline RCP<const Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> >
    convertToXpetra(
      const RCP<const Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> > &a)
    {
      return a;
    }
};

//////////////////////////////////////////////////////////////////////////////
// Tpetra::CrsMatrix
template <typename scalar_t,
          typename lno_t,
          typename gno_t,
          typename node_t>
struct XpetraTraits<Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> >
{
  typedef typename Xpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> 
                   xmatrix_t;
  typedef typename Xpetra::TpetraCrsMatrix<scalar_t,lno_t,gno_t,node_t> 
                   xtmatrix_t;
  typedef typename Tpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> 
                   tmatrix_t;

  static inline RCP<const xmatrix_t> convertToXpetra(
    const RCP<const tmatrix_t> &a)
    {
      return rcp(new xtmatrix_t(rcp_const_cast<tmatrix_t>(a)));
    }
};

//////////////////////////////////////////////////////////////////////////////
// Epetra_CrsMatrix
template < >
struct XpetraTraits<Epetra_CrsMatrix>
{
  typedef InputTraits<Epetra_CrsMatrix>::scalar_t scalar_t;
  typedef InputTraits<Epetra_CrsMatrix>::lno_t    lno_t;
  typedef InputTraits<Epetra_CrsMatrix>::gno_t    gno_t;
  typedef InputTraits<Epetra_CrsMatrix>::node_t   node_t;
  static inline RCP<const Xpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> >
    convertToXpetra(const RCP<const Epetra_CrsMatrix> &a)
    {
      return rcp(new Xpetra::EpetraCrsMatrix(
                             rcp_const_cast<Epetra_CrsMatrix>(a)));
    }
};

//////////////////////////////////////////////////////////////////////////////
// Xpetra::CrsGraph
// KDDKDD Do we need specializations for Xpetra::TpetraCrsGraph and
// KDDKDD Xpetra::EpetraCrsGraph?
template <typename lno_t,
          typename gno_t,
          typename node_t>
struct XpetraTraits<Xpetra::CrsGraph<lno_t, gno_t, node_t> >
{
  static inline RCP<const Xpetra::CrsGraph<lno_t, gno_t, node_t> >
    convertToXpetra(
      const RCP<const Xpetra::CrsGraph<lno_t, gno_t, node_t> > &a)
    {
      return a;
    }
};

//////////////////////////////////////////////////////////////////////////////
// Tpetra::CrsGraph
template <typename lno_t,
          typename gno_t,
          typename node_t>
struct XpetraTraits<Tpetra::CrsGraph<lno_t, gno_t, node_t> >
{
  typedef typename Xpetra::CrsGraph<lno_t, gno_t, node_t> xgraph_t;
  typedef typename Xpetra::TpetraCrsGraph<lno_t, gno_t, node_t> xtgraph_t;
  typedef typename Tpetra::CrsGraph<lno_t, gno_t, node_t> tgraph_t;

  static inline RCP<const xgraph_t> convertToXpetra(
    const RCP<const tgraph_t> &a)
    {
      return rcp(new xtgraph_t(rcp_const_cast<tgraph_t>(a)));
    }
};

//////////////////////////////////////////////////////////////////////////////
// Epetra_CrsGraph
template < >
struct XpetraTraits<Epetra_CrsGraph>
{
  typedef InputTraits<Epetra_CrsMatrix>::lno_t    lno_t;
  typedef InputTraits<Epetra_CrsMatrix>::gno_t    gno_t;
  typedef InputTraits<Epetra_CrsMatrix>::node_t   node_t;
  static inline RCP<const Xpetra::CrsGraph<lno_t,gno_t,node_t> >
    convertToXpetra(const RCP<const Epetra_CrsGraph> &a)
    {
      return rcp(new Xpetra::EpetraCrsGraph(
                             rcp_const_cast<Epetra_CrsGraph>(a)));
    }
};

}  //namespace Zoltan2
  
#endif
