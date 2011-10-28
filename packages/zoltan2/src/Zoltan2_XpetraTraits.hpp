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
#include <Xpetra_EpetraUtils.hpp>

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

  typedef typename InputTraits<User>::gno_t gno_t;
  typedef typename InputTraits<User>::lno_t lno_t;

  /*! Given a user object and a new row distribution, create and
   *  return a *  new user object with the new distribution.
   */

  static RCP<User> doImport(const RCP<const User> &from,
      gno_t numGlobalRows, lno_t numLocalRows,
      const gno_t *myNewRows, gno_t base, const Teuchos::Comm<int> &comm)
  {
    return from;
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

  static RCP<tmatrix_t> doImport(const RCP<const tmatrix_t> &from,
      gno_t numGlobalRows, lno_t numLocalRows,
      const gno_t *myNewRows, gno_t base, const Teuchos::Comm<int> &comm)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;

    // source map
    const RCP<const map_t> &smap = from->getRowMap();
    int oldNumElts = smap->getNodeNumElements();

    // target map
    ArrayView<const gno_t> rowList(myNewRows, numLocalRows);
    RCP<const Teuchos::Comm<int> > commPtr(&comm);
    RCP<const map_t> tmap = rcp(
      new map_t(numGlobalRows, rowList, base, commPtr));
    int newNumElts = numLocalRows;

    // importer
    Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

    // number of non zeros in my new rows
    typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> vector_t;
    vector_t numOld(smap);
    vector_t numNew(tmap);
    for (int lid=0; lid < oldNumElts; lid++){
      numOld.replaceGlobalValue(smap->getGlobalElement(lid), 
        double(from->getNumEntriesInLocalRow(lid)));
    }
    numNew.doImport(numOld, importer, Tpetra::INSERT);

    ArrayRCP<size_t> nnz(newNumElts);
    for (int lid=0; lid < newNumElts; lid++){
      nnz[lid] = static_cast<size_t>(*numNew.getDataNonConst(lid));
    }

    // target matrix
    RCP<tmatrix_t> M = rcp(new tmatrix_t(tmap, nnz, Tpetra::StaticProfile));
    M->doImport(*from, importer, Tpetra::INSERT);

    return M;
  }
};

//////////////////////////////////////////////////////////////////////////////
// Epetra_CrsMatrix

template <>
struct XpetraTraits<Epetra_CrsMatrix>
{
  typedef double scalar_t;
  typedef int lno_t;
  typedef int gno_t;
  typedef Zoltan2::default_node_t node_t;

  static inline RCP<const Xpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> >
    convertToXpetra(const RCP<const Epetra_CrsMatrix> &a)
    {
      return rcp(new Xpetra::EpetraCrsMatrix(
                             rcp_const_cast<Epetra_CrsMatrix>(a)));
    }


  static RCP<Epetra_CrsMatrix> doImport(const RCP<const Epetra_CrsMatrix> &from,
      gno_t numGlobalRows, lno_t numLocalRows,
      const gno_t *myNewRows, gno_t base, const Teuchos::Comm<int> &comm)
  {
    RCP<const Comm<int> > commptr = 
      rcp_const_cast<const Comm<int> >(rcp(&comm));
    const RCP<const Epetra_Comm> ecomm = Xpetra::toEpetra(commptr);
    commptr.release();
 
    // source map
    const Epetra_Map &smap = from->RowMap();
    int oldNumElts = smap.NumMyElements();

    // target map
    Epetra_Map tmap(numGlobalRows, numLocalRows, myNewRows, base, *ecomm);
    int newNumElts = numLocalRows;

    // importer
    Epetra_Import importer(tmap, smap);

    // number of non zeros in my new rows
    Epetra_MultiVector numOld(smap, 1);
    Epetra_MultiVector numNew(tmap, 1);
    for (int lid=0; lid < oldNumElts; lid++){
      numOld.ReplaceGlobalValue(
        smap.GID(lid), 1, double(from->NumMyEntries(lid)));
    }
    numNew.Import(numOld, importer, Insert);

    Array<int> nnz(newNumElts);
    for (int lid=0; lid < newNumElts; lid++){
      nnz[lid] = static_cast<int>(*(numNew[lid]));
    }

    // target matrix
    RCP<Epetra_CrsMatrix> M = rcp(
      new Epetra_CrsMatrix(Copy, tmap, nnz.getRawPtr(), true));
    M->Import(*from, importer, Insert);

    return M;
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
  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> x_matrix_t;
  typedef Xpetra::TpetraCrsMatrix<scalar_t, lno_t, gno_t, node_t> xt_matrix_t;
  typedef Xpetra::EpetraCrsMatrix xe_matrix_t;
  typedef Tpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> t_matrix_t; 
  typedef Epetra_CrsMatrix e_matrix_t; 

  static inline RCP<const x_matrix_t>
    convertToXpetra( const RCP<const x_matrix_t > &a)
    {
      return a;
    }

  static RCP<x_matrix_t> doImport(const RCP<const x_matrix_t> &from,
      gno_t numGlobalRows, lno_t numLocalRows,
      const gno_t *myNewRows, gno_t base, const Teuchos::Comm<int> &comm)
  {
    Xpetra::UnderlyingLib l = from.getRowMap()->lib();

    if (l == Xpetra::UseEpetra){
      // Do the import with the Epetra_CrsMatrix traits object
      RCP< const xe_matrix_t> xem = 
       rcp_implicit_cast<const xe_matrix_t>(from);
      RCP<const e_matrix_t> em = xem->getEpetra_CrsMatrix();
      RCP<const e_matrix_t> & emnew = 
        XpetraTraits<e_matrix_t>::doImport(em, numGlobalRows, 
          numLocalRows, myNewRows, base, comm);
      RCP<const xe_matrix_t> xemnew = 
        XpetraTraits<xe_matrix_t>::convertToXpetra(emnew);

      RCP<const x_matrix_t> xmnew = rcp_implicit_cast<const x_matrix_t>(xemnew);

      return xmnew;

    } else{
      // Do the import with the Tpetra::CrsMatrix traits object
      RCP< const xt_matrix_t> xtm = rcp_implicit_cast<const xt_matrix_t>(from);
      RCP<const t_matrix_t> tm = xtm->getTpetra_CrsMatrix();

      RCP<const t_matrix_t> &tmnew = XpetraTraits<t_matrix_t>::doImport(
        tm, numGlobalRows, numLocalRows, myNewRows, base, comm);

      RCP<const xt_matrix_t> xtmnew = 
        XpetraTraits<xt_matrix_t>::convertToXpetra(tmnew);
      RCP<const x_matrix_t> xmnew = rcp_implicit_cast(xtmnew);
      return xmnew;
    }
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
