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
#include <Xpetra_EpetraVector.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <Xpetra_EpetraUtils.hpp>
#include <Tpetra_Vector.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_Standards.hpp>

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

  typedef long gno_t;
  typedef int lno_t;

  //TODO generate a compile time error if we end up here

  /*! Given a user object and a new row distribution, create and
   *  return a new user object with the new distribution.
   */

  static RCP<User> doImport(const RCP<const User> &from,
      lno_t numLocalRows, const gno_t *myNewRows, gno_t base)
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
      lno_t numLocalRows, const gno_t *myNewRows, gno_t base)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;

    // source map
    const RCP<const map_t> &smap = from->getRowMap();
    int oldNumElts = smap->getNodeNumElements();
    gno_t numGlobalRows = smap->getGlobalNumElements();

    // target map
    ArrayView<const gno_t> rowList(myNewRows, numLocalRows);
    const RCP<const Teuchos::Comm<int> > &comm = from->getComm();
    RCP<const map_t> tmap = rcp(
      new map_t(numGlobalRows, rowList, base, comm));
    int newNumElts = numLocalRows;

    // importer
    Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

    // number of non zeros in my new rows
    typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> vector_t;
    vector_t numOld(smap);
    vector_t numNew(tmap);
    for (int lid=0; lid < oldNumElts; lid++){
      numOld.replaceGlobalValue(smap->getGlobalElement(lid), 
        scalar_t(from->getNumEntriesInLocalRow(lid)));
    }
    numNew.doImport(numOld, importer, Tpetra::INSERT);

    ArrayRCP<size_t> nnz(newNumElts);
    ArrayRCP<scalar_t> ptr = numNew.getDataNonConst(0);
    for (int lid=0; lid < newNumElts; lid++){
      nnz[lid] = static_cast<size_t>(ptr[lid]);
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
  typedef InputTraits<Epetra_CrsMatrix>::scalar_t scalar_t;
  typedef InputTraits<Epetra_CrsMatrix>::lno_t lno_t;
  typedef InputTraits<Epetra_CrsMatrix>::gno_t gno_t;
  typedef InputTraits<Epetra_CrsMatrix>::node_t node_t;

  static inline RCP<const Xpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> >
    convertToXpetra(const RCP<const Epetra_CrsMatrix> &a)
    {
      return rcp(new Xpetra::EpetraCrsMatrix(
                             rcp_const_cast<Epetra_CrsMatrix>(a)));
    }


  static RCP<Epetra_CrsMatrix> doImport(
      const RCP<const Epetra_CrsMatrix> &from,
      lno_t numLocalRows, const gno_t *myNewRows, gno_t base)
  {
    // source map
    const Epetra_Map &smap = from->RowMap();
    int oldNumElts = smap.NumMyElements();
    gno_t numGlobalRows = smap.NumGlobalElements();

    // target map
    const Epetra_Comm &comm = from->Comm();
    Epetra_Map tmap(numGlobalRows, numLocalRows, myNewRows, base, comm);
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
      lno_t numLocalRows, const gno_t *myNewRows, gno_t base)
  {
    Xpetra::UnderlyingLib l = from.getRowMap()->lib();

    if (l == Xpetra::UseEpetra){
      // Do the import with the Epetra_CrsMatrix traits object
      RCP< const xe_matrix_t> xem = 
       rcp_implicit_cast<const xe_matrix_t>(from);
      RCP<const e_matrix_t> em = xem->getEpetra_CrsMatrix();
      RCP<const e_matrix_t> & emnew = 
        XpetraTraits<e_matrix_t>::doImport(em,
          numLocalRows, myNewRows, base);
      RCP<const xe_matrix_t> xemnew = 
        XpetraTraits<xe_matrix_t>::convertToXpetra(emnew);

      RCP<const x_matrix_t> xmnew = rcp_implicit_cast<const x_matrix_t>(xemnew);

      return xmnew;

    } else{
      // Do the import with the Tpetra::CrsMatrix traits object
      RCP< const xt_matrix_t> xtm = rcp_implicit_cast<const xt_matrix_t>(from);
      RCP<const t_matrix_t> tm = xtm->getTpetra_CrsMatrix();

      RCP<const t_matrix_t> &tmnew = XpetraTraits<t_matrix_t>::doImport(
        tm, numLocalRows, myNewRows, base);

      RCP<const xt_matrix_t> xtmnew = 
        XpetraTraits<xt_matrix_t>::convertToXpetra(tmnew);
      RCP<const x_matrix_t> xmnew = rcp_implicit_cast(xtmnew);
      return xmnew;
    }
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

  static RCP<tgraph_t> doImport(const RCP<const tgraph_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows, gno_t base)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;

    // source map
    const RCP<const map_t> &smap = from->getRowMap();
    int oldNumElts = smap->getNodeNumElements();
    gno_t numGlobalRows = smap->getGlobalNumElements();

    // target map
    ArrayView<const gno_t> rowList(myNewRows, numLocalRows);
    const RCP<const Teuchos::Comm<int> > &comm = from->getComm();
    RCP<const map_t> tmap = rcp(
      new map_t(numGlobalRows, rowList, base, comm));
    int newNumElts = numLocalRows;

    // importer
    Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

    // number of entries in my new rows
    typedef Tpetra::Vector<size_t, lno_t, gno_t, node_t> vector_t;
    vector_t numOld(smap);
    vector_t numNew(tmap);
    for (int lid=0; lid < oldNumElts; lid++){
      numOld.replaceGlobalValue(smap->getGlobalElement(lid),
        from->getNumEntriesInLocalRow(lid));
    }
    numNew.doImport(numOld, importer, Tpetra::INSERT);

    ArrayRCP<size_t> nnz = numNew.getDataNonConst(0);

    // target graph
    RCP<tgraph_t> G = rcp(new tgraph_t(tmap, nnz, Tpetra::StaticProfile));
    G->doImport(*from, importer, Tpetra::INSERT);

    return G;
  }

};

//////////////////////////////////////////////////////////////////////////////
// Epetra_CrsGraph
template < >
struct XpetraTraits<Epetra_CrsGraph>
{
  typedef InputTraits<Epetra_CrsGraph>::lno_t    lno_t;
  typedef InputTraits<Epetra_CrsGraph>::gno_t    gno_t;
  typedef InputTraits<Epetra_CrsGraph>::node_t   node_t;
  static inline RCP<const Xpetra::CrsGraph<lno_t,gno_t,node_t> >
    convertToXpetra(const RCP<const Epetra_CrsGraph> &a)
    {
      return rcp(new Xpetra::EpetraCrsGraph(
                             rcp_const_cast<Epetra_CrsGraph>(a)));
    }
    // TODO doImport
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
    // TODO doImport
};

//////////////////////////////////////////////////////////////////////////////
// Tpetra::Vector
template <typename scalar_t,
          typename lno_t,
          typename gno_t,
          typename node_t>
struct XpetraTraits<Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> >
{
  typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> t_vector_t;
  typedef Xpetra::TpetraVector<scalar_t, lno_t, gno_t, node_t> xt_vector_t;
  typedef Xpetra::Vector<scalar_t, lno_t, gno_t, node_t> x_vector_t;

  static inline RCP<const x_vector_t>
    convertToXpetra(const RCP<const t_vector_t> &a)
    {
      return rcp(new xt_vector_t(rcp_const_cast<t_vector_t>(a)));
    }

  static RCP<t_vector_t> doImport(const RCP<const t_vector_t> &from,
      lno_t numLocalElts, const gno_t *myNewElts, gno_t base)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;

    // source map
    const RCP<const map_t> &smap = from->getMap();
    gno_t numGlobalElts = smap->getGlobalNumElements();

    // target map
    ArrayView<const gno_t> eltList(myNewElts, numLocalElts);
    const RCP<const Teuchos::Comm<int> > comm = from->getMap()->getComm();
    RCP<const map_t> tmap = rcp(
      new map_t(numGlobalElts, eltList, base, comm));

    // importer
    Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

    // target vector 
    RCP<t_vector_t> V = 
      Tpetra::createVector<scalar_t,lno_t,gno_t,node_t>(tmap);
    V->doImport(*from, importer, Tpetra::INSERT);

    return V;
  }
};

//////////////////////////////////////////////////////////////////////////////
// Epetra_Vector
template < >
struct XpetraTraits<Epetra_Vector>
{
  typedef InputTraits<Epetra_Vector>::lno_t    lno_t;
  typedef InputTraits<Epetra_Vector>::gno_t    gno_t;
  typedef InputTraits<Epetra_Vector>::node_t   node_t;
  typedef InputTraits<Epetra_Vector>::scalar_t   scalar_t;
  
  typedef Xpetra::Vector<scalar_t, lno_t, gno_t, node_t> x_vector_t;

  static inline RCP<const x_vector_t>
    convertToXpetra(const RCP<const Epetra_Vector> &a)
    {
      RCP<Xpetra::EpetraVector> xev = rcp(
        new Xpetra::EpetraVector(rcp_const_cast<Epetra_Vector>(a)));
      return rcp_implicit_cast<x_vector_t>(xev);
    }

  static RCP<Epetra_Vector> doImport(const RCP<const Epetra_Vector> &from,
      lno_t numLocalElts, const gno_t *myNewElts, gno_t base)
  {
    // source map
    const Epetra_BlockMap &smap = from->Map();
    gno_t numGlobalElts = smap.NumGlobalElements();

    // target map
    const Epetra_Comm &comm = from->Comm();
    const Epetra_BlockMap tmap(numGlobalElts, numLocalElts, myNewElts, 
      1, base, comm);

    // importer
    Epetra_Import importer(tmap, smap);

    // target vector 
    RCP<Epetra_Vector> V = rcp(new Epetra_Vector(tmap, true));
    Epetra_CombineMode c = Insert;
    V->Import(*from, importer, c);

    return V;
  }
};

//////////////////////////////////////////////////////////////////////////////
// Xpetra::Vector
template <typename scalar_t,
          typename lno_t,
          typename gno_t,
          typename node_t>
struct XpetraTraits<Xpetra::Vector<scalar_t, lno_t, gno_t, node_t> >
{
  typedef Xpetra::Vector<scalar_t, lno_t, gno_t, node_t> x_vector_t;
  typedef Xpetra::TpetraVector<scalar_t, lno_t, gno_t, node_t> xt_vector_t;
  typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> t_vector_t;
  typedef Xpetra::EpetraVector xe_vector_t;

  static inline RCP<const x_vector_t>
    convertToXpetra(const RCP<const x_vector_t> &a)
    {
      return a;
    }

  static RCP<x_vector_t> doImport(const RCP<const x_vector_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows, gno_t base)
  {
    Xpetra::UnderlyingLib l = from.getMap()->lib();

    if (l == Xpetra::UseEpetra){
      // Do the import with the Epetra_Vector traits object
      RCP< const xe_vector_t> xem =
       rcp_implicit_cast<const xe_vector_t>(from);
      RCP<const Epetra_Vector> em = xem->getEpetra_Vector();
      RCP<const Epetra_Vector> & emnew =
        XpetraTraits<Epetra_Vector>::doImport(em, 
          numLocalRows, myNewRows, base);
      RCP<const xe_vector_t> xemnew =
        XpetraTraits<xe_vector_t>::convertToXpetra(emnew);

      RCP<const x_vector_t> xmnew = 
        rcp_implicit_cast<const x_vector_t>(xemnew);

      return xmnew;

    } else{
      // Do the import with the Tpetra::Vector traits object
      RCP< const xt_vector_t> xtm = 
        rcp_implicit_cast<const xt_vector_t>(from);
      RCP<const t_vector_t> tm = xtm->getTpetra_Vector();

      RCP<const t_vector_t> &tmnew = XpetraTraits<t_vector_t>::doImport(
        tm, numLocalRows, myNewRows, base);

      RCP<const xt_vector_t> xtmnew =
        XpetraTraits<xt_vector_t>::convertToXpetra(tmnew);
      RCP<const x_vector_t> xmnew = rcp_implicit_cast(xtmnew);
      return xmnew;
    }
  }
};

//////////////////////////////////////////////////////////////////////////////
// Tpetra::MultiVector
template <typename scalar_t,
          typename lno_t,
          typename gno_t,
          typename node_t>
struct XpetraTraits<Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> >
{
  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> t_vector_t;
  typedef Xpetra::TpetraMultiVector<scalar_t, lno_t, gno_t, node_t> xt_vector_t;
  typedef Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> x_vector_t;

  static inline RCP<const x_vector_t>
    convertToXpetra(const RCP<const t_vector_t> &a)
    {
      return rcp(new xt_vector_t(rcp_const_cast<t_vector_t>(a)));
    }

  static RCP<t_vector_t> doImport(const RCP<const t_vector_t> &from,
      lno_t numLocalElts, const gno_t *myNewElts, gno_t base)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;

    // source map
    const RCP<const map_t> &smap = from->getMap();
    gno_t numGlobalElts = smap->getGlobalNumElements();

    // target map
    ArrayView<const gno_t> eltList(myNewElts, numLocalElts);
    const RCP<const Teuchos::Comm<int> > comm = from->getMap()->getComm();
    RCP<const map_t> tmap = rcp(
      new map_t(numGlobalElts, eltList, base, comm));

    // importer
    Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

    // target vector 
    RCP<t_vector_t> MV = rcp(
      new t_vector_t(tmap, from->getNumVectors(), true));
    MV->doImport(*from, importer, Tpetra::INSERT);

    return MV;
  }
};

//////////////////////////////////////////////////////////////////////////////
// Epetra_MultiVector
template < >
struct XpetraTraits<Epetra_MultiVector>
{
  typedef InputTraits<Epetra_MultiVector>::lno_t    lno_t;
  typedef InputTraits<Epetra_MultiVector>::gno_t    gno_t;
  typedef InputTraits<Epetra_MultiVector>::node_t   node_t;
  typedef InputTraits<Epetra_MultiVector>::scalar_t   scalar_t;
  typedef Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> x_mvector_t;

  static inline RCP<const x_mvector_t>
    convertToXpetra(const RCP<const Epetra_MultiVector> &a)
    {
      RCP<Xpetra::EpetraMultiVector> xemv = rcp(
        new Xpetra::EpetraMultiVector(
          rcp_const_cast<Epetra_MultiVector>(a)));
      return rcp_implicit_cast<x_mvector_t>(xemv);
    }

  static RCP<Epetra_MultiVector> doImport(
    const RCP<const Epetra_MultiVector> &from,
    lno_t numLocalElts, const gno_t *myNewElts, gno_t base)
  {
    // source map
    const Epetra_BlockMap &smap = from->Map();
    gno_t numGlobalElts = smap.NumGlobalElements();

    // target map
    const Epetra_Comm &comm = from->Comm();
    const Epetra_BlockMap tmap(numGlobalElts, numLocalElts, myNewElts, 
      1, base, comm);

    // importer
    Epetra_Import importer(tmap, smap);

    // target vector 
    RCP<Epetra_MultiVector> MV = rcp(
      new Epetra_MultiVector(tmap, from->NumVectors(), true));
    Epetra_CombineMode c = Insert;
    MV->Import(*from, importer, c);

    return MV;
  }
};

//////////////////////////////////////////////////////////////////////////////
// Xpetra::MultiVector
template <typename scalar_t,
          typename lno_t,
          typename gno_t,
          typename node_t>
struct XpetraTraits<Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> >
{
  typedef Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> x_mvector_t;
  typedef Xpetra::TpetraMultiVector<scalar_t, lno_t, gno_t, node_t> xt_mvector_t;
  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> t_mvector_t;
  typedef Xpetra::EpetraMultiVector xe_mvector_t;

  static inline RCP<const x_mvector_t>
    convertToXpetra(const RCP<const x_mvector_t> &a)
    {
      return a;
    }

  static RCP<x_mvector_t> doImport(const RCP<const x_mvector_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows, gno_t base)
  {
    Xpetra::UnderlyingLib l = from.getMap()->lib();

    if (l == Xpetra::UseEpetra){
      // Do the import with the Epetra_MultiVector traits object
      RCP< const xe_mvector_t> xem =
       rcp_implicit_cast<const xe_mvector_t>(from);
      RCP<const Epetra_MultiVector> em = xem->getEpetra_MultiVector();
      RCP<const Epetra_MultiVector> & emnew =
        XpetraTraits<Epetra_MultiVector>::doImport(em,
          numLocalRows, myNewRows, base);
      RCP<const xe_mvector_t> xemnew =
        XpetraTraits<xe_mvector_t>::convertToXpetra(emnew);

      RCP<const x_mvector_t> xmnew = 
        rcp_implicit_cast<const x_mvector_t>(xemnew);

      return xmnew;

    } else{
      // Do the import with the Tpetra::MultiVector traits object
      RCP< const xt_mvector_t> xtm = 
        rcp_implicit_cast<const xt_mvector_t>(from);
      RCP<const t_mvector_t> tm = xtm->getTpetra_MultiVector();

      RCP<const t_mvector_t> &tmnew = XpetraTraits<t_mvector_t>::doImport(
        tm, numLocalRows, myNewRows, base);

      RCP<const xt_mvector_t> xtmnew =
        XpetraTraits<xt_mvector_t>::convertToXpetra(tmnew);
      RCP<const x_mvector_t> xmnew = rcp_implicit_cast(xtmnew);
      return xmnew;
    }
  }
};
}  //namespace Zoltan2

#endif // _ZOLTAN2_XPETRATRAITS_HPP_
