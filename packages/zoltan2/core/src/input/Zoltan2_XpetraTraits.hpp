// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_XpetraTraits.hpp
    \brief Traits of Xpetra classes, including migration method.
*/

#ifndef _ZOLTAN2_XPETRATRAITS_HPP_
#define _ZOLTAN2_XPETRATRAITS_HPP_

#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_Standards.hpp>

#include <Xpetra_TpetraCrsGraph.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <Tpetra_Vector.hpp>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////////
// Extra traits needed only for Xpetra matrices and graphs

/*! \brief  Defines the traits required for Tpetra, Eptra and Xpetra objects.
 *
 *    Definitions are provided for:
 *
 * \li Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t>
 * \li Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t>
 * \li Xpetra::CrsMatrix<double, int, int, node_t>
 * \li Tpetra::CrsGraph<lno_t, gno_t, node_t>
 * \li Xpetra::CrsGraph<lno_t, gno_t, node_t>
 * \li Xpetra::CrsGraph<int, int, node_t>
 * \li Tpetra::Vector<scalar_t, lno_t, gno_t, node_t>
 * \li Xpetra::Vector<scalar_t, lno_t, gno_t, node_t>
 * \li Xpetra::Vector<double, int, int, node_t>
 * \li Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t>
 * \li Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t>
 * \li Xpetra::MultiVector<double, int, int, node_t>
 */

template <typename User>
struct XpetraTraits
{
  /*! \brief Convert the object to its Xpetra wrapped version.
   */
  static inline RCP<User> convertToXpetra(const RCP<User> &a)
  {
    return a;
  }

  /*! \brief The objects global ordinal data type.
   */
  typedef default_gno_t gno_t;

  /*! \brief The objects local ordinal data type.
   */
  typedef default_lno_t lno_t;

  /*! \brief Migrate the object
   *  Given a user object and a new row distribution, create and
   *  return a new user object with the new distribution.
   */

  static RCP<User> doMigration(const User &from,
      size_t numLocalRows, const gno_t *myNewRows)
  {
    return Teuchos::null;
  }
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS

//////////////////////////////////////////////////////////////////////////////
// Tpetra::CrsMatrix
template <typename scalar_t,
          typename lno_t,
          typename gno_t,
          typename node_t>
struct XpetraTraits<Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> >
{
  typedef typename Xpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> xmatrix_t;
  typedef typename Xpetra::TpetraCrsMatrix<scalar_t,lno_t,gno_t,node_t> xtmatrix_t;
  typedef typename Tpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> tmatrix_t;

  static inline RCP<xmatrix_t> convertToXpetra(const RCP<tmatrix_t> &a)
  {
    return rcp(new xtmatrix_t(a));
  }

  static RCP<tmatrix_t> doMigration(const tmatrix_t &from,
      size_t numLocalRows, const gno_t *myNewRows)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;

    // source map
    const RCP<const map_t> &smap = from.getRowMap();
    gno_t numGlobalRows = smap->getGlobalNumElements();
    gno_t base = smap->getMinAllGlobalIndex();

    // target map
    ArrayView<const gno_t> rowList(myNewRows, numLocalRows);
    const RCP<const Teuchos::Comm<int> > &comm = from.getComm();
    RCP<const map_t> tmap = rcp(new map_t(numGlobalRows, rowList, base, comm));

    // importer
    Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

    // target matrix
    // Chris Siefert proposed using the following to make migration
    // more efficient.
    // By default, the Domain and Range maps are the same as in "from".
    // As in the original code, we instead set them both to tmap.
    // The assumption is a square matrix.
    // TODO:  what about rectangular matrices?
    // TODO:  Should choice of domain/range maps be an option to this function?

    // KDD 3/7/16:  disabling Chris' new code to avoid dashboard failures;
    // KDD 3/7/16:  can re-enable when issue #114 is fixed.
    // KDD 3/7/16:  when re-enable CSIEFERT code, can comment out
    // KDD 3/7/16:  "Original way" code.
    // KDD 1/27/17: Re-enabling Chris' code, as this issue is resolved.
    RCP<tmatrix_t> M;
    from.importAndFillComplete(M, importer, tmap, tmap);

    //// Original way we did it:
    ////
    //int oldNumElts = smap->getLocalNumElements();
    //int newNumElts = numLocalRows;

    //// number of non zeros in my new rows
    //typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> vector_t;
    //vector_t numOld(smap);  // TODO These vectors should have scalar=size_t,
    //vector_t numNew(tmap);  // but ETI does not yet support that.
    //for (int lid=0; lid < oldNumElts; lid++){
      //numOld.replaceGlobalValue(smap->getGlobalElement(lid),
        //scalar_t(from.getNumEntriesInLocalRow(lid)));
    //}
    //numNew.doImport(numOld, importer, Tpetra::INSERT);

    // TODO Could skip this copy if could declare vector with scalar=size_t.
    //ArrayRCP<size_t> nnz(newNumElts);
    //if (newNumElts > 0){
      //ArrayRCP<scalar_t> ptr = numNew.getDataNonConst(0);
      //for (int lid=0; lid < newNumElts; lid++){
        //nnz[lid] = static_cast<size_t>(ptr[lid]);
      //}
    //}

    //RCP<tmatrix_t> M = rcp(new tmatrix_t(tmap, nnz, Tpetra::StaticProfile));
    //M->doImport(from, importer, Tpetra::INSERT);
    //M->fillComplete();

    // End of original way we did it.
    return M;
  }
};


//////////////////////////////////////////////////////////////////////////////
// Xpetra::CrsMatrix
template <typename scalar_t,
          typename lno_t,
          typename gno_t,
          typename node_t>
struct XpetraTraits<Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> >
{
  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> x_matrix_t;
  typedef Xpetra::TpetraCrsMatrix<scalar_t, lno_t, gno_t, node_t> xt_matrix_t;
  typedef Tpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> t_matrix_t;

  static inline RCP<x_matrix_t> convertToXpetra(const RCP<x_matrix_t > &a)
  {
    return a;
  }

  static RCP<x_matrix_t> doMigration(const x_matrix_t &from,
      size_t numLocalRows, const gno_t *myNewRows)
  {
    {
      // Do the import with the Tpetra::CrsMatrix traits object
      const xt_matrix_t *xtm = dynamic_cast<const xt_matrix_t *>(&from);
      RCP<const t_matrix_t> tm = xtm->getTpetra_CrsMatrix();

      RCP<t_matrix_t> tmnew = XpetraTraits<t_matrix_t>::doMigration(
        *tm, numLocalRows, myNewRows);

      RCP<x_matrix_t> xmnew = XpetraTraits<t_matrix_t>::convertToXpetra(tmnew);

      return xmnew;
    }
  }
};

//////////////////////////////////////////////////////////////////////////////
// Xpetra::CrsMatrix specialization

template <typename node_t>
struct XpetraTraits<Xpetra::CrsMatrix<double, int, int, node_t> >
{
  typedef double scalar_t;
  typedef int lno_t;
  typedef int gno_t;
  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> x_matrix_t;
  typedef Xpetra::TpetraCrsMatrix<scalar_t, lno_t, gno_t, node_t> xt_matrix_t;
  typedef Tpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> t_matrix_t;

  static inline RCP<x_matrix_t> convertToXpetra(const RCP<x_matrix_t > &a)
  {
    return a;
  }

  static RCP<x_matrix_t> doMigration(const x_matrix_t &from,
      size_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from.getRowMap()->lib();

    {
      // Do the import with the Tpetra::CrsMatrix traits object
      const xt_matrix_t *xtm = dynamic_cast<const xt_matrix_t *>(&from);
      RCP<const t_matrix_t> tm = xtm->getTpetra_CrsMatrix();

      RCP<t_matrix_t> tmnew = XpetraTraits<t_matrix_t>::doMigration(
        *tm, numLocalRows, myNewRows);

      RCP<x_matrix_t> xmnew = XpetraTraits<t_matrix_t>::convertToXpetra(tmnew);

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

  static inline RCP<xgraph_t> convertToXpetra(const RCP<tgraph_t> &a)
    {
      return rcp(new xtgraph_t(a));
    }

  static RCP<tgraph_t> doMigration(const tgraph_t &from,
      size_t numLocalRows, const gno_t *myNewRows)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;

    // source map
    const RCP<const map_t> &smap = from.getRowMap();
    int oldNumElts = smap->getLocalNumElements();
    gno_t numGlobalRows = smap->getGlobalNumElements();
    gno_t base = smap->getMinAllGlobalIndex();

    // target map
    ArrayView<const gno_t> rowList(myNewRows, numLocalRows);
    const RCP<const Teuchos::Comm<int> > &comm = from.getComm();
    RCP<const map_t> tmap = rcp(new map_t(numGlobalRows, rowList, base, comm));

    // importer
    Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

    // number of entries in my new rows
    typedef Tpetra::Vector<gno_t, lno_t, gno_t, node_t> vector_t;
    vector_t numOld(smap);
    vector_t numNew(tmap);
    for (int lid=0; lid < oldNumElts; lid++){
      numOld.replaceGlobalValue(smap->getGlobalElement(lid),
        from.getNumEntriesInLocalRow(lid));
    }
    numNew.doImport(numOld, importer, Tpetra::INSERT);

    size_t numElts = tmap->getLocalNumElements();
    ArrayRCP<const gno_t> nnz;
    if (numElts > 0)
      nnz = numNew.getData(0);    // hangs if vector len == 0

    ArrayRCP<const size_t> nnz_size_t;

    if (numElts && sizeof(gno_t) != sizeof(size_t)){
      size_t *vals = new size_t [numElts];
      nnz_size_t = arcp(vals, 0, numElts, true);
      for (size_t i=0; i < numElts; i++){
        vals[i] = static_cast<size_t>(nnz[i]);
      }
    }
    else{
      nnz_size_t = arcp_reinterpret_cast<const size_t>(nnz);
    }

    // target graph
    RCP<tgraph_t> G = rcp(new tgraph_t(tmap, nnz_size_t()));

    G->doImport(from, importer, Tpetra::INSERT);
    G->fillComplete();

    return G;
  }

};


//////////////////////////////////////////////////////////////////////////////
// Xpetra::CrsGraph
template <typename lno_t,
          typename gno_t,
          typename node_t>
struct XpetraTraits<Xpetra::CrsGraph<lno_t, gno_t, node_t> >
{
  typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> x_graph_t;
  typedef Xpetra::TpetraCrsGraph<lno_t, gno_t, node_t> xt_graph_t;
  typedef Tpetra::CrsGraph<lno_t,gno_t,node_t> t_graph_t;

  static inline RCP<x_graph_t> convertToXpetra(const RCP<x_graph_t> &a)
  {
    return a;
  }

  static RCP<x_graph_t> doMigration(const x_graph_t &from,
      size_t numLocalRows, const gno_t *myNewRows)
  {
    {
      // Do the import with the Tpetra::CrsGraph traits object
      const xt_graph_t *xtg = dynamic_cast<const xt_graph_t *>(&from);
      RCP<const t_graph_t> tg = xtg->getTpetra_CrsGraph();

      RCP<t_graph_t> tgnew = XpetraTraits<t_graph_t>::doMigration(
        *tg, numLocalRows, myNewRows);

      RCP<x_graph_t> xgnew = XpetraTraits<t_graph_t>::convertToXpetra(tgnew);
      return xgnew;
    }
  }
};


//////////////////////////////////////////////////////////////////////////////
// Xpetra::CrsGraph specialization
template < typename node_t>
struct XpetraTraits<Xpetra::CrsGraph<int, int, node_t> >
{
  typedef int lno_t;
  typedef int gno_t;
  typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> x_graph_t;
  typedef Xpetra::TpetraCrsGraph<lno_t, gno_t, node_t> xt_graph_t;
  typedef Tpetra::CrsGraph<lno_t,gno_t,node_t> t_graph_t;

  static inline RCP<x_graph_t> convertToXpetra(const RCP<x_graph_t> &a)
  {
    return a;
  }

  static RCP<x_graph_t> doMigration(const x_graph_t &from,
      size_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from.getRowMap()->lib();

    {
      // Do the import with the Tpetra::CrsGraph traits object
      const xt_graph_t *xtg = dynamic_cast<const xt_graph_t *>(&from);
      RCP<const t_graph_t> tg = xtg->getTpetra_CrsGraph();

      RCP<t_graph_t> tgnew = XpetraTraits<t_graph_t>::doMigration(
        *tg, numLocalRows, myNewRows);

      RCP<x_graph_t> xgnew = XpetraTraits<t_graph_t>::convertToXpetra(tgnew);

      return xgnew;
    }
  }
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

  static inline RCP<x_vector_t> convertToXpetra(const RCP<t_vector_t> &a)
  {
    return rcp(new xt_vector_t(a));
  }

  static RCP<t_vector_t> doMigration(const t_vector_t &from,
      size_t numLocalElts, const gno_t *myNewElts)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;

    // source map
    const RCP<const map_t> &smap = from.getMap();
    gno_t numGlobalElts = smap->getGlobalNumElements();
    gno_t base = smap->getMinAllGlobalIndex();

    // target map
    ArrayView<const gno_t> eltList(myNewElts, numLocalElts);
    const RCP<const Teuchos::Comm<int> > comm = from.getMap()->getComm();
    RCP<const map_t> tmap = rcp(new map_t(numGlobalElts, eltList, base, comm));

    // importer
    Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

    // target vector
    RCP<t_vector_t> V =
      Tpetra::createVector<scalar_t,lno_t,gno_t,node_t>(tmap);
    V->doImport(from, importer, Tpetra::INSERT);

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

  static inline RCP<x_vector_t> convertToXpetra(const RCP<x_vector_t> &a)
  {
    return a;
  }

  static RCP<x_vector_t> doMigration(const x_vector_t &from,
      size_t numLocalRows, const gno_t *myNewRows)
  {
    {
      // Do the import with the Tpetra::Vector traits object
      const xt_vector_t *xtv = dynamic_cast<const xt_vector_t *>(&from);
      RCP<const t_vector_t> tv = xtv->getTpetra_Vector();

      RCP<t_vector_t> tvnew = XpetraTraits<t_vector_t>::doMigration(
        *tv, numLocalRows, myNewRows);

      RCP<x_vector_t> xvnew = XpetraTraits<t_vector_t>::convertToXpetra(tvnew);

      return xvnew;
    }
  }
};

//////////////////////////////////////////////////////////////////////////////
// Xpetra::Vector specialization
template <typename node_t>
struct XpetraTraits<Xpetra::Vector<double, int, int, node_t> >
{
  typedef double scalar_t;
  typedef int lno_t;
  typedef int gno_t;
  typedef Xpetra::Vector<scalar_t, lno_t, gno_t, node_t> x_vector_t;
  typedef Xpetra::TpetraVector<scalar_t, lno_t, gno_t, node_t> xt_vector_t;
  typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> t_vector_t;

  static inline RCP<x_vector_t> convertToXpetra(const RCP<x_vector_t> &a)
  {
    return a;
  }

  static RCP<x_vector_t> doMigration(const x_vector_t &from,
      size_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from.getMap()->lib();

    {
      // Do the import with the Tpetra::Vector traits object
      const xt_vector_t *xtv = dynamic_cast<const xt_vector_t *>(&from);
      RCP<t_vector_t> tv = xtv->getTpetra_Vector();

      RCP<t_vector_t> tvnew = XpetraTraits<t_vector_t>::doMigration(
        *tv, numLocalRows, myNewRows);

      RCP<x_vector_t> xvnew = XpetraTraits<t_vector_t>::convertToXpetra(tvnew);

      return xvnew;
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

  static inline RCP<x_vector_t> convertToXpetra(const RCP<t_vector_t> &a)
  {
    return rcp(new xt_vector_t(a));
  }

  static RCP<t_vector_t> doMigration(const t_vector_t &from,
      size_t numLocalElts, const gno_t *myNewElts)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;

    // source map
    const RCP<const map_t> &smap = from.getMap();
    gno_t numGlobalElts = smap->getGlobalNumElements();
    gno_t base = smap->getMinAllGlobalIndex();

    // target map
    ArrayView<const gno_t> eltList(myNewElts, numLocalElts);
    const RCP<const Teuchos::Comm<int> > comm = from.getMap()->getComm();
    RCP<const map_t> tmap = rcp(new map_t(numGlobalElts, eltList, base, comm));

    // importer
    Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

    // target vector
    RCP<t_vector_t> MV = rcp(
      new t_vector_t(tmap, from.getNumVectors(), true));
    MV->doImport(from, importer, Tpetra::INSERT);

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

  static inline RCP<x_mvector_t> convertToXpetra(const RCP<x_mvector_t> &a)
  {
    return a;
  }

  static RCP<x_mvector_t> doMigration(const x_mvector_t &from,
      size_t numLocalRows, const gno_t *myNewRows)
  {
    {
      // Do the import with the Tpetra::MultiVector traits object
      const xt_mvector_t *xtv = dynamic_cast<const xt_mvector_t *>(&from);
      RCP<t_mvector_t> tv = xtv->getTpetra_MultiVector();

      RCP<t_mvector_t> tvnew = XpetraTraits<t_mvector_t>::doMigration(
        *tv, numLocalRows, myNewRows);

      RCP<x_mvector_t> xvnew = XpetraTraits<t_mvector_t>::convertToXpetra(tvnew);

      return xvnew;
    }
  }
};

//////////////////////////////////////////////////////////////////////////////
// Xpetra::MultiVector specialization
template <typename node_t>
struct XpetraTraits<Xpetra::MultiVector<double, int, int, node_t> >
{
  typedef double scalar_t;
  typedef int lno_t;
  typedef int gno_t;
  typedef Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> x_mvector_t;
  typedef Xpetra::TpetraMultiVector<scalar_t, lno_t, gno_t, node_t> xt_mvector_t;
  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> t_mvector_t;

  static inline RCP<x_mvector_t> convertToXpetra(const RCP<x_mvector_t> &a)
  {
    return a;
  }

  static RCP<x_mvector_t> doMigration(const x_mvector_t &from,
      size_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from.getMap()->lib();

    {
      // Do the import with the Tpetra::MultiVector traits object
      const xt_mvector_t *xtv = dynamic_cast<const xt_mvector_t *>(&from);
      RCP<t_mvector_t> tv = xtv->getTpetra_MultiVector();

      RCP<t_mvector_t> tvnew = XpetraTraits<t_mvector_t>::doMigration(
        *tv, numLocalRows, myNewRows);

      RCP<x_mvector_t> xvnew = XpetraTraits<t_mvector_t>::convertToXpetra(tvnew);

      return xvnew;
    }
  }
};

#endif // DOXYGEN_SHOULD_SKIP_THIS

}  //namespace Zoltan2

#endif // _ZOLTAN2_XPETRATRAITS_HPP_
