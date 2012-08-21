// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_XpetraTraits.hpp
    \brief Traits of Xpetra classes, including migration method.
*/

#ifndef _ZOLTAN2_XPETRATRAITS_HPP_
#define _ZOLTAN2_XPETRATRAITS_HPP_

#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_Standards.hpp>

#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#include <Xpetra_EpetraVector.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <Xpetra_EpetraUtils.hpp>
#include <Tpetra_Vector.hpp>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////////
// Extra traits needed only for Xpetra matrices and graphs

/*! \brief  Defines the traits required for Tpetra, Eptra and Xpetra objects.
 *
 *    Definitions are provided for:
 *
 * \li Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> 
 * \li Epetra_CrsMatrix
 * \li Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> 
 * \li Xpetra::CrsMatrix<double, int, int, node_t> 
 * \li Tpetra::CrsGraph<lno_t, gno_t, node_t> 
 * \li Epetra_CrsGraph
 * \li Xpetra::CrsGraph<lno_t, gno_t, node_t> 
 * \li Xpetra::CrsGraph<int, int, node_t> 
 * \li Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> 
 * \li Epetra_Vector
 * \li Xpetra::Vector<scalar_t, lno_t, gno_t, node_t> 
 * \li Xpetra::Vector<double, int, int, node_t> 
 * \li Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> 
 * \li Epetra_MultiVector
 * \li Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> 
 * \li Xpetra::MultiVector<double, int, int, node_t> 
 */

template <typename User>
struct XpetraTraits
{
  /*! \brief Convert the object to its Xpetra wrapped version.
   */
  static inline RCP<const User> convertToXpetra(const RCP<const User> &a) 
  {
    return a;
  }

  /*! \brief The objects global ordinal data type.
   */
  typedef long gno_t;

  /*! \brief The objects local ordinal data type.
   */
  typedef int lno_t;

  /*! \brief Migrate the object
   *  Given a user object and a new row distribution, create and
   *  return a new user object with the new distribution.
   */

  static RCP<const User> doMigration(const RCP<const User> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    return from;
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

  static RCP<const tmatrix_t> doMigration(const RCP<const tmatrix_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;
    lno_t base = 0;

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
    vector_t numOld(smap);  // TODO These vectors should have scalar = size_t, 
    vector_t numNew(tmap);  // but explicit instantiation does not yet support that.
    for (int lid=0; lid < oldNumElts; lid++){
      numOld.replaceGlobalValue(smap->getGlobalElement(lid), 
        scalar_t(from->getNumEntriesInLocalRow(lid)));
    }
    numNew.doImport(numOld, importer, Tpetra::INSERT);

    // TODO Could skip this copy if could declare vector with scalar = size_t.
    ArrayRCP<size_t> nnz(newNumElts);
    if (newNumElts > 0){
      ArrayRCP<scalar_t> ptr = numNew.getDataNonConst(0);
      for (int lid=0; lid < newNumElts; lid++){
        nnz[lid] = static_cast<size_t>(ptr[lid]);
      }
    }

    // target matrix
    RCP<tmatrix_t> M = rcp(new tmatrix_t(tmap, nnz, Tpetra::StaticProfile));
    M->doImport(*from, importer, Tpetra::INSERT);
    M->fillComplete();

    return rcp_const_cast<const tmatrix_t>(M);
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


  static RCP<Epetra_CrsMatrix> doMigration(
      const RCP<const Epetra_CrsMatrix> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    lno_t base = 0;

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
    Epetra_Vector numOld(smap);
    Epetra_Vector numNew(tmap);

    for (int lid=0; lid < oldNumElts; lid++){
      numOld[lid] = from->NumMyEntries(lid);
    }
    numNew.Import(numOld, importer, Insert);

    Array<int> nnz(newNumElts);
    for (int lid=0; lid < newNumElts; lid++){
      nnz[lid] = static_cast<int>(numNew[lid]);
    }

    // target matrix
    RCP<Epetra_CrsMatrix> M = rcp(
      new Epetra_CrsMatrix(Copy, tmap, nnz.getRawPtr(), true));
    M->Import(*from, importer, Insert);
    M->FillComplete();

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
  typedef Tpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> t_matrix_t; 

  static inline RCP<const x_matrix_t>
    convertToXpetra( const RCP<const x_matrix_t > &a)
    {
      return a;
    }

  static RCP<const x_matrix_t> doMigration(const RCP<const x_matrix_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from->getRowMap()->lib();

    if (lib == Xpetra::UseEpetra){
       throw logic_error("compiler should have used specialization");
    } else{
      // Do the import with the Tpetra::CrsMatrix traits object
      const x_matrix_t *xm = from.get();
      const xt_matrix_t *xtm = dynamic_cast<const xt_matrix_t *>(xm);
      RCP<const t_matrix_t> tm = xtm->getTpetra_CrsMatrix();

      RCP<const t_matrix_t> tmnew = XpetraTraits<t_matrix_t>::doMigration(
        tm, numLocalRows, myNewRows);

      RCP<const x_matrix_t> xmnew = 
        XpetraTraits<t_matrix_t>::convertToXpetra(tmnew);

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
  typedef Xpetra::EpetraCrsMatrix xe_matrix_t;
  typedef Epetra_CrsMatrix e_matrix_t; 

  static inline RCP<const x_matrix_t>
    convertToXpetra( const RCP<const x_matrix_t > &a)
    {
      return a;
    }

  static RCP<const x_matrix_t> doMigration(const RCP<const x_matrix_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from->getRowMap()->lib();
    const x_matrix_t *xm = from.get();

    if (lib == Xpetra::UseEpetra){
      // Do the import with the Epetra_CrsMatrix traits object
      const xe_matrix_t *xem = dynamic_cast<const xe_matrix_t *>(xm);
      RCP<const e_matrix_t> em = xem->getEpetra_CrsMatrix();

      RCP<const e_matrix_t> emnew = XpetraTraits<e_matrix_t>::doMigration(
        em, numLocalRows, myNewRows);

      RCP<const x_matrix_t> xmnew = 
        XpetraTraits<e_matrix_t>::convertToXpetra(emnew);

      return xmnew;
    } else{
      // Do the import with the Tpetra::CrsMatrix traits object
      const xt_matrix_t *xtm = dynamic_cast<const xt_matrix_t *>(xm);
      RCP<const t_matrix_t> tm = xtm->getTpetra_CrsMatrix();

      RCP<const t_matrix_t> tmnew = XpetraTraits<t_matrix_t>::doMigration(
        tm, numLocalRows, myNewRows);

      RCP<const x_matrix_t> xmnew = 
        XpetraTraits<t_matrix_t>::convertToXpetra(tmnew);

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

  static RCP<const tgraph_t> doMigration(const RCP<const tgraph_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;
    lno_t base = 0;

    // source map
    const RCP<const map_t> &smap = from->getRowMap();
    int oldNumElts = smap->getNodeNumElements();
    gno_t numGlobalRows = smap->getGlobalNumElements();

    // target map
    ArrayView<const gno_t> rowList(myNewRows, numLocalRows);
    const RCP<const Teuchos::Comm<int> > &comm = from->getComm();
    RCP<const map_t> tmap = rcp(
      new map_t(numGlobalRows, rowList, base, comm));

    // importer
    Tpetra::Import<lno_t, gno_t, node_t> importer(smap, tmap);

    // number of entries in my new rows
    typedef Tpetra::Vector<gno_t, lno_t, gno_t, node_t> vector_t;
    vector_t numOld(smap);
    vector_t numNew(tmap);
    for (int lid=0; lid < oldNumElts; lid++){
      numOld.replaceGlobalValue(smap->getGlobalElement(lid),
        from->getNumEntriesInLocalRow(lid));
    }
    numNew.doImport(numOld, importer, Tpetra::INSERT);

    size_t numElts = tmap->getNodeNumElements();
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
    RCP<tgraph_t> G = rcp(new tgraph_t(tmap, nnz_size_t, Tpetra::StaticProfile));

    G->doImport(*from, importer, Tpetra::INSERT);
    G->fillComplete();

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

  static RCP<const Epetra_CrsGraph> doMigration(
      const RCP<const Epetra_CrsGraph> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    lno_t base = 0;

    // source map
    const Epetra_BlockMap &smap = from->RowMap();
    gno_t numGlobalRows = smap.NumGlobalElements();
    lno_t oldNumElts = smap.NumMyElements();

    // target map
    const Epetra_Comm &comm = from->Comm();
    Epetra_BlockMap tmap(numGlobalRows, numLocalRows, 
       myNewRows, 1, base, comm);
    lno_t newNumElts = tmap.NumMyElements();

    // importer
    Epetra_Import importer(tmap, smap);

    // number of non zeros in my new rows
    Epetra_Vector numOld(smap);
    Epetra_Vector numNew(tmap);

    for (int lid=0; lid < oldNumElts; lid++){
      numOld[lid] = from->NumMyIndices(lid);
    }
    numNew.Import(numOld, importer, Insert);

    Array<int> nnz(newNumElts);
    for (int lid=0; lid < newNumElts; lid++){
      nnz[lid] = static_cast<int>(numNew[lid]);
    }

    // target graph
    RCP<Epetra_CrsGraph> G = rcp(
      new Epetra_CrsGraph(Copy, tmap, nnz.getRawPtr(), true));
    G->Import(*from, importer, Insert);
    G->FillComplete();

    return G;
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
  typedef Xpetra::CrsGraph<lno_t, gno_t, node_t> x_graph_t;
  typedef Xpetra::TpetraCrsGraph<lno_t, gno_t, node_t> xt_graph_t;
  typedef Tpetra::CrsGraph<lno_t,gno_t,node_t> t_graph_t; 

  static inline RCP<const x_graph_t>
    convertToXpetra(const RCP<const x_graph_t> &a)
    {
      return a;
    }

  static RCP<const x_graph_t> doMigration(const RCP<const x_graph_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from->getRowMap()->lib();

    if (lib == Xpetra::UseEpetra){
       throw logic_error("compiler should have used specialization");
    } else{
      // Do the import with the Tpetra::CrsGraph traits object
      const x_graph_t *xg = from.get();
      const xt_graph_t *xtg = dynamic_cast<const xt_graph_t *>(xg);
      RCP<const t_graph_t> tg = xtg->getTpetra_CrsGraph();

      RCP<const t_graph_t> tgnew = XpetraTraits<t_graph_t>::doMigration(
        tg, numLocalRows, myNewRows);

      RCP<const x_graph_t> xgnew =
        XpetraTraits<t_graph_t>::convertToXpetra(tgnew);
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
  typedef Xpetra::EpetraCrsGraph xe_graph_t;
  typedef Epetra_CrsGraph e_graph_t; 

  static inline RCP<const x_graph_t>
    convertToXpetra(const RCP<const x_graph_t> &a)
    {
      return a;
    }

  static RCP<const x_graph_t> doMigration(const RCP<const x_graph_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from->getRowMap()->lib();
    const x_graph_t *xg = from.get();

    if (lib == Xpetra::UseEpetra){
      // Do the import with the Epetra_CrsGraph traits object
      const xe_graph_t *xeg = dynamic_cast<const xe_graph_t *>(xg);
      RCP<const e_graph_t> eg = xeg->getEpetra_CrsGraph();

      RCP<const e_graph_t> egnew = XpetraTraits<e_graph_t>::doMigration(
        eg, numLocalRows, myNewRows);

      RCP<const x_graph_t> xgnew =
        XpetraTraits<e_graph_t>::convertToXpetra(egnew);

      return xgnew;
    } else{
      // Do the import with the Tpetra::CrsGraph traits object
      const xt_graph_t *xtg = dynamic_cast<const xt_graph_t *>(xg);
      RCP<const t_graph_t> tg = xtg->getTpetra_CrsGraph();

      RCP<const t_graph_t> tgnew = XpetraTraits<t_graph_t>::doMigration(
        tg, numLocalRows, myNewRows);

      RCP<const x_graph_t> xgnew =
        XpetraTraits<t_graph_t>::convertToXpetra(tgnew);

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

  static inline RCP<const x_vector_t>
    convertToXpetra(const RCP<const t_vector_t> &a)
    {
      return rcp(new xt_vector_t(rcp_const_cast<t_vector_t>(a)));
    }

  static RCP<const t_vector_t> doMigration(const RCP<const t_vector_t> &from,
      lno_t numLocalElts, const gno_t *myNewElts)
  {
    lno_t base = 0;
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

  static RCP<Epetra_Vector> doMigration(const RCP<const Epetra_Vector> &from,
      lno_t numLocalElts, const gno_t *myNewElts)
  {
    lno_t base = 0;
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

  static inline RCP<const x_vector_t>
    convertToXpetra(const RCP<const x_vector_t> &a)
    {
      return a;
    }

  static RCP<const x_vector_t> doMigration(const RCP<const x_vector_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from->getMap()->lib();

    if (lib == Xpetra::UseEpetra){
       throw logic_error("compiler should have used specialization");
    } else{
      // Do the import with the Tpetra::Vector traits object
      const x_vector_t *xv = from.get();
      const xt_vector_t *xtv = dynamic_cast<const xt_vector_t *>(xv);
      RCP<const t_vector_t> tv = xtv->getTpetra_Vector();

      RCP<const t_vector_t> tvnew = XpetraTraits<t_vector_t>::doMigration(
        tv, numLocalRows, myNewRows);

      RCP<const x_vector_t> xvnew =
        XpetraTraits<t_vector_t>::convertToXpetra(tvnew);

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
  typedef Xpetra::EpetraVector xe_vector_t;
  typedef Epetra_Vector e_vector_t;

  static inline RCP<const x_vector_t>
    convertToXpetra(const RCP<const x_vector_t> &a)
    {
      return a;
    }

  static RCP<const x_vector_t> doMigration(const RCP<const x_vector_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from->getMap()->lib();
    const x_vector_t *vec = from.get();

    if (lib == Xpetra::UseEpetra){
      // Do the import with the Epetra_Vector traits object
      const xe_vector_t *xev = dynamic_cast<const xe_vector_t *>(vec);
      RCP<const e_vector_t> ev = rcp(xev->getEpetra_Vector());

      RCP<const e_vector_t> evnew = XpetraTraits<e_vector_t>::doMigration(
        ev, numLocalRows, myNewRows);

      RCP<const x_vector_t> xvnew =
        XpetraTraits<e_vector_t>::convertToXpetra(evnew);
          
      return xvnew;
    } else{
      // Do the import with the Tpetra::Vector traits object
      const xt_vector_t *xtv = dynamic_cast<const xt_vector_t *>(vec);
      RCP<t_vector_t> tv = xtv->getTpetra_Vector();
      RCP<const t_vector_t> ctv = rcp_const_cast<const t_vector_t>(tv);

      RCP<const t_vector_t> tvnew = XpetraTraits<t_vector_t>::doMigration(
        ctv, numLocalRows, myNewRows);

      RCP<const x_vector_t> xvnew =
        XpetraTraits<t_vector_t>::convertToXpetra(tvnew);

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

  static inline RCP<const x_vector_t>
    convertToXpetra(const RCP<const t_vector_t> &a)
    {
      return rcp(new xt_vector_t(rcp_const_cast<t_vector_t>(a)));
    }

  static RCP<const t_vector_t> doMigration(const RCP<const t_vector_t> &from,
      lno_t numLocalElts, const gno_t *myNewElts)
  {
    typedef Tpetra::Map<lno_t, gno_t, node_t> map_t;
    lno_t base = 0;

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

  static RCP<Epetra_MultiVector> doMigration(
    const RCP<const Epetra_MultiVector> &from,
    lno_t numLocalElts, const gno_t *myNewElts)
  {
    lno_t base = 0;
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

  static inline RCP<const x_mvector_t>
    convertToXpetra(const RCP<const x_mvector_t> &a)
    {
      return a;
    }

  static RCP<const x_mvector_t> doMigration(const RCP<const x_mvector_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from->getMap()->lib();

    if (lib == Xpetra::UseEpetra){
       throw logic_error("compiler should have used specialization");
    } else{
      // Do the import with the Tpetra::MultiVector traits object
      const x_mvector_t *xmv = from.get();
      const xt_mvector_t *xtv = dynamic_cast<const xt_mvector_t *>(xmv);
      RCP<t_mvector_t> tv = xtv->getTpetra_MultiVector();
      RCP<const t_mvector_t> ctv = rcp_const_cast<const t_mvector_t>(tv);

      RCP<const t_mvector_t> tvnew = XpetraTraits<t_mvector_t>::doMigration(
        ctv, numLocalRows, myNewRows);

      RCP<const x_mvector_t> xvnew =
        XpetraTraits<t_mvector_t>::convertToXpetra(tvnew);

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
  typedef Xpetra::EpetraMultiVector xe_mvector_t;
  typedef Epetra_MultiVector e_mvector_t;

  static inline RCP<const x_mvector_t>
    convertToXpetra(const RCP<const x_mvector_t> &a)
    {
      return a;
    }

  static RCP<const x_mvector_t> doMigration(const RCP<const x_mvector_t> &from,
      lno_t numLocalRows, const gno_t *myNewRows)
  {
    Xpetra::UnderlyingLib lib = from->getMap()->lib();
    const x_mvector_t *xmv = from.get();

    if (lib == Xpetra::UseEpetra){
      // Do the import with the Epetra_MultiVector traits object
      const xe_mvector_t *xev = dynamic_cast<const xe_mvector_t *>(xmv);
      RCP<e_mvector_t> ev = xev->getEpetra_MultiVector();
      RCP<const e_mvector_t> cev = rcp_const_cast<const e_mvector_t>(ev);

      RCP<const e_mvector_t> evnew = XpetraTraits<e_mvector_t>::doMigration(
        cev, numLocalRows, myNewRows);

      RCP<const x_mvector_t> xvnew =
        XpetraTraits<e_mvector_t>::convertToXpetra(evnew);

      return xvnew;

    } else{
      // Do the import with the Tpetra::MultiVector traits object
      const xt_mvector_t *xtv = dynamic_cast<const xt_mvector_t *>(xmv);
      RCP<t_mvector_t> tv = xtv->getTpetra_MultiVector();
      RCP<const t_mvector_t> ctv = rcp_const_cast<const t_mvector_t>(tv);

      RCP<const t_mvector_t> tvnew = XpetraTraits<t_mvector_t>::doMigration(
        ctv, numLocalRows, myNewRows);

      RCP<const x_mvector_t> xvnew =
        XpetraTraits<t_mvector_t>::convertToXpetra(tvnew);

      return xvnew;
    }
  }
};

#endif // DOXYGEN_SHOULD_SKIP_THIS

}  //namespace Zoltan2

#endif // _ZOLTAN2_XPETRATRAITS_HPP_
