// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_ModelHelpers.hpp
    \brief Defines helper functions for use in the models
*/

#include <Zoltan2_MeshAdapter.hpp>

#ifndef _ZOLTAN2_MODELHELPERS_HPP_
#define _ZOLTAN2_MODELHELPERS_HPP_

namespace Zoltan2 {

//GFD this declaration is really messy is there a better way? I couldn't typedef outside since
//    there is no user type until the function.
template <typename User>
RCP<Tpetra::CrsMatrix<int,
                      typename MeshAdapter<User>::lno_t,
                      typename MeshAdapter<User>::gno_t,
                      typename MeshAdapter<User>::node_t> >
get2ndAdjsMatFromAdjs(const Teuchos::RCP<const MeshAdapter<User> > &ia,
                      const RCP<const Comm<int> > comm,
                      Zoltan2::MeshEntityType sourcetarget,
                      Zoltan2::MeshEntityType through) {
  typedef typename MeshAdapter<User>::gno_t gno_t;
  typedef typename MeshAdapter<User>::lno_t lno_t;
  typedef typename MeshAdapter<User>::node_t node_t;
  typedef typename MeshAdapter<User>::offset_t offset_t;

  typedef int nonzero_t;  // adjacency matrix doesn't need scalar_t
  typedef Tpetra::CrsMatrix<nonzero_t,lno_t,gno_t,node_t>   sparse_matrix_type;
  typedef Tpetra::Map<lno_t, gno_t, node_t>                 map_type;
  typedef Tpetra::global_size_t GST;
  const GST dummy = Teuchos::OrdinalTraits<GST>::invalid ();

/* Find the adjacency for a nodal based decomposition */
  if (ia->availAdjs(sourcetarget, through)) {
    using Teuchos::Array;
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Get node-element connectivity

    const offset_t *offsets=NULL;
    const gno_t *adjacencyIds=NULL;
    ia->getAdjsView(sourcetarget, through, offsets, adjacencyIds);

    const gno_t *Ids=NULL;
    ia->getIDsViewOf(sourcetarget, Ids);

    const gno_t *throughIds=NULL;
    ia->getIDsViewOf(through, throughIds);

    size_t LocalNumIDs = ia->getLocalNumOf(sourcetarget);

    /***********************************************************************/
    /************************* BUILD MAPS FOR ADJS *************************/
    /***********************************************************************/

    RCP<const map_type> sourcetargetMapG;
    RCP<const map_type> throughMapG;

    // count owned nodes
    size_t LocalNumOfThrough = ia->getLocalNumOf(through);

    // Build a list of the global sourcetarget ids...
    gno_t min[2];
    size_t maxcols = 0;
    min[0] = std::numeric_limits<gno_t>::max();
    for (size_t i = 0; i < LocalNumIDs; ++i) {
      if (Ids[i] < min[0]) {
        min[0] = Ids[i];
      }
      size_t ncols = offsets[i+1] - offsets[i];
      if (ncols > maxcols) maxcols = ncols;
    }

    // min(throughIds[i])
    min[1] = std::numeric_limits<gno_t>::max();
    for (size_t i = 0; i < LocalNumOfThrough; ++i) {
      if (throughIds[i] < min[1]) {
        min[1] = throughIds[i];
      }
    }

    gno_t gmin[2];
    Teuchos::reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_MIN, 2, min, gmin);

    //Generate Map for sourcetarget.
    ArrayView<const gno_t> sourceTargetGIDs(Ids, LocalNumIDs);
    sourcetargetMapG = rcp(new map_type(dummy,
                                        sourceTargetGIDs, gmin[0], comm));

    //Create a new map with IDs uniquely assigned to ranks (oneToOneSTMap)
    /*RCP<const map_type> oneToOneSTMap =
      Tpetra::createOneToOne<lno_t, gno_t, node_t>(sourcetargetMapG);*/

    //Generate Map for through.
// TODO
// TODO Could check for max through id as well, and if all through ids are
// TODO in gmin to gmax, then default constructors works below.
// TODO Otherwise, may need a constructor that is not one-to-one containing
// TODO all through entities on processor, followed by call to createOneToOne
// TODO

    ArrayView<const gno_t> throughGIDs(throughIds, LocalNumOfThrough);
    throughMapG = rcp (new map_type(dummy,
                                    throughGIDs, gmin[1], comm));

    //Create a new map with IDs uniquely assigned to ranks (oneToOneTMap)
    RCP<const map_type> oneToOneTMap =
      Tpetra::createOneToOne<lno_t, gno_t, node_t>(throughMapG);

    /***********************************************************************/
    /************************* BUILD GRAPH FOR ADJS ************************/
    /***********************************************************************/

    RCP<sparse_matrix_type> adjsMatrix;

    // Construct Tpetra::CrsGraph objects.
    Array<size_t> rowlens(LocalNumIDs);
    for (size_t localElement = 0; localElement < LocalNumIDs; localElement++)
      rowlens[localElement] = offsets[localElement+1] - offsets[localElement];
    adjsMatrix = rcp (new sparse_matrix_type (sourcetargetMapG,//oneToOneSTMap,
                                              rowlens()));

    Array<nonzero_t> justOneA(maxcols, 1);
    ArrayView<const gno_t> adjacencyIdsAV(adjacencyIds, offsets[LocalNumIDs]);

    for (size_t localElement=0; localElement<LocalNumIDs; ++localElement){
      // Insert all columns for global row Ids[localElement]
      size_t ncols = offsets[localElement+1] - offsets[localElement];
      adjsMatrix->insertGlobalValues(Ids[localElement],
                              adjacencyIdsAV(offsets[localElement], ncols),
                              justOneA(0, ncols));
    }// *** source loop ***

    //Fill-complete adjs Graph
    adjsMatrix->fillComplete (oneToOneTMap, //throughMapG,
                              adjsMatrix->getRowMap());

    // Form 2ndAdjs
    RCP<sparse_matrix_type> secondAdjs =
      rcp (new sparse_matrix_type(adjsMatrix->getRowMap(),0));
    Tpetra::MatrixMatrix::Multiply(*adjsMatrix,false,*adjsMatrix,
                                     true,*secondAdjs);
    return secondAdjs;
  }
  return RCP<sparse_matrix_type>();
}

template <typename User>
void get2ndAdjsViewFromAdjs(
    const Teuchos::RCP<const MeshAdapter<User> > &ia,
    const RCP<const Comm<int> > comm,
    Zoltan2::MeshEntityType sourcetarget,
    Zoltan2::MeshEntityType through,
    Teuchos::ArrayRCP<typename MeshAdapter<User>::offset_t> &offsets,
    Teuchos::ArrayRCP<typename MeshAdapter<User>::gno_t> &adjacencyIds)
{
  typedef typename MeshAdapter<User>::gno_t gno_t;
  typedef typename MeshAdapter<User>::lno_t lno_t;
  typedef typename MeshAdapter<User>::offset_t offset_t;
  typedef typename MeshAdapter<User>::node_t node_t;

  typedef int nonzero_t;  // adjacency matrix doesn't need scalar_t
  typedef Tpetra::CrsMatrix<nonzero_t,lno_t,gno_t,node_t>   sparse_matrix_type;

  RCP<sparse_matrix_type> secondAdjs = get2ndAdjsMatFromAdjs(ia,comm,sourcetarget,through);

  /* Allocate memory necessary for the adjacency */
  size_t LocalNumIDs = ia->getLocalNumOf(sourcetarget);
  lno_t *start = new lno_t [LocalNumIDs+1];

  if (secondAdjs!=RCP<sparse_matrix_type>()) {

    size_t nadj = 0;

    gno_t const *Ids=NULL;
    ia->getIDsViewOf(sourcetarget, Ids);

    std::vector<gno_t> adj;

    for (size_t localElement=0; localElement<LocalNumIDs; ++localElement){
      start[localElement] = nadj;
      const gno_t globalRow = Ids[localElement];
      size_t NumEntries = secondAdjs->getNumEntriesInGlobalRow (globalRow);
      typename sparse_matrix_type::nonconst_global_inds_host_view_type  Indices("Indices", NumEntries);
      typename sparse_matrix_type::nonconst_values_host_view_type Values("Values", NumEntries);
      secondAdjs->getGlobalRowCopy (globalRow,Indices,Values,NumEntries);

      for (size_t j = 0; j < NumEntries; ++j) {
        if(globalRow != Indices[j]) {
          adj.push_back(Indices[j]);
          nadj++;;
        }
      }
    }

    Ids = NULL;
    start[LocalNumIDs] = nadj;

    gno_t *adj_ = new gno_t [nadj];

    for (size_t i=0; i < nadj; i++) {
      adj_[i] = adj[i];
    }

    offsets = arcp<offset_t>(start, 0, LocalNumIDs+1, true);
    adjacencyIds = arcp<gno_t>(adj_, 0, nadj, true);
  }
  else {
    // No adjacencies could be computed; return no edges and valid offsets array
    for (size_t i = 0; i <= LocalNumIDs; i++)
      start[i] = 0;

    offsets = arcp<offset_t>(start, 0, LocalNumIDs+1, true);
    adjacencyIds = Teuchos::null;
  }

  //return nadj;
}

}

#endif
