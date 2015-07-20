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

/*! \file Zoltan2_ModelHelpers.hpp
    \brief Defines helper functions for use in the models
*/

#include <Zoltan2_MeshAdapter.hpp>

#ifndef _ZOLTAN2_MODELHELPERS_HPP_
#define _ZOLTAN2_MODELHELPERS_HPP_

namespace Zoltan2 {
template <typename User>
void get2ndAdjsViewFromAdjs(const Teuchos::RCP<const MeshAdapter<User> > &ia,
                            Zoltan2::MeshEntityType sourcetarget,
                            Zoltan2::MeshEntityType through,
                            const typename MeshAdapter<User>::lno_t *&offsets,
                            const typename MeshAdapter<User>::zgid_t *&adjacencyIds)
{
  typedef typename MeshAdapter<User>::zgid_t zgid_t;
  typedef typename MeshAdapter<User>::gno_t gno_t;
  typedef typename MeshAdapter<User>::lno_t lno_t;
  typedef typename MeshAdapter<User>::node_t node_t;

  typedef int nonzero_t;  // adjacency matrix doesn't need scalar_t
  typedef Tpetra::CrsMatrix<nonzero_t,lno_t,gno_t,node_t>   sparse_matrix_type;
  typedef Tpetra::Map<lno_t, gno_t, node_t>                 map_type;
  //typedef Tpetra::global_size_t GST;
  //const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

  /* Find the adjacency for a nodal based decomposition */
  size_t nadj = 0;
  if (ia->availAdjs(sourcetarget, through)) {
    using Tpetra::DefaultPlatform;
    using Teuchos::Array;
    using Teuchos::as;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Get the default communicator and Kokkos Node instance
    // TODO:  Default communicator is not correct here; need to get
    // TODO:  communicator from the problem
    RCP<const Comm<int> > comm =
      DefaultPlatform::getDefaultPlatform ().getComm ();

    // Get node-element connectivity

    offsets=NULL;
    adjacencyIds=NULL;
    ia->getAdjsView(sourcetarget, through, offsets, adjacencyIds);

    zgid_t const *Ids=NULL;
    ia->getIDsViewOf(sourcetarget, Ids);

    zgid_t const *throughIds=NULL;
    ia->getIDsViewOf(through, throughIds);

    size_t LocalNumIDs = ia->getLocalNumOf(sourcetarget);

    /***********************************************************************/
    /************************* BUILD MAPS FOR ADJS *************************/
    /***********************************************************************/

    Array<gno_t> sourcetargetGIDs;
    RCP<const map_type> sourcetargetMapG;
    RCP<const map_type> throughMapG;

    // count owned nodes
    size_t LocalNumOfThrough = ia->getLocalNumOf(through);

    // Build a list of the global sourcetarget ids...
    sourcetargetGIDs.resize (LocalNumIDs);
    gno_t min[2];
    min[0] = as<gno_t> (Ids[0]);
    for (size_t i = 0; i < LocalNumIDs; ++i) {
      sourcetargetGIDs[i] = as<gno_t> (Ids[i]);

      if (sourcetargetGIDs[i] < min[0]) {
	min[0] = sourcetargetGIDs[i];
      }
    }

    // min(throughIds[i])
    min[1] = as<gno_t> (throughIds[0]);
    for (size_t i = 0; i < LocalNumOfThrough; ++i) {
      gno_t tmp = as<gno_t> (throughIds[i]);

      if (tmp < min[1]) {
	min[1] = tmp;
      }
    }

    gno_t gmin[2];
    Teuchos::reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_MIN, 2, min, gmin);

    //Generate Map for sourcetarget.
    sourcetargetMapG = rcp(new map_type(ia->getGlobalNumOf(sourcetarget),
					sourcetargetGIDs(), gmin[0], comm));

    //Generate Map for through.
// TODO
// TODO Could check for max through id as well, and if all through ids are
// TODO in gmin to gmax, then default constructors works below.
// TODO Otherwise, may need a constructor that is not one-to-one containing
// TODO all through entities on processor, followed by call to createOneToOne
// TODO

    throughMapG = rcp (new map_type(ia->getGlobalNumOf(through),gmin[1],comm));

    /***********************************************************************/
    /************************* BUILD GRAPH FOR ADJS ************************/
    /***********************************************************************/

    RCP<sparse_matrix_type> adjsMatrix;

    // Construct Tpetra::CrsGraph objects.
    adjsMatrix = rcp (new sparse_matrix_type (sourcetargetMapG, 0));

    nonzero_t justOne = 1;
    ArrayView<nonzero_t> justOneAV = Teuchos::arrayView (&justOne, 1);

    for (size_t localElement=0; localElement<LocalNumIDs; ++localElement){

      //globalRow for Tpetra Graph
      gno_t globalRowT = as<gno_t> (Ids[localElement]);

// KDD can we insert all adjacencies at once instead of one at a time
// (since they are contiguous in adjacencyIds)?
// KDD maybe not until we get rid of zgid_t, as we need the conversion to gno_t.
      for (lno_t j=offsets[localElement]; j<offsets[localElement+1]; ++j){
	gno_t globalCol = as<gno_t> (adjacencyIds[j]);
	//create ArrayView globalCol object for Tpetra
	ArrayView<gno_t> globalColAV = Teuchos::arrayView (&globalCol,1);

	//Update Tpetra adjs Graph
	adjsMatrix->insertGlobalValues(globalRowT,globalColAV,justOneAV);
      }// *** through loop ***
    }// *** source loop ***

    //Fill-complete adjs Graph
    adjsMatrix->fillComplete (throughMapG, adjsMatrix->getRowMap());

    // Form 2ndAdjs
    RCP<sparse_matrix_type> secondAdjs =
      rcp (new sparse_matrix_type(adjsMatrix->getRowMap(),0));
    Tpetra::MatrixMatrix::Multiply(*adjsMatrix,false,*adjsMatrix,
                                     true,*secondAdjs);
    Array<gno_t> Indices;
    Array<nonzero_t> Values;

    /* Allocate memory necessary for the adjacency */
    lno_t *start = new lno_t [LocalNumIDs+1];
    std::vector<gno_t> adj;

    for (size_t localElement=0; localElement<LocalNumIDs; ++localElement){
      start[localElement] = nadj;
      const gno_t globalRow = Ids[localElement];
      size_t NumEntries = secondAdjs->getNumEntriesInGlobalRow (globalRow);
      Indices.resize (NumEntries);
      Values.resize (NumEntries);
      secondAdjs->getGlobalRowCopy (globalRow,Indices(),Values(),NumEntries);

      for (size_t j = 0; j < NumEntries; ++j) {
	if(globalRow != Indices[j]) {
	  adj.push_back(Indices[j]);
	  nadj++;;
	}
      }
    }

    Ids = NULL;
    start[LocalNumIDs] = nadj;

    zgid_t *adj_ = new zgid_t [nadj];

    for (size_t i=0; i < nadj; i++) {
      adj_[i] = adj[i];
    }

    offsets = start;
    adjacencyIds = adj_;
  }

  //return nadj;
}

}

#endif
