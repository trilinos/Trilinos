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
#ifndef _ZOLTAN2_ALGZOLTANCALLBACKS_HPP_
#define _ZOLTAN2_ALGZOLTANCALLBACKS_HPP_

#include <Zoltan2_MeshAdapter.hpp>
#include <Zoltan2_VectorAdapter.hpp>
#include <Zoltan2_GraphAdapter.hpp>
#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_IdentifierAdapter.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_HyperGraphModel.hpp>

#include <Zoltan2_MachineRepresentation.hpp>

#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>
#include <zoltan_cpp.h>

extern "C" {
#include <zz_const.h>
#include <zoltan_mem.h>
}

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgZoltanCallbacks.hpp
//! \brief callback functions for the Zoltan package (templated on Adapter)
//  Callbacks based on Adapter; specializations provided where needed

namespace Zoltan2 {

/////////////////////////////////////////////////////////////////////////////
// CALLBACKS SHARED BY MANY ADAPTERS
/////////////////////////////////////////////////////////////////////////////

////////////////////
// ZOLTAN_NUM_OBJ_FN
template <typename Adapter>
static int zoltanNumObj(void *data, int *ierr) {
  const Adapter *adp = static_cast<Adapter *>(data);
  *ierr = ZOLTAN_OK;
  return int(adp->getLocalNumIDs());
}

/////////////////////
// ZOLTAN_OBJ_LIST_FN
template <typename Adapter>
static void zoltanObjList(void *data, int nGidEnt, int nLidEnt,
                          ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                          int wdim, float *wgts, int *ierr)
{
  const Adapter *adp = static_cast<Adapter *>(data);
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  *ierr = ZOLTAN_OK;

  size_t mynObj = adp->getLocalNumIDs();

  const gno_t *myids = NULL;
  adp->getIDsView(myids);
  for (size_t i = 0; i < mynObj; i++) {
    ZOLTAN_ID_PTR idPtr = &(gids[i*nGidEnt]);
    TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, myids[i]);
    idPtr = &(lids[i*nLidEnt]);
    TPL_Traits<ZOLTAN_ID_PTR,lno_t>::ASSIGN(idPtr, lno_t(i));
  }

  if (wdim) {
    int mywdim = adp->getNumWeightsPerID();
    for (int w = 0; w < wdim; w++) {
      if (w < mywdim) {
        // copy weights from adapter
        const typename Adapter::scalar_t *mywgts;
        int mystride;
        adp->getWeightsView(mywgts, mystride, w);
        for (size_t i = 0; i < mynObj; i++)
          wgts[i*wdim+w] = float(mywgts[i*mystride]);
      }
      else {
        // provide uniform weights
        for (size_t i = 0; i < mynObj; i++)
          wgts[i*wdim+w] = 1.;
      }
    }
  }
}

///////////////////////
// ZOLTAN_PART_MULTI_FN
template <typename Adapter>
static void zoltanParts(void *data, int /* nGidEnt */, int nLidEnt, int nObj,
                        ZOLTAN_ID_PTR /* gids */, ZOLTAN_ID_PTR lids,
                        int *parts, int *ierr)
{
  typedef typename Adapter::lno_t lno_t;
  const Adapter *adp = static_cast<Adapter *>(data);
  *ierr = ZOLTAN_OK;
  const typename Adapter::part_t *myparts;
  adp->getPartsView(myparts);
  // User parts from input adapter
  for (int i = 0; i < nObj; i++) {
    lno_t lidx;
    TPL_Traits<lno_t,ZOLTAN_ID_PTR>::ASSIGN(lidx, &(lids[i*nLidEnt]));
    parts[i] = int(myparts[lidx]);
  }
}

/////////////////////
// ZOLTAN_NUM_GEOM_FN
template <typename Adapter>
static int zoltanNumGeom(void *data, int *ierr)
{
  const Adapter *adp = static_cast<Adapter *>(data);
  *ierr = ZOLTAN_OK;
  return adp->getDimension();
}

///////////////////////
// ZOLTAN_GEOM_MULTI_FN
template <typename Adapter>
static void zoltanGeom(void *data, int /* nGidEnt */, int nLidEnt, int nObj,
                       ZOLTAN_ID_PTR /* gids */, ZOLTAN_ID_PTR lids,
                       int nDim, double *coords, int *ierr)
{
  typedef typename Adapter::lno_t lno_t;
  const Adapter *adp = static_cast<Adapter *>(data);
  *ierr = ZOLTAN_OK;

  for (int d = 0; d < nDim; d++) {
    const typename Adapter::scalar_t *mycoords;
    int mystride;
    adp->getCoordinatesView(mycoords, mystride, d);
    for (int i = 0; i < nObj; i++) {
      lno_t lidx;
      TPL_Traits<lno_t,ZOLTAN_ID_PTR>::ASSIGN(lidx, &(lids[i*nLidEnt]));
      coords[i*nDim+d] = double(mycoords[lidx*mystride]);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// HYPERGRAPH CALLBACKS USING A GRAPH ADAPTER
// Building the most straightforward hypergraph from a graph and, thus,
// avoiding use of HypergraphModel.
// Assuming one hyperedge per vertex, containing vertex and all its nbors
/////////////////////////////////////////////////////////////////////////////

///////////////////////
// ZOLTAN_HG_SIZE_CS_FN
template <typename Adapter>
static void zoltanHGSizeCS_withGraphAdapter(void *data,
                                            int *nLists, int *nPins,
                                            int *format, int *ierr
)
{
  // Assuming one hyperedge per vertex consisting of vertex and its graph nbors.
  const Adapter *adp = static_cast<Adapter *>(data);
  *ierr = ZOLTAN_OK;

  *nLists = Teuchos::as<int>(adp->getLocalNumVertices());
  *nPins = Teuchos::as<int>(adp->getLocalNumEdges()+adp->getLocalNumVertices());
           // number of given edges + self pin for each vertex
  *format = ZOLTAN_COMPRESSED_EDGE;
}

//////////////////
// ZOLTAN_HG_CS_FN
template <typename Adapter>
static void zoltanHGCS_withGraphAdapter(void *data, int nGidEnt, int nLists,
                                        int /* nPins */, int /* format */,
                                        ZOLTAN_ID_PTR listIds, int *listIdx,
                                        ZOLTAN_ID_PTR pinIds, int *ierr
)
{
  // Assuming one hyperedge per vertex consisting of vertex and its graph nbors.
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  const GraphAdapter<user_t, userCoord_t>* adp =
        static_cast<GraphAdapter<user_t, userCoord_t>* >(data);

  *ierr = ZOLTAN_OK;

  const gno_t *ids;
  const gno_t *adjIds;
  const offset_t *offsets;

  try {
    adp->getIDsView(ids);
    adp->getEdgesView(offsets, adjIds);
  }
  catch (std::exception &e) {
    *ierr = ZOLTAN_FATAL;
  }

  if (*ierr == ZOLTAN_OK) {
    // copy into Zoltan's memory
    for (int i=0; i < nLists; i++) {
      ZOLTAN_ID_PTR idPtr = &(listIds[i*nGidEnt]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, ids[i]);
      listIdx[i] = Teuchos::as<int>(offsets[i]+i);  // adding self pin
    }
    listIdx[nLists] = Teuchos::as<int>(offsets[nLists]);
    int pinCnt = 0;
    for (int i=0; i < nLists; i++) {
      ZOLTAN_ID_PTR idPtr = &(pinIds[pinCnt*nGidEnt]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, ids[i]);
      pinCnt++;
      for (offset_t j = offsets[i]; j < offsets[i+1]; j++) {
        idPtr = &(pinIds[pinCnt*nGidEnt]);
        TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, adjIds[j]);
        pinCnt++;
      }
    }
  }
}

//////////////////////////////
// ZOLTAN_HG_SIZE_EDGE_WGTS_FN
template <typename Adapter>
static void zoltanHGSizeEdgeWts_withGraphAdapter(
  void *data,
  int *nEdges,
  int *ierr
)
{
  // Assuming one hyperedge per vertex consisting of vertex and its graph nbors.
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  const GraphAdapter<user_t, userCoord_t>* adp =
        static_cast<GraphAdapter<user_t, userCoord_t>* >(data);
  *ierr = ZOLTAN_OK;
  *nEdges = Teuchos::as<int>(adp->getLocalNumVertices()); // one edge per vertex
}

//////////////////////////////
// ZOLTAN_HG_EDGE_WGTS_FN
template <typename Adapter>
static void zoltanHGEdgeWts_withGraphAdapter(
  void *data,
  int nGidEnt,
  int nLidEnt,
  int nEdges,
  int /* eWgtDim */,
  ZOLTAN_ID_PTR edgeGids,
  ZOLTAN_ID_PTR edgeLids,
  float *edgeWgts,
  int *ierr
)
{
  // Assuming one hyperedge per vertex consisting of vertex and its graph nbors.
  // Hyperedge weight is then sum of edge weights to nbors.
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  const GraphAdapter<user_t, userCoord_t>* adp =
        static_cast<GraphAdapter<user_t, userCoord_t>* >(data);

  *ierr = ZOLTAN_OK;

  const gno_t *ids;
  const gno_t *adjIds;
  const offset_t *offsets;
  const scalar_t *ewgts;
  int stride;
  try {
    adp->getIDsView(ids);
    adp->getEdgesView(offsets, adjIds);
    adp->getEdgeWeightsView(ewgts, stride, 0);  // Use only first weight
  }
  catch (std::exception &e) {
    *ierr = ZOLTAN_FATAL;
  }
  if (ierr == ZOLTAN_OK) {
    for (int i = 0; i < nEdges; i++) {
      float sum = 0;
      for (offset_t j = offsets[i]; j < offsets[i+1]; j++)
        sum += ewgts[j*stride];
      ZOLTAN_ID_PTR idPtr = &(edgeGids[i*nGidEnt]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, ids[i]);
      if (nLidEnt) {
        idPtr = &(edgeLids[i*nLidEnt]);
        TPL_Traits<ZOLTAN_ID_PTR,int>::ASSIGN(idPtr, i);
      }
      edgeWgts[i] = sum;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// HYPERGRAPH CALLBACKS USING A MATRIX ADAPTER
// Building the most straightforward hypergraph from a matrix and, thus,
// avoiding use of HypergraphModel.
// Assuming vertices are rows or columns, and pins are nonzeros.
/////////////////////////////////////////////////////////////////////////////

///////////////////////
// ZOLTAN_HG_SIZE_CS_FN
template <typename Adapter>
static void zoltanHGSizeCS_withMatrixAdapter(void *data,
                                             int *nLists, int *nPins,
                                             int *format, int *ierr
)
{
  *ierr = ZOLTAN_OK;
  typedef typename Adapter::user_t user_t;
  const MatrixAdapter<user_t>* madp = static_cast<MatrixAdapter<user_t>* >(data);

  *nPins = madp->getLocalNumEntries();

  MatrixEntityType etype = madp->getPrimaryEntityType();
  if (etype == MATRIX_ROW && madp->CRSViewAvailable()) {
    *nLists = madp->getLocalNumRows();
    *format = ZOLTAN_COMPRESSED_VERTEX;
  }
  else if (etype == MATRIX_ROW && madp->CCSViewAvailable()) {
    *nLists = madp->getLocalNumColumns();
    *format = ZOLTAN_COMPRESSED_EDGE;
  }
  else if (etype == MATRIX_COLUMN && madp->CRSViewAvailable()) {
    *nLists = madp->getLocalNumRows();
    *format = ZOLTAN_COMPRESSED_EDGE;
  }
  else if (etype == MATRIX_COLUMN && madp->CCSViewAvailable()) {
      *nLists = madp->getLocalNumColumns();
      *format = ZOLTAN_COMPRESSED_VERTEX;
  }
  else {
    // Need either CRSView or CCSView.
    // Also, not yet implemented for matrix nonzeros;
    // may need a hypergraph model.
    std::cout << "For hypergraph partitioning, "
              << "CRSView or CCSView is needed in MatrixAdapter" << std::endl;
    *ierr = ZOLTAN_FATAL;
  }
}

//////////////////
// ZOLTAN_HG_CS_FN
template <typename Adapter>
static void zoltanHGCS_withMatrixAdapter(void *data, int nGidEnt, int nLists,
                                         int nPins, int format,
                                         ZOLTAN_ID_PTR listIds, int *listIdx,
                                         ZOLTAN_ID_PTR pinIds, int *ierr
)
{
  *ierr = ZOLTAN_OK;
  typedef typename Adapter::gno_t gno_t;
  // typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::user_t user_t;
  const MatrixAdapter<user_t>* madp = static_cast<MatrixAdapter<user_t>* >(data);

  const gno_t *Ids;
  const gno_t *pIds;
  const offset_t *offsets;

  // Get the pins and list IDs.
  if (madp->CRSViewAvailable()) {
    try {
      madp->getRowIDsView(Ids);
      madp->getCRSView(offsets, pIds);
    }
    catch (std::exception &e) {
      *ierr = ZOLTAN_FATAL;
    }
  }
  else if (madp->CCSViewAvailable()) {
    try {
      madp->getColumnIDsView(Ids);
      madp->getCCSView(offsets, pIds);
    }
    catch (std::exception &e) {
      *ierr = ZOLTAN_FATAL;
    }
  }
  else {
    // Need either CRSView or CCSView.
    *ierr = ZOLTAN_FATAL;
  }

  if (*ierr == ZOLTAN_OK) {
    // copy into Zoltan's memory
    for (int i=0; i < nLists; i++) {
      ZOLTAN_ID_PTR idPtr = &(listIds[i*nGidEnt]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, Ids[i]);
      listIdx[i] = Teuchos::as<int>(offsets[i]);
    }
    listIdx[nLists] = Teuchos::as<int>(offsets[nLists]);
    for (int i=0; i < nPins; i++) {
      ZOLTAN_ID_PTR idPtr = &(pinIds[i*nGidEnt]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, pIds[i]);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// TODO:  GRAPH ADAPTER CALLBACKS

/////////////////////////////////////////////////////////////////////////////
// HYPERGRAPH CALLBACKS USING A MESH ADAPTER
// Implement Boman/Chevalier's hypergraph mesh model
// Skip explicit construction of a HypergraphModel
// Return either (depending on available adjacencies):
// +  for each primary entity (vtx), the list of assoc adjacency entities (edges)
// +  for each adjacency entity (edge), the list of assoc primary entities (vtx)
/////////////////////////////////////////////////////////////////////////////

///////////////////////
// ZOLTAN_HG_SIZE_CS_FN
template <typename Adapter>
static void zoltanHGSizeCS_withMeshAdapter(void *data, int *nLists, int *nPins,
                                           int *format, int *ierr
)
{
  *ierr = ZOLTAN_OK;
  typedef typename Adapter::user_t user_t;
  const MeshAdapter<user_t>* madp = static_cast<MeshAdapter<user_t>* >(data);
  if (madp->availAdjs(madp->getPrimaryEntityType(),
                      madp->getAdjacencyEntityType()))
  {
    *nLists = madp->getLocalNumOf(madp->getPrimaryEntityType());
    *nPins = madp->getLocalNumAdjs(madp->getPrimaryEntityType(),
                                   madp->getAdjacencyEntityType());
    *format = ZOLTAN_COMPRESSED_VERTEX;
  }
  else if (madp->availAdjs(madp->getAdjacencyEntityType(),
                           madp->getPrimaryEntityType()))
  {
    *nLists = madp->getLocalNumOf(madp->getAdjacencyEntityType());
    *nPins = madp->getLocalNumAdjs(madp->getAdjacencyEntityType(),
                                   madp->getPrimaryEntityType());
    *format = ZOLTAN_COMPRESSED_EDGE;
  }
  else {
    std::cout << "For hypergraph partitioning, need first adjacencies "
              << "(availAdjs, getLocalNumAdjs, getAdjsView) "
              << "in MeshAdapter."  << std::endl;
    *nLists = 0;
    *nPins = 0;
    *format = -1*ZOLTAN_COMPRESSED_VERTEX;
    *ierr = ZOLTAN_FATAL;
  }
}

//////////////////
// ZOLTAN_HG_CS_FN
template <typename Adapter>
static void zoltanHGCS_withMeshAdapter(
  void *data, int nGidEnt, int nLists, int nPins,
  int format, ZOLTAN_ID_PTR listIds,
  int *listIdx, ZOLTAN_ID_PTR pinIds, int *ierr
)
{
  *ierr = ZOLTAN_OK;
  typedef typename Adapter::gno_t gno_t;
  // typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::offset_t offset_t;
  const MeshAdapter<user_t>* madp = static_cast<MeshAdapter<user_t>*>(data);

  // Select listType and pinType based on format specified in ZOLTAN_HG_CS_SIZE_FN
  MeshEntityType listType, pinType;
  if (format == ZOLTAN_COMPRESSED_VERTEX)
  {
    listType = madp->getPrimaryEntityType();
    pinType = madp->getAdjacencyEntityType();
  }
  else if (format == ZOLTAN_COMPRESSED_EDGE)
  {
    listType = madp->getAdjacencyEntityType();
    pinType = madp->getPrimaryEntityType();
  }
  else {
    *ierr = ZOLTAN_FATAL;
  }

  if (*ierr == ZOLTAN_OK) {

    // get list IDs
    const gno_t *Ids;
    try {
      madp->getIDsViewOf(listType,Ids);
    }
    catch (std::exception &e) {
      *ierr = ZOLTAN_FATAL;
    }

    // get pins
    const offset_t* offsets;
    const gno_t* adjIds;
    try {
      madp->getAdjsView(listType, pinType, offsets, adjIds);
    }
    catch (std::exception &e) {
      *ierr = ZOLTAN_FATAL;
    }

    // copy into Zoltan's memory
    for (int i=0; i < nLists; i++) {
      ZOLTAN_ID_PTR idPtr = &(listIds[i*nGidEnt]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, Ids[i]);
      listIdx[i] = Teuchos::as<int>(offsets[i]);
    }
    listIdx[nLists] = Teuchos::as<int>(offsets[nLists]);
    for (int i=0; i < nPins; i++) {
      ZOLTAN_ID_PTR idPtr = &(pinIds[i*nGidEnt]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, adjIds[i]);
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// HYPERGRAPH CALLBACKS FROM A HYPERGRAPH MODEL
/////////////////////////////////////////////////////////////////////////////

////////////////////
// ZOLTAN_NUM_OBJ_FN
template <typename Adapter>
static int zoltanHGNumObj_withModel(void *data, int *ierr) {
  const HyperGraphModel<Adapter>* mdl =
                                  static_cast<HyperGraphModel<Adapter>* >(data);
  *ierr = ZOLTAN_OK;
  return int(mdl->getLocalNumOwnedVertices());
}

/////////////////////
// ZOLTAN_OBJ_LIST_FN
template <typename Adapter>
static void zoltanHGObjList_withModel(void *data, int nGidEnt, int nLidEnt,
                                      ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                                      int wdim, float *wgts, int *ierr)
{
  const HyperGraphModel<Adapter>* mdl =
                                  static_cast<HyperGraphModel<Adapter>* >(data);
  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::scalar_t    scalar_t;
  typedef StridedData<lno_t, scalar_t>  input_t;

  *ierr = ZOLTAN_OK;
  ArrayView<const gno_t> Ids;
  ArrayView<input_t> model_wgts;
  size_t num_verts = mdl->getVertexList(Ids,model_wgts);
  ArrayView<bool> isOwner;
  mdl->getOwnedList(isOwner);
  int j=0;
  for (size_t i=0;i<num_verts;i++) {
    if (isOwner[i]) {
      ZOLTAN_ID_PTR idPtr = &(gids[j*nGidEnt]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, Ids[i]);
      idPtr = &(lids[j*nLidEnt]);
      TPL_Traits<ZOLTAN_ID_PTR,lno_t>::ASSIGN(idPtr, lno_t(i));
      j++;
    }
  }
  if (wdim) {
    int mywdim = mdl->getNumWeightsPerVertex();
    for (int w = 0; w < wdim; w++) {
      j=0;
      if (w < mywdim) {
        for (size_t i = 0; i < num_verts; i++)  {
          if (isOwner[i]) {
            wgts[j*wdim+w] = float(model_wgts[w][i]);
            j++;
          }
        }
      }
      else {
        // provide uniform weights
        for (size_t i = 0; i < num_verts; i++) {
          if (isOwner[i]) {
            wgts[j*wdim+w] = 1.;
            j++;
          }
        }
      }
    }
  }
}

///////////////////////
// ZOLTAN_HG_SIZE_CS_FN
template <typename Adapter>
static void zoltanHGSizeCS_withModel(void *data, int *nEdges, int *nPins,
                                     int *format, int *ierr
)
{
  *ierr = ZOLTAN_OK;
  const HyperGraphModel<Adapter>* mdl =
                                  static_cast<HyperGraphModel<Adapter>* >(data);
  *nEdges = mdl->getLocalNumHyperEdges();
  *nPins = mdl->getLocalNumPins();
  if (mdl->getCentricView()==HYPEREDGE_CENTRIC)
    *format = ZOLTAN_COMPRESSED_EDGE;
  else
    *format = ZOLTAN_COMPRESSED_VERTEX;
}

//////////////////
// ZOLTAN_HG_CS_FN
template <typename Adapter>
static void zoltanHGCS_withModel(void *data, int nGidEnt, int nEdges, int nPins,
                                 int format, ZOLTAN_ID_PTR edgeIds,
                                 int *edgeIdx, ZOLTAN_ID_PTR pinIds, int *ierr
)
{
  *ierr = ZOLTAN_OK;
  const HyperGraphModel<Adapter>* mdl =
                                  static_cast<HyperGraphModel<Adapter>* >(data);
  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::offset_t    offset_t;
  typedef typename Adapter::scalar_t    scalar_t;
  typedef StridedData<lno_t, scalar_t>  input_t;

  ArrayView<const gno_t> Ids;
  ArrayView<input_t> wgts;
  mdl->getEdgeList(Ids,wgts);
  ArrayView<const gno_t> pinIds_;
  ArrayView<const offset_t> offsets;
  ArrayView<input_t> pin_wgts;
  mdl->getPinList(pinIds_,offsets,pin_wgts);
  for (int i=0;i<nEdges;i++) {
    ZOLTAN_ID_PTR idPtr = &(edgeIds[i*nGidEnt]);
    TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, Ids[i]);
    edgeIdx[i] = Teuchos::as<int>(offsets[i]);
  }

  for (int i=0;i<nPins;i++) {
    ZOLTAN_ID_PTR idPtr = &(pinIds[i*nGidEnt]);
    TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, pinIds_[i]);
  }
}

////////////////////////////
// ZOLTAN_NUM_EDGES_MULTI_FN
template <typename Adapter>
static void zoltanNumEdgesMulti_withGraphModel(void *data,
                                                 int nGidEnt,
                                                 int nLigEnt,
                                                 int nObjs,
                                                 ZOLTAN_ID_PTR gids,
                                                 ZOLTAN_ID_PTR lids,
                                                 int *nEdges,
                                                 int *ierr)
{
/*  // Graph Adapter
  const Adapter *adp = static_cast<Adapter *>(data);
  *ierr = ZOLTAN_OK;

  std::cout << "\n\nHello NUM EDGES MULTI" << std::endl;


  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::offset_t    offset_t;
  typedef typename Adapter::scalar_t    scalar_t;

  ArrayView<const offset_t> offsets;
  ArrayView<const gno_t> adjIds;

  adp->getEdgesView(offsets, adjIds);

  std::cout << "\nOffsets length: " << offsets."

  for (int i = 0; i < nObjs; ++i) {

    int obj_edges = Teuchos::as<int>(offsets[i+1]) - Teuchos::as<int>(offsets[i]);

    nEdges[i] = obj_edges;
  }

  // Output is nEdges

  return;
*/


  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::offset_t    offset_t;
  typedef typename Adapter::scalar_t    scalar_t;
//  typedef typename Adapter::input_t     input_t;
  typedef typename Adapter::user_t      user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  typedef StridedData<lno_t, scalar_t>  input_t;

  const GraphModel<GraphAdapter<user_t, userCoord_t> >* graph_model =
        static_cast<GraphModel<GraphAdapter<user_t, userCoord_t>>* >(data);

  *ierr = ZOLTAN_OK;

  ArrayView<const gno_t>    adjIds;
  ArrayView<const offset_t> offsets;
  ArrayView<input_t>        ewgts;
  size_t nLocalEdges;

  nLocalEdges = graph_model->getEdgeList(adjIds, offsets, ewgts);

  int tot_edges = 0;

  for (int i = 0; i < nObjs; ++i) {
    nEdges[i] = Teuchos::as<int>(offsets[i + 1] - offsets[i]);

    tot_edges += nEdges[i];

  }

  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  std::cout << "\nRank: " << myrank << ", Local NUM_EDGES: " << tot_edges << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);

/*
    for (int i = 0; i < nEdges; i++) {
      float sum = 0;
      for (offset_t j = offsets[i]; j < offsets[i+1]; j++)
        sum += ewgts[j*stride];
      ZOLTAN_ID_PTR idPtr = &(edgeGids[i*nGidEnt]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, ids[i]);
      if (nLidEnt) {
        idPtr = &(edgeLids[i*nLidEnt]);
        TPL_Traits<ZOLTAN_ID_PTR,int>::ASSIGN(idPtr, i);
      }
      edgeWgts[i] = sum;
    }
  }
*/

}



////////////////////////////
// ZOLTAN_EDGE_LIST_MULTI_FN
template <typename Adapter>
static void zoltanEdgeListMulti_withGraphModel(void *data,
                                                 int nGidEnt,
                                                 int nLigEnt,
                                                 int nObjs,
                                                 ZOLTAN_ID_PTR gids,
                                                 ZOLTAN_ID_PTR lids,
                                                 int *nEdges,
                                                 ZOLTAN_ID_PTR nbor_gids,
                                                 int *nbor_procs,
                                                 int wdim,
                                                 float *nbor_ewgts,
                                                 int *ierr)
{



  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::offset_t    offset_t;
  typedef typename Adapter::scalar_t    scalar_t;
  typedef typename Adapter::user_t      user_t;
  typedef typename Adapter::userCoord_t userCoord_t;

  typedef StridedData<lno_t, scalar_t>  input_t;

  const GraphModel<GraphAdapter<user_t, userCoord_t> >* graph_model =
        static_cast<GraphModel<GraphAdapter<user_t, userCoord_t> >* >(data);

  *ierr = ZOLTAN_OK;

  ArrayView<const gno_t> Ids;
  ArrayView<input_t> vwgts;

  size_t nLocalVerts = graph_model->getVertexList(Ids, vwgts);

//  std::cout << "\nNLocalVerts: " << nLocalVerts << " Ids: ";
//  for (int i = 0; i < Teuchos::as<int>(nLocalVerts); ++i) {
//      std::cout << " " << Ids[i];
//  }
//  std::cout << std::endl;

  ArrayView<const gno_t> adjIds;
  ArrayView<const offset_t> offsets;
  ArrayView<input_t> ewgts;


  size_t nLocalEdges = graph_model->getEdgeList(adjIds, offsets, ewgts);



  ArrayView<size_t> vtxdist;

  // vtxdist size: comm->getSize() + 1
  graph_model->getVertexDist(vtxdist);

//  std::cout << "\nVTXDIST: ";
//  for (int i = 0; i < vtxdist.size(); ++i) {
//    std::cout << " " << Teuchos::as<int>(vtxdist[i]);
//  }
//  std::cout << std::endl;


  int edge_idx = 0;
  int proc_idx = 0;

  std::cout << "Size vtxdist: " << vtxdist.size() << std::endl;

//  std::cout << "nObjs: " << nObjs << std::endl;

//  std::cout << "NELEM: ";
//  for (int i = 0; i < nObjs; ++i) {
//    std::cout << " " << nEdges[i];
//  }
//  std::cout << std::endl;

  for (int i = 0; i < nObjs; ++i) {

    // If we pass the next offset, we are on to the next proc
//    if (i >= Teuchos::as<int>(vtxdist[proc_idx + 1])) {
//      proc_idx++;
//    }

    for (int j = 0;  j < nEdges[i]; ++j) {

//      if (edge_idx >= Teuchos::as<int>(vtxdist[proc_idx + 1])) {
//        proc_idx++;
//        std::cout << "PROCIDX: " << proc_idx << std::endl;
//      }

      nbor_gids[edge_idx] = Teuchos::as<int>(adjIds[edge_idx]);
//      nbor_procs[edge_idx] = proc_idx;


      for (int k = 0; k < vtxdist.size() - 1; ++k) {

        if (nbor_gids[edge_idx] >= vtxdist[k] &&
            nbor_gids[edge_idx] < vtxdist[k + 1]) {

          proc_idx = k;
        }

      }

      nbor_procs[edge_idx] = proc_idx;



//      std::cout << "NUM_GIDS: e" << edge_idx << " o" << i << " of " << nObjs
//          << " nE" << nEdges[i] << " ng" << nbor_gids[edge_idx] << std::endl;
//      std::cout << "NUM_PROCS: e" << edge_idx << " o" << i << " of " << nObjs
//          << " nE" << nEdges[i] << " NP" << nbor_procs[edge_idx] << std::endl;

      edge_idx++;
    }
  }

  int myrank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  std::cout << "Rank: " << myrank << ", Max EDGE_IDX: " << edge_idx << std::endl;


/*
  int proc_idx = 0;
  std::cout << "nGIDENT: " << nGidEnt << std::endl;

  for (int i = 0; i < nGidEnt; ++i) {

    if (i >= Teuchos::as<int>(vtxdist[proc_idx + 1])) {
      proc_idx++;
    }

    nbor_gids[i] = Teuchos::as<int>(adjIds[i]);
    nbor_procs[i] = proc_idx;
  }
*/

/*
  if (wdim) {
    int mywdim = graph_model->getNumWeightsPerEdge();
    for (int w = 0; w < wdim; w++) {
      int edge_idx = 0;
      for (size_t obj = 0; obj < nObjs; ++obj) {
          if (w < mywdim) {
            for (size_t i = 0; i < nEdges[]; i++)  {
                nbor_ewgts[edge_idx*wdim+w] = Teuchos::as<float>(ewgts[w][i]);
                edge_idx++;
            }
          }
          else {
            // provide uniform weights
            for (size_t i = 0; i < num_verts; i++) {
                nbor_ewgts[edge_idx*wdim+w] = 1.;
                edge_idx++;
            }
          }
      }
    }
  }
*/

  if (wdim) {

    //std::cout << "wdim" << std::endl;
    int mywdim = graph_model->getNumWeightsPerEdge();

    edge_idx = 0;

    for (int obj = 0; obj < nObjs; ++obj) {

      for (int i = 0; i < nEdges[obj]; ++i) {

        for (int w = 0; w < wdim; ++w) {
          if (w < mywdim) {
            nbor_ewgts[edge_idx + w] = Teuchos::as<float>(ewgts[w][i]);
          }
          else {
            nbor_ewgts[edge_idx + w] = 1.;

          }
//          std::cout << "NUM_EWGTS: e" << edge_idx << " " << w << " with wd" << wdim
//                << " o" << obj << " of " << nObjs
//                << " nE" << nEdges[obj] << " EW" << nbor_ewgts[edge_idx] << std::endl;
        }

        edge_idx += wdim;
      }
    }
  }
//  else {
//    std::cout << "nowdim" << std::endl;
//  }





/*
  // Graph Adapater
  const Adapter *adp = static_cast<Adapter *>(data);
  *ierr = ZOLTAN_OK;
  std::cout << "\n\nHello EDGE LIST MULTI" << std::endl;
  // Output is nbor_gids, nbor_procs, nbor_ewgts
  return;
*/
/*
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  const GraphAdapter<user_t, userCoord_t>* adp =
        static_cast<GraphAdapter<user_t, userCoord_t>* >(data);

  *ierr = ZOLTAN_OK;

  const gno_t *ids;
  const gno_t *adjIds;
  const offset_t *offsets;
  const scalar_t *ewgts;
  int stride;
  try {
    adp->getIDsView(ids);
    adp->getEdgesView(offsets, adjIds);
    adp->getEdgeWeightsView(ewgts, stride, 0);  // Use only first weight
  }
  catch (std::exception &e) {
    *ierr = ZOLTAN_FATAL;
  }
  if (ierr == ZOLTAN_OK) {

    int edge_idx = 0;
    for (int i = 0; i < nObjs; ++i) {
        for (int j = 0; j < nEdges[i]; ++j) {

            nbor_gids[edge_idx] = Teuchos::as<int>(adjIds[edge_idx]);

            nbor_procs = ;

            for (int k = 0; k < stride; ++k) {
                nbor_ewgts[stride * edge_idx + k] = Teuchos::as<float>(ewgts[stride * edge_idx + k]);
            }

            edge_idx++;
        }
    }

*/

/*
    for (int i = 0; i < nEdges; i++) {
      float sum = 0;
      for (offset_t j = offsets[i]; j < offsets[i+1]; j++)
        sum += ewgts[j*stride];
      ZOLTAN_ID_PTR idPtr = &(edgeGids[i*nGidEnt]);
      TPL_Traits<ZOLTAN_ID_PTR,gno_t>::ASSIGN(idPtr, ids[i]);
      if (nLidEnt) {
        idPtr = &(edgeLids[i*nLidEnt]);
        TPL_Traits<ZOLTAN_ID_PTR,int>::ASSIGN(idPtr, i);
      }
      edgeWgts[i] = sum;
    }

  }
*/


}





////////////////////////////
// ZOLTAN_HIER_NUM_LEVELS_FN
template <typename Adapter>
static int zoltanHierNumLevels(void *data,
                               //int /* nGidEnt */, int nLidEnt, int nObj,
                               //ZOLTAN_ID_PTR /* gids */, ZOLTAN_ID_PTR lids,
                               //int *parts,
                               int *ierr)
{


  // Find out how many levels of hierarchy this proc will participate in


  typedef typename Adapter::scalar_t   pcoord_t;
  typedef typename Adapter::part_t     part_t;

  const MachineRepresentation<pcoord_t, part_t> *machine
    = static_cast<MachineRepresentation<pcoord_t, part_t> *>(data);
  *ierr = ZOLTAN_OK;

  int rank = machine->getMyRank();
//  int ranks = machine->getNumRanks();

//  std::cout << "\nRank: " << rank
//    << " HierNumLevels: " << machine->getNumNonuniformLevels()
//    << " (+ 1?)" << std::endl;


  // Group, Subgroup, Rack
  return machine->getNumNonuniformLevels() + 1;
//  return machine->getNumNonuniformLevels();
}

//////////////////////
// ZOLTAN_HIER_PART_FN
template <typename Adapter>
static int zoltanHierPart(void *data,
                           //int /* nGidEnt */, int nLidEnt, int nObj,
                           //ZOLTAN_ID_PTR /* gids */, ZOLTAN_ID_PTR lids,
                           int level,
                           int *ierr)
{


  // Determines parts to compute at this level


  typedef typename Adapter::scalar_t   pcoord_t;
  typedef typename Adapter::part_t     part_t;

  const MachineRepresentation<pcoord_t, part_t> *machine
    = static_cast<MachineRepresentation<pcoord_t, part_t> *>(data);
  *ierr = ZOLTAN_OK;

  int rank = machine->getMyRank();
//  int ranks = machine->getNumRanks();

//  std::cout << "\n-\nRank: " << rank
//    << " zoltanHierPart, level: " << level
//    << "\n-\n"<< std::endl;

  part_t num_unique_groups = machine->getNumUniqueGroups();

//  std::cout << "\nrank: " << rank
//    << " NUG: " << num_unique_groups << std::endl;

//  std::cout << "\nrank: " << rank << ", 1" << std::endl;

//  part_t *group_count = new part_t[num_unique_groups];

  std::vector<part_t> group_count;
  machine->getGroupCount2(group_count);

  int group_idx;

  int upper_cdf = 0;
  int lower_cdf = 0;

  for (int i = 0; i < num_unique_groups; ++i) {

    lower_cdf = upper_cdf;
    upper_cdf += int(group_count[i]);

    if (rank < upper_cdf && rank >= lower_cdf) {

//      std::cout << "\nRank: " << rank << " solving part " << i << std::endl;

      group_idx = i;
      break;
    }
  }


  if (level == 0)
    return group_idx;

  std::vector<part_t> num_unique_subgroups;
  std::vector<std::vector<part_t>> subgroup_counts;

  machine->getNumUniqueSubgroups(num_unique_subgroups);
  machine->getSubgroupCounts(subgroup_counts);


  int subgroup_idx = 0;
  int group_rank = rank - lower_cdf;

  lower_cdf = 0;
  upper_cdf = 0;

  for (int i = 0; i < int(num_unique_subgroups[group_idx]); ++i) {

    lower_cdf = upper_cdf;
    upper_cdf += int(subgroup_counts[group_idx][i]);

    if (group_rank < upper_cdf && group_rank >= lower_cdf) {

//      std::cout << "\nRank: " << rank
//        << " solving part " << i
//        << " of Group " << group_idx
//        << " with group rank " << group_rank << std::endl;

      subgroup_idx = i;
      break;
    }
  }

  if (level == 1)
    return subgroup_idx;



//  for (int i = 0; i < num_unique_groups; ++i) {

//    for (int j = 0; j < num_unique_subgroups[i]; ++j) {
//      std::cout << "\nRank: " << rank
//        << " sgc[" << i << "][" << j << "]: "
//        << subgroup_counts[i][j];
//    }
//  }



//  std::cout << "\nRank: " << rank
//    << " HierPart final level: " << level << std::endl;

  return group_rank - lower_cdf;

}

////////////////////////
// ZOLTAN_HIER_METHOD_FN
template <typename Adapter>
static void zoltanHierMethod(void *data,
                             //int /* nGidEnt */, int nLidEnt, int nObj,
                             //ZOLTAN_ID_PTR /* gids */, ZOLTAN_ID_PTR lids,
                             //int *parts,
                             int level,
                             struct Zoltan_Struct *zz,
                             int *ierr)
{

  // Let the application specify any balancing params for this level


  typedef typename Adapter::scalar_t    pcoord_t;
  typedef typename Adapter::part_t      part_t;

  const MachineRepresentation<pcoord_t, part_t> *machine
    = static_cast<MachineRepresentation<pcoord_t, part_t> *>(data);
  *ierr = ZOLTAN_OK;

//  std::cout << "\nRank: " << machine->getMyRank() << " HierMethod 4, machine: " << machine->getMachineDim() << std::endl;

//  int ranks = machine->getNumRanks();

  int rank = machine->getMyRank();

//  std::cout << "\n--\nRank: " << rank
//    << " zoltanHierMethod, level: " << level
//    << "\n--\n"<< std::endl;

  Zoltan_Set_Param(zz, "LB_Approach", "repartition");
//  Zoltan_Set_Param(zz, "LB_Method", "RCB");
  Zoltan_Set_Param(zz, "LB_Method", "PARMETIS");

  Zoltan_Set_Param(zz, "CHECK_GRAPH", "0");
  Zoltan_Set_Param(zz, "GRAPH_BUILD_TYPE", "FAST_NO_DUP");

//  Zoltan_Set_Param(zz, "TFLOPS_SPECIAL", "1");
//  Zoltan_Set_Param(zz, "", "1");

  part_t num_unique_groups = machine->getNumUniqueGroups();
  std::vector<part_t> group_count;

  if (level == 0) {
    zz->Num_Unique_Groups = num_unique_groups;
    machine->getGroupCount2(group_count);
  }
  else if (level == 1) {
    std::vector<std::vector<part_t>> subgroup_counts;
    std::vector<part_t> num_unique_subgroups;

    machine->getGroupCount2(group_count);
    machine->getSubgroupCounts(subgroup_counts);
    machine->getNumUniqueSubgroups(num_unique_subgroups);

//    pcoord_t * xyz = new pcoord_t[3];

//    machine->getMyMachineCoordinate(xyz);
//    int group_idx = xyz[0];



//    std::cout << "\nRank: " << rank << " NUG: " << num_unique_groups << ",   ";
//    for (int i = 0; i < num_unique_groups; ++i) {
//        std::cout << " " << group_count[i];
//    }
//    std::cout << std::endl;


    int group_idx = 0;
    int upper_cdf = 0;
    int lower_cdf = 0;

    for (int i = 0; i < num_unique_groups; ++i) {

      lower_cdf = upper_cdf;
      upper_cdf += int(group_count[i]);

      if (rank < upper_cdf && rank >= lower_cdf) {
        group_idx = i;
        break;
      }
    }


//    std::cout << "Rank: " << rank << " mach coord: "
//      << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;

//    std::cout << "\n\nSGC Rank: " << rank << std::endl;
//    for (size_t i = 0; i < subgroup_counts.size(); ++i) {

//      std::cout << "Rank: " << rank << " Group: " << i << "  ";
//      for (size_t j = 0; j < subgroup_counts[i].size(); ++j) {
//        std::cout << " " << subgroup_counts[i][j];
//      }
//      std::cout << std::endl;
//    }

//    std::cout << "\n\nRank: " << rank
//      << " Level 1 Group idx: " << group_idx << std::endl;
    //  << " nug: " << zz->Num_Unique_Groups << std::endl;
   // std::cout << "\nRank: " << rank
   //   << " Level 1 Group idx: " << group_idx << std::endl;

//    std::cout << "\nRank: " << rank << " hi inside 0" << std::endl;
    zz->Num_Unique_Groups = int(num_unique_subgroups[group_idx]);
//    std::cout << "\nRank: " << rank << " hi inside 1" << std::endl;
    group_count.resize(zz->Num_Unique_Groups);
//    std::cout << "\nRank: " << rank << " hi inside 2" << std::endl;
    group_count = subgroup_counts[group_idx];
//    std::cout << "\nRank: " << rank << " hi inside 3" << std::endl;

//    std::cout << "\nRank: " << rank
//      << " Level 1 Group idx: " << group_idx
//      << " nug: " << zz->Num_Unique_Groups << std::endl;

    group_count.resize(zz->Num_Unique_Groups);

//    delete [] xyz;
  }
  // level > 1
  else {
    std::vector<std::vector<part_t>> subgroup_counts;
    std::vector<part_t> num_unique_subgroups;

    machine->getGroupCount2(group_count);
    machine->getSubgroupCounts(subgroup_counts);
    machine->getNumUniqueSubgroups(num_unique_subgroups);

//    pcoord_t * xyz = new pcoord_t[3];

//    if (rank == 0) {
//        std::cout << "\nNUG: " << num_unique_groups << ",   ";
//        for (int i = 0; i < num_unique_groups; ++i) {
//            std::cout << " " << group_count[i];
//        }
//        std::cout << std::endl;
//    }

//    machine->getMyMachineCoordinate(xyz);
//    int group_idx = int(xyz[0]);
//    int subgroup_idx = int(xyz[1]);
    int group_idx = 0;
    int upper_cdf = 0;
    int lower_cdf = 0;

    for (int i = 0; i < num_unique_groups; ++i) {

      lower_cdf = upper_cdf;
      upper_cdf += int(group_count[i]);

      if (rank < upper_cdf && rank >= lower_cdf) {

//      std::cout << "\nRank: " << rank << " solving part " << i << std::endl;

        group_idx = i;
        break;
      }
    }



    int subgroup_idx = 0;
    int group_rank = rank - lower_cdf;

    lower_cdf = 0;
    upper_cdf = 0;

    for (int i = 0; i < int(num_unique_subgroups[group_idx]); ++i) {

      lower_cdf = upper_cdf;
      upper_cdf += int(subgroup_counts[group_idx][i]);

      if (group_rank < upper_cdf && group_rank >= lower_cdf) {

//        std::cout << "\nPartMethod Rank: " << rank
//          << " solving part " << i
//          << " of Group " << group_idx
//          << " with group rank " << group_rank << std::endl;

        subgroup_idx = i;
        break;
      }
    }







//    std::cout << "\nRank: " << rank
//      << " Level 2+ Group idx: " << group_idx
//      << " Subgroup idx: " << subgroup_idx
//      << std::endl;

    zz->Num_Unique_Groups = int(subgroup_counts[group_idx][subgroup_idx]);

    group_count.resize(zz->Num_Unique_Groups);
//    std::iota(group_count.begin(), group_count.end(), 0);

    for (size_t i = 0; i < group_count.size(); ++i) {
        group_count[i] = i;
    }


//    std::cout << "\nRank: " << rank
//      << " Level 2+ Group idx: " << group_idx
//      << " Subgroup idx: " << subgroup_idx
//      << " nug " << zz->Num_Unique_Groups
//      << std::endl;

//    delete [] xyz;
  }
//  zz->Groups =


//  std::cout << "\nRank: " << rank
//    << " group_count size!: "
//    << group_count.size() << std::endl;

  std::vector<int> group_count_ints(group_count.begin(), group_count.end());

//  zz->Groups = group_count_ints.data();

//  zz->Group_Count = new int[zz->Num_Unique_Groups];
  zz->Group_Count =
//    (int *) std::malloc(zz->Num_Unique_Groups * sizeof(int));
    (int *) ZOLTAN_MALLOC(zz->Num_Unique_Groups * sizeof(int));


//  if (machine->getMyRank() == 3)
//    group_count_ints[5] = 777;

  std::copy(group_count_ints.begin(),
            group_count_ints.end(),
            zz->Group_Count);


//  std::cout << "\nRank: " << machine->getMyRank()
//   << " Callback groups address: " << (void *) zz->Group_Count
//   << std::endl;


//  zz->Current_Hier_Level = level;

//  char msg[16];

//  sprintf(msg, "%d", 1);
//  Zoltan_Set_Param(zz, "Num_Local_Parts", msg);

//  sprintf(msg, "%d", ranks);
//  Zoltan_Set_Param(zz, "Num_Global_Parts", msg);

//  zz->LB.Approach = "repartition"
//  zz->LB.Method = "RCB"
//  zz->Imb_Tol_Len = machine->getNumRanks();

//  Zoltan_Set_Param(zz, "Imb_Tol_Len", machine->getNumRanks());

//  ZOLTAN_FREE(&msg);

  return;
}

}


#endif
