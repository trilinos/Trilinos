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

#include <Zoltan2_HyperGraphModel.hpp>

#include <Zoltan2_Util.hpp>
#include <Zoltan2_TPLTraits.hpp>
#include <zoltan_cpp.h>

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

}


#endif
