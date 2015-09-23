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

#ifdef HAVE_ZOLTAN2_HYPERGRAPHMODEL
#include <Zoltan2_HyperGraphModel.hpp>
#endif

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
  *ierr = ZOLTAN_OK;

  size_t mynObj = adp->getLocalNumIDs();
   
  const typename Adapter::gno_t *myids = NULL;
  adp->getIDsView(myids);
  for (size_t i = 0; i < mynObj; i++) {
    gids[i] = ZOLTAN_ID_TYPE(myids[i]); // TODO TRAITS CONVERSION MAY BE NEEDED
    lids[i] = ZOLTAN_ID_TYPE(i);      // TODO TRAITS CONVERSION MAY BE NEEDED
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
static void zoltanParts(void *data, int nGidEnt, int nLidEnt, int nObj,
                        ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                        int *parts, int *ierr)
{
  const Adapter *adp = static_cast<Adapter *>(data);
  *ierr = ZOLTAN_OK;
  const typename Adapter::part_t *myparts;
  adp->getPartsView(myparts);
  // User parts from input adapter
  for (int i = 0; i < nObj; i++)
    parts[i] = int(myparts[lids[i]]);
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
static void zoltanGeom(void *data, int nGidEnt, int nLidEnt, int nObj,
                       ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                       int nDim, double *coords, int *ierr)
{
  const Adapter *adp = static_cast<Adapter *>(data);
  *ierr = ZOLTAN_OK;

  for (int d = 0; d < nDim; d++) {
    const typename Adapter::scalar_t *mycoords;
    int mystride;
    adp->getCoordinatesView(mycoords, mystride, d);
    for (int i = 0; i < nObj; i++)
      coords[i*nDim+d] = double(mycoords[lids[i]*mystride]);
  }
}

/////////////////////////////////////////////////////////////////////////////
// MATRIX ADAPTER CALLBACKS
/////////////////////////////////////////////////////////////////////////////

///////////////////////
// ZOLTAN_HG_SIZE_CS_FN
template <typename Adapter>
static void zoltanHGSizeCSForMatrixAdapter(
  void *data, int *nEdges, int *nPins,
  int *format, int *ierr
) 
{
  std::cout << "HELLO FROM HGSizeCS with MATRIX ADAPTER" << std::endl;
  *ierr = ZOLTAN_FATAL;
}

//////////////////
// ZOLTAN_HG_CS_FN
template <typename Adapter>
static void zoltanHGCSForMatrixAdapter(
  void *data, int nGidEnt, int nEdges, int nPins,
  int format, ZOLTAN_ID_PTR edgeIds, 
  int *edgeIdx, ZOLTAN_ID_PTR pinIds, int *ierr
)
{
  std::cout << "HELLO FROM HGCS with MATRIX ADAPTER" << std::endl;
  *ierr = ZOLTAN_FATAL;
}

/////////////////////////////////////////////////////////////////////////////
// TODO:  GRAPH ADAPTER CALLBACKS


/////////////////////////////////////////////////////////////////////////////
// HYPERGRAPH MODEL CALLBACKS
/////////////////////////////////////////////////////////////////////////////

#ifdef HAVE_ZOLTAN2_HYPERGRAPHMODEL
////////////////////
// ZOLTAN_NUM_OBJ_FN
template <typename Adapter>
static int zoltanHGModelNumObj(void *data, int *ierr) {
  const HyperGraphModel<Adapter>* mdl = static_cast<HyperGraphModel<Adapter>* >(data);
  *ierr = ZOLTAN_OK;
  return int(mdl->getLocalNumOwnedVertices());
}

/////////////////////
// ZOLTAN_OBJ_LIST_FN
template <typename Adapter>
static void zoltanHGModelObjList(void *data, int nGidEnt, int nLidEnt, 
                          ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                          int wdim, float *wgts, int *ierr) 
{
  const HyperGraphModel<Adapter>* mdl = static_cast<HyperGraphModel<Adapter>* >(data);
  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::scalar_t    scalar_t;
  typedef StridedData<lno_t, scalar_t>  input_t;

  *ierr = ZOLTAN_OK;
  ArrayView<const gno_t> Ids;
  ArrayView<input_t> model_wgts;
  ArrayView<input_t> xyz;
  size_t num_verts = mdl->getVertexList(Ids,xyz,model_wgts);
  ArrayView<bool> isOwner;
  mdl->getOwnedList(isOwner);
  int j=0;
  for (size_t i=0;i<num_verts;i++) {
    if (isOwner[i]) {
      lids[j] = i;
      gids[j] = Ids[i];
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
static void zoltanHGModelSizeCSForMeshAdapter(
  void *data, int *nEdges, int *nPins,
  int *format, int *ierr
) 
{
  *ierr = ZOLTAN_OK;
  const HyperGraphModel<Adapter>* mdl = static_cast<HyperGraphModel<Adapter>* >(data);
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
static void zoltanHGModelCSForMeshAdapter(
  void *data, int nGidEnt, int nEdges, int nPins,
  int format, ZOLTAN_ID_PTR edgeIds, 
  int *edgeIdx, ZOLTAN_ID_PTR pinIds, int *ierr
)
{
  *ierr = ZOLTAN_OK;
  const HyperGraphModel<Adapter>* mdl = static_cast<HyperGraphModel<Adapter>* >(data);
  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::scalar_t    scalar_t;
  typedef StridedData<lno_t, scalar_t>  input_t;

  ArrayView<const gno_t> Ids;
  ArrayView<input_t> wgts;
  mdl->getEdgeList(Ids,wgts);
  ArrayView<const gno_t> pinIds_;
  ArrayView<const lno_t> offsets;
  ArrayView<input_t> pin_wgts;
  mdl->getPinList(pinIds_,offsets,pin_wgts);
  for (int i=0;i<nEdges;i++) {
    edgeIds[i]=Ids[i];
    edgeIdx[i]=offsets[i];
  }
  
  for (int i=0;i<nPins;i++)
    pinIds[i] = pinIds_[i];
}

/////////////////////////////////////////////////////////////////////////////
// MESH ADAPTER CALLBACKS
/////////////////////////////////////////////////////////////////////////////

///////////////////////
// ZOLTAN_HG_SIZE_CS_FN
template <typename Adapter>
static void zoltanHGSizeCSForMeshAdapter(
  void *data, int *nEdges, int *nPins,
  int *format, int *ierr
) 
{
  *ierr = ZOLTAN_OK;  
  typedef typename Adapter::user_t user_t;
  const MeshAdapter<user_t>* madp = static_cast<MeshAdapter<user_t>* >(data);
  *nEdges = madp->getLocalNumOf(madp->getAdjacencyEntityType());
  *nPins = madp->getLocalNumAdjs(madp->getAdjacencyEntityType(),madp->getPrimaryEntityType());
  *format = ZOLTAN_COMPRESSED_EDGE;
}

//////////////////
// ZOLTAN_HG_CS_FN
template <typename Adapter>
static void zoltanHGCSForMeshAdapter(
  void *data, int nGidEnt, int nEdges, int nPins,
  int format, ZOLTAN_ID_PTR edgeIds, 
  int *edgeIdx, ZOLTAN_ID_PTR pinIds, int *ierr
)
{
  *ierr = ZOLTAN_OK;
  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::lno_t       lno_t;  
  typedef typename Adapter::user_t      user_t;
  const MeshAdapter<user_t>* madp = static_cast<MeshAdapter<user_t>*>(data);
  const gno_t *Ids;
  madp->getIDsViewOf(madp->getAdjacencyEntityType(),Ids);
  const lno_t* offsets;
  const gno_t* adjIds;
  madp->getAdjsView(madp->getAdjacencyEntityType(), madp->getPrimaryEntityType(),
                    offsets, adjIds);
  for (int i=0;i<nEdges;i++) {
    edgeIds[i]=Ids[i];
    edgeIdx[i]=offsets[i];
  }
  for (int i=0;i<nPins;i++)
    pinIds[i] = adjIds[i];
}
#endif

}


#endif
