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

/*! \file Zoltan2_GraphMetricValuesUtility.hpp
 */

#ifndef ZOLTAN2_GRAPHICMETRICVALUESUTILITY_HPP
#define ZOLTAN2_GRAPHICMETRICVALUESUTILITY_HPP

#include <Zoltan2_ImbalanceMetrics.hpp>
#include <Zoltan2_MetricUtility.hpp>
#include <zoltan_dd.h>
#include <Zoltan2_TPLTraits.hpp>
#include <map>
#include <vector>

namespace Zoltan2{


template <typename Adapter, typename MachineRep>
void globalWeightedCutsMessagesHopsByPart(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graph,
    const ArrayView<const typename Adapter::part_t> &parts,
    typename Adapter::part_t &numParts,
    ArrayRCP<RCP<BaseClassMetrics<typename Adapter::scalar_t> > > &metrics,
    ArrayRCP<typename Adapter::scalar_t> &globalSums,
    const RCP <const MachineRep> machine)
{
  env->debug(DETAILED_STATUS, "Entering globalWeightedCutsMessagesHopsByPart");
  //////////////////////////////////////////////////////////
  // Initialize return values

  typedef typename Adapter::lno_t t_lno_t;
  typedef typename Adapter::gno_t t_gno_t;
  typedef typename Adapter::offset_t t_offset_t;
  typedef typename Adapter::scalar_t t_scalar_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::node_t t_node_t;


  typedef typename Zoltan2::GraphModel<typename Adapter::base_adapter_t>::input_t t_input_t;

  t_lno_t localNumVertices = graph->getLocalNumVertices();
  t_lno_t localNumEdges = graph->getLocalNumEdges();

  ArrayView<const t_gno_t> Ids;
  ArrayView<t_input_t> v_wghts;
  graph->getVertexList(Ids, v_wghts);

  typedef GraphMetrics<t_scalar_t> gm_t;

  //get the edge ids, and weights
  ArrayView<const t_gno_t> edgeIds;
  ArrayView<const t_offset_t> offsets;
  ArrayView<t_input_t> e_wgts;
  graph->getEdgeList(edgeIds, offsets, e_wgts);


  std::vector <t_scalar_t> edge_weights;
  int numWeightPerEdge = graph->getNumWeightsPerEdge();

  int numMetrics = 4;                   // "edge cuts", messages, hops, weighted hops
  if (numWeightPerEdge) numMetrics += numWeightPerEdge * 2;   // "weight n", weighted hops per weight n

  // add some more metrics to the array
  auto next = metrics.size(); // where we begin filling
  for (auto n = 0; n < numMetrics; ++n)  {
    addNewMetric<gm_t, t_scalar_t>(env, metrics);
  }

  std::vector <part_t> e_parts (localNumEdges);
#ifdef HAVE_ZOLTAN2_MPI
  if (comm->getSize() > 1)
  {
    Zoltan_DD_Struct *dd = NULL;

    MPI_Comm mpicomm = Teuchos::getRawMpiComm(*comm);
    int size_gnot = Zoltan2::TPL_Traits<ZOLTAN_ID_PTR, t_gno_t>::NUM_ID;

    int debug_level = 0;
    Zoltan_DD_Create(&dd, mpicomm,
        size_gnot, 0,
        sizeof(part_t), localNumVertices, debug_level);

    ZOLTAN_ID_PTR ddnotneeded = NULL;  // Local IDs not needed
    Zoltan_DD_Update(
        dd,
        (localNumVertices ? (ZOLTAN_ID_PTR) Ids.getRawPtr() : NULL),
        ddnotneeded,
        (localNumVertices ? (char *) &(parts[0]) : NULL),
        NULL,
        int(localNumVertices));

    Zoltan_DD_Find(
        dd,
        (localNumEdges ? (ZOLTAN_ID_PTR) edgeIds.getRawPtr() : NULL),
        ddnotneeded,
        (localNumEdges ? (char *)&(e_parts[0]) : NULL),
        NULL,
        localNumEdges,
        NULL
        );
    Zoltan_DD_Destroy(&dd);
  } else
#endif
  {

    std::map<t_gno_t,t_lno_t> global_id_to_local_index;

    //else everything is local.
    //we need a globalid to local index conversion.
    //this does not exists till this point, so we need to create one.
    for (t_lno_t i = 0; i < localNumVertices; ++i){
      //at the local index i, we have the global index Ids[i].
      //so write i, to Ids[i] index of the vector.
      global_id_to_local_index[Ids[i]] = i;
    }

    for (t_lno_t i = 0; i < localNumEdges; ++i){
      t_gno_t ei = edgeIds[i];
      //ei is the global index of the neighbor one.
      part_t p = parts[global_id_to_local_index[ei]];
      e_parts[i] = p;
    }
  }

  RCP<const Teuchos::Comm<int> > tcomm = comm;

  env->timerStart(MACRO_TIMERS, "Communication Graph Create");
  {
    //get the vertices in each part in my part.
    std::vector <t_lno_t> part_begins(numParts, -1);
    std::vector <t_lno_t> part_nexts(localNumVertices, -1);

    //cluster vertices according to their parts.
    //create local part graph.
    for (t_lno_t i = 0; i < localNumVertices; ++i){
      part_t ap = parts[i];
      part_nexts[i] = part_begins[ap];
      part_begins[ap] = i;
    }


    for (int weight_index = -1; weight_index < numWeightPerEdge ; ++weight_index){

      //MD: these two should be part_t.
      //but we dont want to compile tpetra from the beginning.
      //This can be changed when directory is updated.
      typedef t_lno_t local_part_type;
      typedef t_gno_t global_part_type;

      typedef Tpetra::Map<local_part_type, global_part_type, t_node_t> map_t;
      Teuchos::RCP<const map_t> map = Teuchos::rcp (new map_t (numParts, 0, tcomm));

      typedef Tpetra::CrsMatrix<t_scalar_t, local_part_type, global_part_type, t_node_t> tcrsMatrix_t;
      Teuchos::RCP<tcrsMatrix_t> tMatrix(new tcrsMatrix_t (map, 0));


      std::vector <global_part_type> part_neighbors (numParts);

      std::vector <t_scalar_t> part_neighbor_weights(numParts, 0);
      std::vector <t_scalar_t> part_neighbor_weights_ordered(numParts);

      //coarsen for all vertices in my part in order with parts.
      for (global_part_type i = 0; i < (global_part_type) numParts; ++i){
        part_t num_neighbor_parts = 0;
        t_lno_t v = part_begins[i];
        //get part i, and first vertex in this part v.
        while (v != -1){
          //now get the neightbors of v.
          for (t_offset_t j = offsets[v]; j < offsets[v+1]; ++j){
            //get the part of the second vertex.
            part_t ep = e_parts[j];

            t_scalar_t ew = 1;
            if (weight_index > -1){
              ew = e_wgts[weight_index][j];
            }
            //add it to my local part neighbors for part i.
            if (part_neighbor_weights[ep] < 0.00001){
              part_neighbors[num_neighbor_parts++] = ep;
            }
            part_neighbor_weights[ep] += ew;
          }
          v = part_nexts[v];
        }

        //now get the part list.
        for (t_lno_t j = 0; j < num_neighbor_parts; ++j){
          part_t neighbor_part = part_neighbors[j];
          part_neighbor_weights_ordered[j] = part_neighbor_weights[neighbor_part];
          part_neighbor_weights[neighbor_part] = 0;
        }

        //insert it to tpetra crsmatrix.
        if (num_neighbor_parts > 0){
          Teuchos::ArrayView<const global_part_type> destinations(&(part_neighbors[0]), num_neighbor_parts);
          Teuchos::ArrayView<const t_scalar_t> vals(&(part_neighbor_weights_ordered[0]), num_neighbor_parts);
          tMatrix->insertGlobalValues (i,destinations, vals);
        }
      }
      tMatrix->fillComplete ();
      local_part_type num_local_parts = map->getNodeNumElements();

      Array<global_part_type> Indices;
      Array<t_scalar_t> Values;

      t_scalar_t max_edge_cut = 0;
      t_scalar_t total_edge_cut = 0;
      global_part_type max_message = 0;
      global_part_type total_message = 0;

      global_part_type total_hop_count = 0;
      t_scalar_t total_weighted_hop_count = 0;
      global_part_type max_hop_count = 0;
      t_scalar_t max_weighted_hop_count = 0;

      for (local_part_type i=0; i < num_local_parts; i++) {

        const global_part_type globalRow = map->getGlobalElement(i);
        size_t NumEntries = tMatrix->getNumEntriesInGlobalRow (globalRow);
        Indices.resize (NumEntries);
        Values.resize (NumEntries);
        tMatrix->getGlobalRowCopy (globalRow,Indices(),Values(),NumEntries);

        t_scalar_t part_edge_cut = 0;
        global_part_type part_messages = 0;

        for (size_t j=0; j < NumEntries; j++){
          if (Indices[j] != globalRow){
            part_edge_cut += Values[j];
            part_messages += 1;

            typename MachineRep::machine_pcoord_t hop_count = 0;
            machine->getHopCount(globalRow, Indices[j], hop_count);

            global_part_type hop_counts = hop_count;
            t_scalar_t weighted_hop_counts = hop_count * Values[j];

            total_hop_count += hop_counts;
            total_weighted_hop_count += weighted_hop_counts;

            if (hop_counts > max_hop_count ){
              max_hop_count = hop_counts;
            }
            if (weighted_hop_counts > max_weighted_hop_count ){
              max_weighted_hop_count = weighted_hop_counts;
            }
          }
        }
        if (part_edge_cut > max_edge_cut){
          max_edge_cut = part_edge_cut;
        }
        total_edge_cut += part_edge_cut;

        if (part_messages > max_message){
          max_message = part_messages;
        }
        total_message += part_messages;

      }
      t_scalar_t g_max_edge_cut = 0;
      t_scalar_t g_total_edge_cut = 0;
      global_part_type g_max_message = 0;
      global_part_type g_total_message = 0;



      global_part_type g_total_hop_count = 0;
      t_scalar_t g_total_weighted_hop_count = 0;
      global_part_type g_max_hop_count = 0;
      t_scalar_t g_max_weighted_hop_count = 0;

      try{

        Teuchos::reduceAll<int, t_scalar_t>(*comm,Teuchos::REDUCE_MAX,1,&max_edge_cut,&g_max_edge_cut);
        Teuchos::reduceAll<int, global_part_type>(*comm,Teuchos::REDUCE_MAX,1,&max_message,&g_max_message);

        Teuchos::reduceAll<int, global_part_type>(*comm,Teuchos::REDUCE_MAX,1,&max_hop_count,&g_max_hop_count);
        Teuchos::reduceAll<int, t_scalar_t>(*comm,Teuchos::REDUCE_MAX,1,&max_weighted_hop_count,&g_max_weighted_hop_count);

        Teuchos::reduceAll<int, t_scalar_t>(*comm,Teuchos::REDUCE_SUM,1,&total_edge_cut,&g_total_edge_cut);
        Teuchos::reduceAll<int, global_part_type>(*comm,Teuchos::REDUCE_SUM,1,&total_message,&g_total_message);

        Teuchos::reduceAll<int, global_part_type>(*comm,Teuchos::REDUCE_SUM,1,&total_hop_count,&g_total_hop_count);
        Teuchos::reduceAll<int, t_scalar_t>(*comm,Teuchos::REDUCE_SUM,1,&total_weighted_hop_count,&g_total_weighted_hop_count);

      }
      Z2_THROW_OUTSIDE_ERROR(*env);


      if (weight_index == -1){
        metrics[next]->setName("md edge cuts");
      }
      else {
        std::ostringstream oss;
        oss << "md weighted edge cuts" << weight_index;
        metrics[next]->setName( oss.str());
      }

      metrics[next]->setMetricValue("global maximum", g_max_edge_cut);
      metrics[next]->setMetricValue("global sum", g_total_edge_cut);
      next++;

      if (weight_index == -1){
        metrics[next]->setName("message");
        metrics[next]->setMetricValue("global maximum", g_max_message);
        metrics[next]->setMetricValue("global sum", g_total_message);
        next++;
      }


      if (weight_index == -1){
        metrics[next]->setName("hops (No Weight)");
        metrics[next]->setMetricValue("global maximum", g_max_hop_count);
        metrics[next]->setMetricValue("global sum", g_total_hop_count);
        next++;
      }

      std::ostringstream oss;
      oss << "weighted hops" << weight_index;
      metrics[next]->setName( oss.str());
      metrics[next]->setMetricValue("global maximum", g_max_weighted_hop_count);
      metrics[next]->setMetricValue("global sum", g_total_weighted_hop_count);
      next++;

    }
  }
  env->timerStop(MACRO_TIMERS, "Communication Graph Create");

  env->debug(DETAILED_STATUS, "Exiting globalWeightedCutsMessagesHopsByPart");
}


template <typename Adapter>
void globalWeightedCutsMessagesByPart(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graph,
    const ArrayView<const typename Adapter::part_t> &parts,
    typename Adapter::part_t &numParts,
    ArrayRCP<RCP<BaseClassMetrics<typename Adapter::scalar_t> > > &metrics,
    ArrayRCP<typename Adapter::scalar_t> &globalSums)
{
  env->debug(DETAILED_STATUS, "Entering globalWeightedCutsMessagesByPart");
  //////////////////////////////////////////////////////////
  // Initialize return values

  typedef typename Adapter::lno_t t_lno_t;
  typedef typename Adapter::gno_t t_gno_t;
  typedef typename Adapter::offset_t t_offset_t;
  typedef typename Adapter::scalar_t t_scalar_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::node_t t_node_t;


  typedef typename Zoltan2::GraphModel<typename Adapter::base_adapter_t>::input_t t_input_t;

  t_lno_t localNumVertices = graph->getLocalNumVertices();
  t_lno_t localNumEdges = graph->getLocalNumEdges();

  ArrayView<const t_gno_t> Ids;
  ArrayView<t_input_t> v_wghts;
  graph->getVertexList(Ids, v_wghts);

  typedef GraphMetrics<t_scalar_t> gm_t;

  //get the edge ids, and weights
  ArrayView<const t_gno_t> edgeIds;
  ArrayView<const t_offset_t> offsets;
  ArrayView<t_input_t> e_wgts;
  graph->getEdgeList(edgeIds, offsets, e_wgts);


  std::vector <t_scalar_t> edge_weights;
  int numWeightPerEdge = graph->getNumWeightsPerEdge();

  int numMetrics = 2;                   // "edge cuts", messages
  if (numWeightPerEdge) numMetrics += numWeightPerEdge;   // "weight n"

  // add some more metrics to the array
  auto next = metrics.size(); // where we begin filling
  for (auto n = 0; n < numMetrics; ++n)  {
    addNewMetric<gm_t, t_scalar_t>(env, metrics);
  }

  std::vector <part_t> e_parts (localNumEdges);
#ifdef HAVE_ZOLTAN2_MPI
  if (comm->getSize() > 1)
  {
    Zoltan_DD_Struct *dd = NULL;

    MPI_Comm mpicomm = Teuchos::getRawMpiComm(*comm);
    int size_gnot = Zoltan2::TPL_Traits<ZOLTAN_ID_PTR, t_gno_t>::NUM_ID;

    int debug_level = 0;
    Zoltan_DD_Create(&dd, mpicomm,
        size_gnot, 0,
        sizeof(part_t), localNumVertices, debug_level);

    ZOLTAN_ID_PTR ddnotneeded = NULL;  // Local IDs not needed
    Zoltan_DD_Update(
        dd,
        (localNumVertices ? (ZOLTAN_ID_PTR) Ids.getRawPtr() : NULL),
        ddnotneeded,
        (localNumVertices ? (char *) &(parts[0]) : NULL),
        NULL,
        int(localNumVertices));

    Zoltan_DD_Find(
        dd,
        (localNumEdges ? (ZOLTAN_ID_PTR) edgeIds.getRawPtr() : NULL),
        ddnotneeded,
        (localNumEdges ? (char *)&(e_parts[0]) : NULL),
        NULL,
        localNumEdges,
        NULL
        );
    Zoltan_DD_Destroy(&dd);
  } else
#endif
  {

    std::map<t_gno_t,t_lno_t> global_id_to_local_index;

    //else everything is local.
    //we need a globalid to local index conversion.
    //this does not exists till this point, so we need to create one.
    for (t_lno_t i = 0; i < localNumVertices; ++i){
      //at the local index i, we have the global index Ids[i].
      //so write i, to Ids[i] index of the vector.
      global_id_to_local_index[Ids[i]] = i;
    }

    for (t_lno_t i = 0; i < localNumEdges; ++i){
      t_gno_t ei = edgeIds[i];
      //ei is the global index of the neighbor one.
      part_t p = parts[global_id_to_local_index[ei]];
      e_parts[i] = p;
    }
  }

  RCP<const Teuchos::Comm<int> > tcomm = comm;

  env->timerStart(MACRO_TIMERS, "Communication Graph Create");
  {
    //get the vertices in each part in my part.
    std::vector <t_lno_t> part_begins(numParts, -1);
    std::vector <t_lno_t> part_nexts(localNumVertices, -1);

    //cluster vertices according to their parts.
    //create local part graph.
    for (t_lno_t i = 0; i < localNumVertices; ++i){
      part_t ap = parts[i];
      part_nexts[i] = part_begins[ap];
      part_begins[ap] = i;
    }

    for (int weight_index = -1; weight_index < numWeightPerEdge ; ++weight_index){

      //MD: these two should be part_t.
      //but we dont want to compile tpetra from the beginning.
      //This can be changed when directory is updated.
      typedef t_lno_t local_part_type;
      typedef t_gno_t global_part_type;

      typedef Tpetra::Map<local_part_type, global_part_type, t_node_t> map_t;
      Teuchos::RCP<const map_t> map = Teuchos::rcp (new map_t (numParts, 0, tcomm));

      typedef Tpetra::CrsMatrix<t_scalar_t, local_part_type, global_part_type, t_node_t> tcrsMatrix_t;
      Teuchos::RCP<tcrsMatrix_t> tMatrix(new tcrsMatrix_t (map, 0));


      std::vector <global_part_type> part_neighbors (numParts);

      std::vector <t_scalar_t> part_neighbor_weights(numParts, 0);
      std::vector <t_scalar_t> part_neighbor_weights_ordered(numParts);

      //coarsen for all vertices in my part in order with parts.
      for (global_part_type i = 0; i < (global_part_type) numParts; ++i){
        part_t num_neighbor_parts = 0;
        t_lno_t v = part_begins[i];
        //get part i, and first vertex in this part v.
        while (v != -1){
          //now get the neightbors of v.
          for (t_offset_t j = offsets[v]; j < offsets[v+1]; ++j){
            //get the part of the second vertex.
            part_t ep = e_parts[j];

            t_scalar_t ew = 1;
            if (weight_index > -1){
              ew = e_wgts[weight_index][j];
            }
            //add it to my local part neighbors for part i.
            if (part_neighbor_weights[ep] < 0.00001){
              part_neighbors[num_neighbor_parts++] = ep;
            }
            part_neighbor_weights[ep] += ew;
          }
          v = part_nexts[v];
        }

        //now get the part list.
        for (t_lno_t j = 0; j < num_neighbor_parts; ++j){
          part_t neighbor_part = part_neighbors[j];
          part_neighbor_weights_ordered[j] = part_neighbor_weights[neighbor_part];
          part_neighbor_weights[neighbor_part] = 0;
        }

        //insert it to tpetra crsmatrix.
        if (num_neighbor_parts > 0){
          Teuchos::ArrayView<const global_part_type> destinations(&(part_neighbors[0]), num_neighbor_parts);
          Teuchos::ArrayView<const t_scalar_t> vals(&(part_neighbor_weights_ordered[0]), num_neighbor_parts);
          tMatrix->insertGlobalValues (i,destinations, vals);
        }
      }
      tMatrix->fillComplete ();
      local_part_type num_local_parts = map->getNodeNumElements();

      Array<global_part_type> Indices;
      Array<t_scalar_t> Values;

      t_scalar_t max_edge_cut = 0;
      t_scalar_t total_edge_cut = 0;
      global_part_type max_message = 0;
      global_part_type total_message = 0;

      for (local_part_type i=0; i < num_local_parts; i++) {

        const global_part_type globalRow = map->getGlobalElement(i);
        size_t NumEntries = tMatrix->getNumEntriesInGlobalRow (globalRow);
        Indices.resize (NumEntries);
        Values.resize (NumEntries);
        tMatrix->getGlobalRowCopy (globalRow,Indices(),Values(),NumEntries);

        t_scalar_t part_edge_cut = 0;
        global_part_type part_messages = 0;
        for (size_t j=0; j < NumEntries; j++){
          if (Indices[j] != globalRow){
            part_edge_cut += Values[j];
            part_messages += 1;
          }
        }
        if (part_edge_cut > max_edge_cut){
          max_edge_cut = part_edge_cut;
        }
        total_edge_cut += part_edge_cut;

        if (part_messages > max_message){
          max_message = part_messages;
        }
        total_message += part_messages;

      }
      t_scalar_t g_max_edge_cut = 0;
      t_scalar_t g_total_edge_cut = 0;
      global_part_type g_max_message = 0;
      global_part_type g_total_message = 0;

      try{

        Teuchos::reduceAll<int, t_scalar_t>(*comm,Teuchos::REDUCE_MAX,1,&max_edge_cut,&g_max_edge_cut);
        Teuchos::reduceAll<int, global_part_type>(*comm,Teuchos::REDUCE_MAX,1,&max_message,&g_max_message);

        Teuchos::reduceAll<int, t_scalar_t>(*comm,Teuchos::REDUCE_SUM,1,&total_edge_cut,&g_total_edge_cut);
        Teuchos::reduceAll<int, global_part_type>(*comm,Teuchos::REDUCE_SUM,1,&total_message,&g_total_message);
      }
      Z2_THROW_OUTSIDE_ERROR(*env);

      if (weight_index == -1){
        metrics[next]->setName("md edge cuts");
      }
      else {
        std::ostringstream oss;
        oss << "md weight " << weight_index;
        metrics[next]->setName( oss.str());
      }
      metrics[next]->setMetricValue("global maximum", g_max_edge_cut);
      metrics[next]->setMetricValue("global sum", g_total_edge_cut);
      next++;
      if (weight_index == -1){
        metrics[next]->setName("md message");
        metrics[next]->setMetricValue("global maximum", g_max_message);
        metrics[next]->setMetricValue("global sum", g_total_message);
        next++;
      }
    }
  }
  env->timerStop(MACRO_TIMERS, "Communication Graph Create");

  env->debug(DETAILED_STATUS, "Exiting globalWeightedCutsMessagesByPart");
}

/*! \brief Given the local partitioning, compute the global weighted cuts in each part.
 *
 *   \param env   Environment for error handling
 *   \param comm   communicator
 *   \param graph Graph model
 *   \param part   \c part[i] is the part ID for local object \c i
 *   \param numParts  on return this is the global number of parts.
 *   \param metrics on return points to a list of named GraphMetricValues cuts
 *     that each contains the global max and sum over parts of
 *     the item being measured. The list may contain "edge cuts", or
 *     "weight 0", "weight 1" and so on in that order.
 *     If uniform weights were given, then only "edge cuts" appears.
 *     If one set of non-uniform weights were given, then
 *     "weight 0" appear.  Finally, if multiple
 *     weights were given, we have
 *     the individual weights "weight 0", "weight 1", and so on.
 *   \param globalSums If weights are uniform, the globalSums is the
 *      \c numParts totals of global number of cuts in each part.
 *     Suppose the number of weights is \c W.  If
 *     W is 1, then on return this is an array of length \c numParts .
 *     The \c numParts entries are the total weight in each part.
 *     If \c W is greater than one, then the length of this array is
 *     \c W*numParts .
 *     The entries are the sum of the individual weights in each part,
 *     by weight index by part number.  The array is allocated here.
 *
 * globalWeightedCutsByPart() must be called by all processes in \c comm.
 */
template <typename Adapter>
  void globalWeightedCutsByPart(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graph,
    const ArrayView<const typename Adapter::part_t> &part,
    typename Adapter::part_t &numParts,
    ArrayRCP<RCP<BaseClassMetrics<typename Adapter::scalar_t> > > &metrics,
    ArrayRCP<typename Adapter::scalar_t> &globalSums)
{
  env->debug(DETAILED_STATUS, "Entering globalWeightedCutsByPart");
  //////////////////////////////////////////////////////////
  // Initialize return values

  numParts = 0;

  int ewgtDim = graph->getNumWeightsPerEdge();

  int numMetrics = 1;                   // "edge cuts"
  if (ewgtDim) numMetrics += ewgtDim;   // "weight n"

  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::node_t node_t;
  typedef typename Adapter::part_t part_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  typedef Tpetra::CrsMatrix<part_t,lno_t,gno_t,node_t>  sparse_matrix_type;
  typedef Tpetra::Vector<part_t,lno_t,gno_t,node_t>     vector_t;
  typedef Tpetra::Map<lno_t, gno_t, node_t>             map_type;
  typedef Tpetra::global_size_t GST;
  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

  using Teuchos::as;

  auto next = metrics.size(); // where we begin filling
  typedef GraphMetrics<scalar_t> gm_t;
  for (auto n = 0; n < numMetrics; ++n)  {
    addNewMetric<gm_t, scalar_t>(env, metrics);
  }

  //////////////////////////////////////////////////////////
  // Figure out the global number of parts in use.
  // Verify number of vertex weights is the same everywhere.

  lno_t localNumObj = part.size();
  part_t localNum[2], globalNum[2];
  localNum[0] = static_cast<part_t>(ewgtDim);
  localNum[1] = 0;

  for (lno_t i=0; i < localNumObj; i++)
    if (part[i] > localNum[1]) localNum[1] = part[i];

  try{
    reduceAll<int, part_t>(*comm, Teuchos::REDUCE_MAX, 2,
      localNum, globalNum);
  }
  Z2_THROW_OUTSIDE_ERROR(*env)

  env->globalBugAssertion(__FILE__,__LINE__,
    "inconsistent number of edge weights",
    globalNum[0] == localNum[0], DEBUG_MODE_ASSERTION, comm);

  part_t nparts = globalNum[1] + 1;

  part_t globalSumSize = nparts * numMetrics;
  scalar_t * sumBuf = new scalar_t [globalSumSize];
  env->localMemoryAssertion(__FILE__, __LINE__, globalSumSize, sumBuf);
  globalSums = arcp(sumBuf, 0, globalSumSize);

  //////////////////////////////////////////////////////////
  // Calculate the local totals by part.

  scalar_t *localBuf = new scalar_t [globalSumSize];
  env->localMemoryAssertion(__FILE__,__LINE__,globalSumSize,localBuf);
  memset(localBuf, 0, sizeof(scalar_t) * globalSumSize);

  scalar_t *cut = localBuf;              // # of cuts

  ArrayView<const gno_t> Ids;
  ArrayView<input_t> vwgts;
  //size_t nv =
  graph->getVertexList(Ids, vwgts);

  ArrayView<const gno_t> edgeIds;
  ArrayView<const offset_t> offsets;
  ArrayView<input_t> wgts;
  //size_t numLocalEdges =
  graph->getEdgeList(edgeIds, offsets, wgts);
  // **************************************************************************
  // *************************** BUILD MAP FOR ADJS ***************************
  // **************************************************************************

  RCP<const map_type> vertexMapG;

  // Build a list of the global vertex ids...
  gno_t min = std::numeric_limits<gno_t>::max();
  offset_t maxcols = 0;
  for (lno_t i = 0; i < localNumObj; ++i) {
    if (Ids[i] < min) min = Ids[i];
    offset_t ncols = offsets[i+1] - offsets[i];
    if (ncols > maxcols) maxcols = ncols;
  }

  gno_t gmin;
  Teuchos::reduceAll<int, gno_t>(*comm,Teuchos::REDUCE_MIN,1,&min,&gmin);

  //Generate Map for vertex
  vertexMapG = rcp(new map_type(INVALID, Ids, gmin, comm));

  // **************************************************************************
  // ************************** BUILD GRAPH FOR ADJS **************************
  // **************************************************************************

  //MD:Zoltan Directory could be used instead of adjMatrix.

  RCP<sparse_matrix_type> adjsMatrix;

  // Construct Tpetra::CrsGraph objects.
  adjsMatrix = rcp (new sparse_matrix_type (vertexMapG, 0));

  Array<part_t> justOneA(maxcols, 1);

  for (lno_t localElement=0; localElement<localNumObj; ++localElement){
    // Insert all columns for global row Ids[localElement]
    offset_t ncols = offsets[localElement+1] - offsets[localElement];
    adjsMatrix->insertGlobalValues(Ids[localElement],
                                   edgeIds(offsets[localElement], ncols),
                                   justOneA(0, ncols));
  }

  //Fill-complete adjs Graph
  adjsMatrix->fillComplete ();

  // Compute part
  RCP<vector_t> scaleVec = Teuchos::rcp( new vector_t(vertexMapG,false) );
  for (lno_t localElement=0; localElement<localNumObj; ++localElement) {
    scaleVec->replaceLocalValue(localElement,part[localElement]);
  }

  // Postmultiply adjsMatrix by part
  adjsMatrix->rightScale(*scaleVec);
  Array<gno_t> Indices;
  Array<part_t> Values;

  for (lno_t i=0; i < localNumObj; i++) {
    const gno_t globalRow = Ids[i];
    size_t NumEntries = adjsMatrix->getNumEntriesInGlobalRow (globalRow);
    Indices.resize (NumEntries);
    Values.resize (NumEntries);
    adjsMatrix->getGlobalRowCopy (globalRow,Indices(),Values(),NumEntries);

    for (size_t j=0; j < NumEntries; j++)
      if (part[i] != Values[j])
	cut[part[i]]++;
  }

  if (numMetrics > 1) {

    scalar_t *wgt = localBuf + nparts; // weight 0

    // This code assumes the solution has the part ordered the
    // same way as the user input.  (Bug 5891 is resolved.)
    for (int edim = 0; edim < ewgtDim; edim++){
      for (lno_t i=0; i < localNumObj; i++) {
        const gno_t globalRow = Ids[i];
        size_t NumEntries = adjsMatrix->getNumEntriesInGlobalRow (globalRow);
        Indices.resize (NumEntries);
        Values.resize (NumEntries);
        adjsMatrix->getGlobalRowCopy (globalRow,Indices(),Values(),NumEntries);

        for (size_t j=0; j < NumEntries; j++)
          if (part[i] != Values[j])
            wgt[part[i]] += wgts[edim][offsets[i] + j];
      }
      wgt += nparts;         // individual weights
    }
  }

  //////////////////////////////////////////////////////////
  // Obtain global totals by part.

  try{
    reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_SUM, globalSumSize,
      localBuf, sumBuf);
  }
  Z2_THROW_OUTSIDE_ERROR(*env);

  delete [] localBuf;

  //////////////////////////////////////////////////////////
  // Global max and sum over all parts

  cut = sumBuf;                     // # of cuts
  scalar_t max=0, sum=0;

  ArrayView<scalar_t> cutVec(cut, nparts);
  getStridedStats<scalar_t>(cutVec, 1, 0, max, sum);

  metrics[next]->setName("edge cuts");
  metrics[next]->setMetricValue("global maximum", max);
  metrics[next]->setMetricValue("global sum", sum);

  next++;

  if (numMetrics > 1){
    scalar_t *wgt = sumBuf + nparts;        // weight 0

    for (int edim=0; edim < ewgtDim; edim++){
      ArrayView<scalar_t> fromVec(wgt, nparts);
      getStridedStats<scalar_t>(fromVec, 1, 0, max, sum);

      std::ostringstream oss;
      oss << "weight " << edim;

      metrics[next]->setName(oss.str());
      metrics[next]->setMetricValue("global maximum", max);
      metrics[next]->setMetricValue("global sum", sum);

      next++;
      wgt += nparts;       // individual weights
    }
  }

  numParts = nparts;

  env->debug(DETAILED_STATUS, "Exiting globalWeightedCutsByPart");
}

/*! \brief Print out header info for graph metrics.
 */
template <typename scalar_t, typename part_t>
void printGraphMetricsHeader(std::ostream &os, part_t targetNumParts, part_t numParts )
{
  os << "Graph Metrics:  (" << numParts << " parts)";
  os << std::endl;
  if (targetNumParts != numParts) {
    os << "Target number of parts is: " << targetNumParts << std::endl;
  }
  GraphMetrics<scalar_t>::printHeader(os);
}

/*! \brief Print out list of graph metrics.
 */
template <typename scalar_t, typename part_t>
void printGraphMetrics(std::ostream &os, part_t targetNumParts, part_t numParts, const ArrayView<RCP<BaseClassMetrics<scalar_t> > > &infoList)
{
  printGraphMetricsHeader<scalar_t, part_t>(os, targetNumParts, numParts);
  for (int i=0; i < infoList.size(); i++) {
    if (infoList[i]->getName() != METRICS_UNSET_STRING) {
      infoList[i]->printLine(os);
    }
  }
  os << std::endl;
}

/*! \brief Print out header and a single graph metric.
 */
template <typename scalar_t, typename part_t>
void printGraphMetrics(std::ostream &os, part_t targetNumParts, part_t numParts, RCP<BaseClassMetrics<scalar_t>> metricValue)
{
  printGraphMetricsHeader<scalar_t, part_t>(os, targetNumParts, numParts);
  metricValue->printLine(os);
}

} //namespace Zoltan2

#endif
