// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_GraphMetricValuesUtility.hpp
 */

#ifndef ZOLTAN2_GRAPHICMETRICVALUESUTILITY_HPP
#define ZOLTAN2_GRAPHICMETRICVALUESUTILITY_HPP

#include <Zoltan2_Directory_Impl.hpp>
#include <Zoltan2_ImbalanceMetrics.hpp>
#include <Zoltan2_MetricUtility.hpp>
#include <zoltan_dd.h>
#include <Zoltan2_TPLTraits.hpp>
#include <map>
#include <vector>

namespace Zoltan2 {

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
 * globalWeightedByPart() must be called by all processes in \c comm.
 */
template <typename Adapter,
          typename MachineRep =   // Default MachineRep type
                  MachineRepresentation<typename Adapter::scalar_t,
                                        typename Adapter::part_t> >
void globalWeightedByPart(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graph,
    const ArrayView<const typename Adapter::part_t> &parts,
    typename Adapter::part_t &numParts,
    ArrayRCP<RCP<BaseClassMetrics<typename Adapter::scalar_t> > > &metrics,
    ArrayRCP<typename Adapter::scalar_t> &globalSums,
    bool bMessages = true,
    const RCP <const MachineRep> machine = Teuchos::null) { 

  env->timerStart(MACRO_TIMERS, "globalWeightedByPart");

  // Note we used to have with hops as a separate method but decided to combine
  // both into this single method. machine is an optional parameter to choose
  // between the two methods.
  bool bHops = (machine != Teuchos::null);

  env->debug(DETAILED_STATUS, "Entering globalWeightedByPart");

  //////////////////////////////////////////////////////////
  // Initialize return values

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::part_t part_t;

  typedef typename Zoltan2::GraphModel<typename Adapter::base_adapter_t>::input_t 
    t_input_t;

  lno_t localNumVertices = graph->getLocalNumVertices();
  lno_t localNumEdges = graph->getLocalNumEdges();

  ArrayView<const gno_t> Ids;
  ArrayView<t_input_t> v_wghts;
  graph->getVertexList(Ids, v_wghts);

  typedef GraphMetrics<scalar_t> gm_t;

  // get the edge ids, and weights
  ArrayView<const gno_t> edgeIds;
  ArrayView<const offset_t> offsets;
  ArrayView<t_input_t> e_wgts;
  graph->getEdgeList(edgeIds, offsets, e_wgts);


  std::vector <scalar_t> edge_weights;
  int numWeightPerEdge = graph->getNumWeightsPerEdge();

  int numMetrics = bHops ?
    4 :                      // "edge cuts", messages, hops, weighted hops
    2;                       // "edge cuts", messages

  if (numWeightPerEdge) numMetrics += bHops ?
    numWeightPerEdge * 2 :   // "weight n", weighted hops per weight n
    numWeightPerEdge;        // "weight n"

  // add some more metrics to the array
  auto next = metrics.size(); // where we begin filling
  for (auto n = 0; n < numMetrics; ++n)  {
    addNewMetric<gm_t, scalar_t>(env, metrics);
  }

  std::vector <part_t> e_parts (localNumEdges);

  std::vector<part_t> local_parts;

#ifdef HAVE_ZOLTAN2_MPI
  if (comm->getSize() > 1) {
    const bool bUseLocalIDs = false;  // Local IDs not needed
    typedef Zoltan2_Directory_Simple<gno_t,lno_t,part_t> directory_t;
    int debug_level = 0;
    directory_t directory(comm, bUseLocalIDs, debug_level);

    if (localNumVertices)
      directory.update(localNumVertices, &Ids[0], NULL, &parts[0],
        NULL, directory_t::Update_Mode::Replace);
    else
      directory.update(localNumVertices, NULL, NULL, NULL,
        NULL, directory_t::Update_Mode::Replace);

    if (localNumEdges)
      directory.find(localNumEdges, &edgeIds[0], NULL, &e_parts[0],
        NULL, NULL, false);
    else
      directory.find(localNumEdges, NULL, NULL, NULL,
        NULL, NULL, false);
  } else
#endif
  {
    std::map<gno_t,lno_t> global_id_to_local_index;

    // else everything is local.
    // we need a globalid to local index conversion.
    // this does not exists till this point, so we need to create one.

    for (lno_t i = 0; i < localNumVertices; ++i){
      //at the local index i, we have the global index Ids[i].
      //so write i, to Ids[i] index of the vector.
      global_id_to_local_index[Ids[i]] = i;
    }

    for (lno_t i = 0; i < localNumEdges; ++i){
      gno_t ei = edgeIds[i];
      //ei is the global index of the neighbor one.
      part_t p = parts[global_id_to_local_index[ei]];
      e_parts[i] = p;
    }
  }

  RCP<const Teuchos::Comm<int> > tcomm = comm;

  env->timerStart(MACRO_TIMERS, "Communication Graph Create");
  {
    const bool bUseLocalIDs = false;  // Local IDs not needed
    int debug_level = 0;
    // Create a directory indexed by part_t with values t_scalar_t for weight sums

    // this struct is the user data type for a part
    // each part will have a std::vector<part_info> for its user data
    // representing the list of all neighbors and a weight for each.
    struct part_info {
      part_info() : sum_weights(0) {
      }

      // operator +=
      // this allows the directory to know how to assemble two structs
      // which return true for ==.
      // TODO: Decide if we want directory to work like this for AggregateAdd
      const part_info & operator+=(const part_info & src) {
        sum_weights += src.sum_weights;
        return *this;   // return old value
      }

      // operator>
      // TODO: Decide if we want directory to work like this for AggregateAdd
      bool operator>(const part_info & src) {
        // Note: Currently this doesn't actually do anything except allow this
        // struct to compile. Aggregate mode used operator> to preserve ordering
        // and therefore a custom struct must currently define it. However in
        // this test we will be using AggregateAdd mode which doesn't actually
        // use this operator> . However if we change the test so we require the
        // input data to already be ordered by target_part, then we could change
        // the directory implementation so AggregateAdd and Aggregate are almost
        // identical. The only difference would be that Aggregate mode would
        // check operator== and if true, throw away the duplicate, while
        // AggregateAdd mode would check operator== and if true, call operator+=
        // to combine the values, in this case summing sum_weights.
        return (target_part > src.target_part);
      }

      // operator==
      // this allows the directory to know that two structs should be
      // aggregated into one using the operator +=.
      // TODO: Decide if we want directory to work like this for AggregateAdd
      // This works but seems fussy/complicated - to discuss. I'm not yet sure
      // how to best integrate this so we can aggregate both simple types where
      // we just keep unique elements and more complex structs with a 'rule'
      // for combining them.
      bool operator==(const part_info & src) {
        // if target_part is the same then the values for sum_weights will
        // be summed.
        return (target_part == src.target_part);
      }

      part_t target_part;     // the part this part_info refers to
      scalar_t sum_weights; // the sum of weights
    };

    // get the vertices in each part in my part.
    std::vector <lno_t> part_begins(numParts, -1);
    std::vector <lno_t> part_nexts(localNumVertices, -1);
 
    // cluster vertices according to their parts.
    // create local part graph.
    for (lno_t i = 0; i < localNumVertices; ++i){
      part_t ap = parts[i];
      part_nexts[i] = part_begins[ap];
      part_begins[ap] = i;
    }

    for (int weight_index = -1; 
         weight_index < numWeightPerEdge; ++weight_index) {

      std::vector<part_t> part_data(numParts); // will resize to lower as needed
      std::vector<std::vector<part_info>> user_data(numParts); // also to resize
      int count_total_entries = 0;

      std::vector <part_t> part_neighbors(numParts);
      std::vector <scalar_t> part_neighbor_weights(numParts, 0);
      std::vector <scalar_t> part_neighbor_weights_ordered(numParts);

      // coarsen for all vertices in my part in order with parts.
      for (part_t i = 0; i < numParts; ++i) {
        part_t num_neighbor_parts = 0;
        lno_t v = part_begins[i];
        // get part i, and first vertex in this part v.
        while (v != -1){
          // now get the neightbors of v.
          for (offset_t j = offsets[v]; j < offsets[v+1]; ++j){

            // get the part of the second vertex.
            part_t ep = e_parts[j];

            // TODO: Can we skip condition (i==ep)
            // The self reference set is going to be excluded later anyways
            // so we could make this more efficient.
            scalar_t ew = 1;
            if (weight_index > -1){
              ew = e_wgts[weight_index][j];
            }
            // add it to my local part neighbors for part i.
            if (part_neighbor_weights[ep] < 0.00001){
              part_neighbors[num_neighbor_parts++] = ep;
            }
            part_neighbor_weights[ep] += ew;
          }
          v = part_nexts[v];
        }

        // now get the part list.
        for (lno_t j = 0; j < num_neighbor_parts; ++j) {
          part_t neighbor_part = part_neighbors[j];
          part_neighbor_weights_ordered[j] = 
            part_neighbor_weights[neighbor_part];
          part_neighbor_weights[neighbor_part] = 0;
        }

        if (num_neighbor_parts > 0) {
          // for the new directory note a difference in the logic flow
          // originally we have CrsMatrix which could collect these values
          // as we built each row. For the directory it's probably better to
          // have update called just once so we collect the values and then
          // do all of the update at the end.
          part_data[count_total_entries] = i; // set up for directory
          std::vector<part_info> & add_user_data = 
            user_data[count_total_entries];
          ++count_total_entries;

          add_user_data.resize(num_neighbor_parts);

          for(int n = 0; n < num_neighbor_parts; ++n) {
            part_info & add_data = add_user_data[n];
            add_data.target_part = part_neighbors[n];
            add_data.sum_weights = part_neighbor_weights_ordered[n];
          }
        }
      }
 
      scalar_t max_edge_cut = 0;
      scalar_t total_edge_cut = 0;
      part_t max_message = 0;
      part_t total_message = 0;

      part_t total_hop_count = 0;
      scalar_t total_weighted_hop_count = 0;
      part_t max_hop_count = 0;
      scalar_t max_weighted_hop_count = 0;

      // for serial or comm size 1 we need to fill this from local data
      // TODO: Maybe remove all special casing for serial and make this pipeline
      // uniform always
      if(local_parts.size() == 0) {
        local_parts.resize(numParts);
        for(size_t n = 0; n < local_parts.size(); ++n) {
          local_parts[n] = n;
        }
      }

      std::vector<std::vector<part_info>> find_user_data;

      // directory does not yet support SerialComm because it still has older
      // MPI calls which need to be refactored to Teuchos::comm format. To
      // work around this issue skip the directory calls and just set the
      // find data equal to the input update data. This works because above
      // logic has already summed the weights per process so in the SerialComm
      // case, there won't be duplicates.
      bool bSerialComm =
        (dynamic_cast<const Teuchos::SerialComm<int>*>(&(*comm)) != NULL);

      if(!bSerialComm) {
        typedef Zoltan2_Directory_Vector<part_t,int,std::vector<part_info>>
          directory_t;
        directory_t directory(comm, bUseLocalIDs, debug_level);

        if(count_total_entries) {
          // update
          directory.update(count_total_entries, &part_data[0], 
                           NULL, &user_data[0],
            NULL, directory_t::Update_Mode::AggregateAdd);
        }
        else {
          directory.update(count_total_entries, NULL, NULL, NULL,
            NULL, directory_t::Update_Mode::AggregateAdd);
        }

        // get my local_parts (parts managed on this directory)
        directory.get_locally_managed_gids(local_parts);

        // set up find_user_data to have the right size
        find_user_data.resize(local_parts.size());

        // find
        directory.find(local_parts.size(), &local_parts[0], NULL,
          &find_user_data[0], NULL, NULL, false);
      }
      else {
        find_user_data = user_data;
      }

      for(size_t n = 0; n < local_parts.size(); ++n) {
          scalar_t part_edge_cut = 0;
          part_t part_messages = 0;
          const std::vector<part_info> & data = find_user_data[n];
          for(size_t q = 0; q < data.size(); ++q) {
            const part_t local_part = local_parts[n];
            const part_info & info = data[q];
            if (info.target_part != local_part) {
              part_edge_cut += info.sum_weights;
              part_messages += 1;

              if(bHops) {
                typename MachineRep::machine_pcoord_t hop_count = 0;
                machine->getHopCount(local_part, info.target_part, hop_count);

                part_t hop_counts = hop_count;
                scalar_t weighted_hop_counts = hop_count * info.sum_weights;

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
          }

          if(part_edge_cut > max_edge_cut) {
            max_edge_cut = part_edge_cut;
          }
          total_edge_cut += part_edge_cut;

          if (part_messages > max_message){
            max_message = part_messages;
          }
          total_message += part_messages;
      }


      scalar_t g_max_edge_cut = 0;
      scalar_t g_total_edge_cut = 0;
      part_t g_max_message = 0;
      part_t g_total_message = 0;

      part_t g_total_hop_count = 0;
      scalar_t g_total_weighted_hop_count = 0;
      part_t g_max_hop_count = 0;
      scalar_t g_max_weighted_hop_count = 0;

      try{
        Teuchos::reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_MAX, 1, 
                                          &max_edge_cut, &g_max_edge_cut);
        Teuchos::reduceAll<int, part_t>(*comm, Teuchos::REDUCE_MAX, 1, 
                                        &max_message, &g_max_message);

        Teuchos::reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_SUM, 1, 
                                          &total_edge_cut, &g_total_edge_cut);
        Teuchos::reduceAll<int, part_t>(*comm, Teuchos::REDUCE_SUM, 1, 
                                        &total_message, &g_total_message);

        if(bHops) {
          Teuchos::reduceAll<int, part_t>(*comm, Teuchos::REDUCE_MAX, 1, 
                                          &max_hop_count, &g_max_hop_count);
          Teuchos::reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_MAX, 1, 
                                            &max_weighted_hop_count,
                                            &g_max_weighted_hop_count);

          Teuchos::reduceAll<int, part_t>(*comm, Teuchos::REDUCE_SUM, 1, 
                                          &total_hop_count, &g_total_hop_count);
          Teuchos::reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_SUM, 1,
                                            &total_weighted_hop_count,
                                            &g_total_weighted_hop_count);
        }
      }
      Z2_THROW_OUTSIDE_ERROR(*env);
 
      if (weight_index == -1){
        metrics[next]->setName("edge cuts");
      }
      else {
        std::ostringstream oss;
        oss << "weight " << weight_index;
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

      if(bHops) {
        if (weight_index == -1){
          metrics[next]->setName("hops (No Weight)");
          metrics[next]->setMetricValue("global maximum", g_max_hop_count);
          metrics[next]->setMetricValue("global sum", g_total_hop_count);
          next++;
        }

        std::ostringstream oss;
        oss << "weighted hops" << weight_index;
        metrics[next]->setName( oss.str());
        metrics[next]->
          setMetricValue("global maximum", g_max_weighted_hop_count);
        metrics[next]->
          setMetricValue("global sum", g_total_weighted_hop_count);
        next++;
      }
    }
  }

  env->timerStop(MACRO_TIMERS, "globalWeightedByPart");

  env->debug(DETAILED_STATUS, "Exiting globalWeightedByPart");
}

/*! \brief Print out header info for graph metrics.
 */
template <typename scalar_t, typename part_t>
void printGraphMetricsHeader(std::ostream &os, 
                             part_t targetNumParts, 
                             part_t numParts) {

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
void printGraphMetrics(std::ostream &os, 
                       part_t targetNumParts, 
                       part_t numParts, 
                       const ArrayView<RCP<BaseClassMetrics<scalar_t> > > 
                         &infoList) {

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
void printGraphMetrics(std::ostream &os, 
                       part_t targetNumParts, 
                       part_t numParts, 
                       RCP<BaseClassMetrics<scalar_t>> metricValue) {

  printGraphMetricsHeader<scalar_t, part_t>(os, targetNumParts, numParts);
  metricValue->printLine(os);
}

} //namespace Zoltan2

#endif
