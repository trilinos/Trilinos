// @HEADER
// //
// // ***********************************************************************
// //
// //        MueLu: A package for multigrid based preconditioning
// //                  Copyright 2012 Sandia Corporation
// //
// // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// // the U.S. Government retains certain rights in this software.
// //
// // Redistribution and use in source and binary forms, with or without
// // modification, are permitted provided that the following conditions are
// // met:
// //
// // 1. Redistributions of source code must retain the above copyright
// // notice, this list of conditions and the following disclaimer.
// //
// // 2. Redistributions in binary form must reproduce the above copyright
// // notice, this list of conditions and the following disclaimer in the
// // documentation and/or other materials provided with the distribution.
// //
// // 3. Neither the name of the Corporation nor the names of the
// // contributors may be used to endorse or promote products derived from
// // this software without specific prior written permission.
// //
// // THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// // EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// // PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// // CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// // EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// // PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// // PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// // LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// // NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// // SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// //
// // Questions? Contact
// //                    Jonathan Hu       (jhu@sandia.gov)
// //                    Andrey Prokopenko (aprokop@sandia.gov)
// //                    Ray Tuminaro      (rstumin@sandia.gov)
// //
// // ***********************************************************************
// //
// // @HEADER
/*
 * Xpetra_RegionHandler_def.hpp
 *
 * Created on: August 17, 2017
 * 	Author: Massimiliano Lupo Pasini (massimiliano.lupo.pasini@gmail.com)
 */
#ifndef XPETRA_REGIONHANDLER_DEF_HPP
#define XPETRA_REGIONHANDLER_DEF_HPP

#include "Xpetra_RegionHandler_decl.hpp"
#include <algorithm>

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>::RegionHandler(const std::string &file_name, RCP<const Teuchos::Comm<int> > comm)
  : comm_(comm) {
  ReadFileInfo(file_name);

  // Nodes are shuffled so that regions are sorted in ascending labeling order
  std::sort(nodes_.begin(), nodes_.end(), compareRegions<GlobalOrdinal>);
  nodes_sorted_by_regions_ = true;

  if (comm_->getRank() == 0)
    std::cout << "Started NodesToRegion" << std::endl;
  NodesToRegion();
  if (comm_->getRank() == 0)
    std::cout << "Finished NodesToRegion" << std::endl;
  ComputeProcRegions();
  if (comm_->getRank() == 0)
    std::cout << "Started RowMaps" << std::endl;
  CreateRowMaps();
  if (comm_->getRank() == 0)
    std::cout << "Finished RowMaps" << std::endl;

  num_region_nodes_.clear();
  num_region_nodes_.resize(num_total_regions_);

  // For each region, the following loop counts the number of region nodes and stores them
  for (GlobalOrdinal region_idx = 1; region_idx <= num_total_regions_; ++region_idx) {
    checkerNode<GlobalOrdinal> unaryPredicate(region_idx);
    typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator nodes_iterator1;
    typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator nodes_iterator2;
    nodes_iterator1                   = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerNode<GlobalOrdinal> >(nodes_.begin(), nodes_.end(), unaryPredicate);
    nodes_iterator2                   = std::find_if_not<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerNode<GlobalOrdinal> >(nodes_iterator1, nodes_.end(), unaryPredicate);
    num_region_nodes_[region_idx - 1] = nodes_iterator2 - nodes_iterator1;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadFileInfo(const std::string &file_name) {
  std::ifstream input_file_(file_name, std::ifstream::in);
  std::string line;
  TEUCHOS_TEST_FOR_EXCEPTION(!input_file_.good(), Exceptions::RuntimeError, "Cannot read \"" << file_name << "\"");

  GlobalOrdinal line_index = 0;

  // The information contained in the file is imported and stored in a Teuchos::Array of tuples.
  // The first field of the tuple is the composite node index, the second field of the tuple is the region index
  while (std::getline(input_file_, line)) {
    std::istringstream is(line);
    GlobalOrdinal number;
    Array<GlobalOrdinal> node;
    std::tuple<GlobalOrdinal, GlobalOrdinal> node_region;
    Array<GlobalOrdinal> composite_info;

    node.clear();
    composite_info.clear();

    if (1 == line_index) {
      while (is >> number)
        composite_info.push_back(number);

      TEUCHOS_TEST_FOR_EXCEPTION(composite_info.size() != 2, Exceptions::RuntimeError, "The composite information must be a couple of integers: nTotal of nodes + nTotal of regions \n");
      num_total_nodes_   = composite_info[0];
      num_total_regions_ = composite_info[1];
    } else if (line_index > 2) {
      while (is >> number) {
        node.push_back(number);
      }
      TEUCHOS_TEST_FOR_EXCEPTION(node.size() != 2, Exceptions::RuntimeError, "The node information must be a couple of integers: Node index + Region idnex \n");
      node_region = std::make_tuple(node[0], node[1]);
      nodes_.push_back(node_region);
      node.clear();
    }
    line_index++;
  }
  input_file_.close();
}

// This routines computes the way regions are partitioned across processes
// The partitioning policies of course depends on whether the number of processes
// exceeds the number of regions or not.
// ASSUMPTION: A PROCESS CANNOT OWN CHUNKS OF MULTIPLE REGIONS. EITHER A PROCESS IS CONFINED INSIDE A SINGLE REGION
// OR IT MUST POSSESS ENTIRE REGIONS.
// The distribution of regions (or portions of them) across processes is conductes so to guarantee load balancing
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ComputeProcRegions() {
  int tot_num_proc = comm_->getSize();
  int myPID        = comm_->getRank();

  regions_per_proc_.clear();
  procs_per_region_.clear();

  // If the number of processes instantiate is smaller than the total number of regions,
  // then each process owns entire regions. The number of regions per process is calculates so to guarantee
  // load balancing. After an initial distribution of regions, leftover regions that have not been assigned to any process yet are
  // distributed in a round-robin fashion
  if (tot_num_proc < num_total_regions_) {
    int min_nregions_proc    = std::floor(static_cast<double>(num_total_regions_) / static_cast<double>(tot_num_proc));
    int num_leftover_regions = num_total_regions_ % tot_num_proc;

    for (int i = 1; i <= min_nregions_proc; ++i)
      regions_per_proc_.push_back(myPID * min_nregions_proc + i);

    if (num_leftover_regions >= myPID + 1 && num_leftover_regions != 0)
      regions_per_proc_.push_back(min_nregions_proc * tot_num_proc + (myPID + 1));

    for (int procID = 0; procID < comm_->getSize(); ++procID) {
      Array<GlobalOrdinal> proc;
      proc.clear();
      proc.push_back(procID);
      for (int i = 1; i <= min_nregions_proc; ++i) {
        GlobalOrdinal region_index                                 = (procID)*min_nregions_proc + i;
        std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > tuple_aux = std::make_tuple(region_index, proc);
        procs_per_region_.push_back(tuple_aux);
      }

      if (num_leftover_regions >= procID + 1 && num_leftover_regions != 0) {
        GlobalOrdinal region_index                                 = min_nregions_proc * tot_num_proc + (procID + 1);
        std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > tuple_aux = std::make_tuple(region_index, proc);
        procs_per_region_.push_back(tuple_aux);
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!(procs_per_region_.size() == num_total_regions_), Exceptions::RuntimeError, "PID: " << comm_->getRank() << " - Number of regions detected does not match with the initially declared one \n procs_per_region_ tracks " << procs_per_region_.size() << " regions whereas num_total_regions_ = " << num_total_regions_ << "\n");
  }
  // This is easy: if the number of regions coincides with the total number of processes instantiated, then
  // a one-to-one relation between processes and regions is created
  else if (tot_num_proc == num_total_regions_) {
    regions_per_proc_.push_back(myPID + 1);

    for (int i = 0; i < num_total_regions_; ++i) {
      GlobalOrdinal region_index = i + 1;
      Array<GlobalOrdinal> proc;
      proc.clear();
      proc.push_back(i);
      std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > tuple_aux = std::make_tuple(region_index, proc);
      procs_per_region_.push_back(tuple_aux);
    }
  }
  // If the number of processes exceeds the number of regions in the domain,
  // then each process is given a subset of a region.
  // N.B.: A SINGLE PROCESS IS NOT ALLOWED TO OWN CHUNCKS OF MULTIPLE REGIONS.
  // IN THIS CONFIGURATION EACH PROCESS IS CONFINED TO A SINGLE REGION
  else if (tot_num_proc > num_total_regions_) {
    int num_procs_region       = std::ceil(static_cast<double>(tot_num_proc) / static_cast<double>(num_total_regions_));
    int num_regions_extra_proc = tot_num_proc % num_total_regions_;
    int proc_count             = 0;
    std::tuple<int, Array<GlobalOrdinal> > region_tuple;

    for (int i = 1; i <= num_total_regions_; ++i) {
      Array<GlobalOrdinal> procs;
      procs.clear();
      if (i <= num_regions_extra_proc || num_regions_extra_proc == 0)
        for (int j = 1; j <= num_procs_region; ++j) {
          procs.push_back(proc_count);
          proc_count++;
        }
      else
        for (int j = 1; j <= num_procs_region - 1; ++j) {
          procs.push_back(proc_count);
          proc_count++;
        }
      std::sort(procs.begin(), procs.end());
      region_tuple = std::make_tuple(i, procs);
      procs_per_region_.push_back(region_tuple);
    }
    regions_per_proc_.clear();
  }
}

// This routine associates a globally indexed node with the list of regions it belongs to
// This is helpful to spot which nodes lie on a interregion interface. In fact, these nodes must have
// the list of regions with more than one element
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NodesToRegion() {
  nodesToRegion_.clear();
  interfaceNodes_.clear();
  Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > nodes_reordered;
  nodes_reordered = nodes_;
  std::sort(nodes_reordered.begin(), nodes_reordered.end(), compareNodes<GlobalOrdinal>);

  typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator node_iterator;
  node_iterator = nodes_reordered.begin();

  while (node_iterator != nodes_reordered.end()) {
    GlobalOrdinal current_node = std::get<0>(*(node_iterator));
    Array<GlobalOrdinal> regions;
    regions.clear();
    regions.push_back(std::get<1>(*node_iterator));

    typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator next_node_iterator = node_iterator + 1;

    while (next_node_iterator != nodes_reordered.end()) {
      GlobalOrdinal next_node = std::get<0>(*(next_node_iterator));
      if (current_node == next_node) {
        // As long as the information spanned regards the same node,
        // the algorithm keeps on increasing the list of regions the given mesh node belong to
        regions.push_back(std::get<1>(*(next_node_iterator)));
        next_node_iterator++;
      } else {
        // When the mesh node label changes, then the algorithm
        // stops recording information about the previous node and it starts recording information for the new one
        node_iterator = next_node_iterator;
        break;
      }
    }
    std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > new_tuple;
    std::sort(regions.begin(), regions.end());
    new_tuple = std::make_tuple(current_node, regions);
    nodesToRegion_.push_back(new_tuple);

    if (regions.size() > 1)
      interfaceNodes_.push_back(new_tuple);

    if (next_node_iterator == nodes_reordered.end())
      break;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(!(nodesToRegion_.size() == num_total_nodes_), Exceptions::RuntimeError, "Number of nodes detected does not match with the initially declared one \n"
                                                                                                         << "nodesToRegion tracks " << nodesToRegion_.size() << " nodes whereas num_total_nodes_ =" << num_total_nodes_ << "\n");
}

// This routine creates row maps for composite and region matrices
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>::CreateRowMaps() {
  TEUCHOS_TEST_FOR_EXCEPTION((procs_per_region_.empty() && regions_per_proc_.empty()), Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Information about region partitioning across processors is not consistent: incorrect values for number of processors or number of regions \n");
  Array<GlobalOrdinal> elements;
  Array<GlobalOrdinal> region_elements;
  Array<Array<GlobalOrdinal> > elements_per_region;
  Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > regionToAll;
  int myPID = comm_->getRank();

  elements.clear();
  region_elements.clear();
  elements_per_region.clear();
  elements_per_region.resize(num_total_regions_);
  regionToAll.clear();
  regionToAll.resize(num_total_regions_);

  TEUCHOS_TEST_FOR_EXCEPTION(!nodes_sorted_by_regions_, Exceptions::RuntimeError, "Nodes are not sorted by regions in ascending order \n");
  TEUCHOS_TEST_FOR_EXCEPTION(num_total_nodes_ > nodes_.size(), Exceptions::RuntimeError, "Number of nodes declared in input file does not match with the effective number of nodes provided\n"
                                                                                             << "num_total_nodes_ =" << num_total_nodes_ << " whereas nodes_ tracks " << nodes_.size() << " nodes \n");
  TEUCHOS_TEST_FOR_EXCEPTION(num_total_regions_ != std::get<1>(*(nodes_.end() - 1)), Exceptions::RuntimeError, "Number of regions declared in input file does not match with the effective number of regions provided\n");

  if (!(regions_per_proc_.empty())) {
    typename Array<GlobalOrdinal>::iterator iter_array;
    for (iter_array = regions_per_proc_.begin(); iter_array != regions_per_proc_.end(); ++iter_array) {
      region_elements.clear();
      checkerNode<GlobalOrdinal> unaryPredicate(*iter_array);
      typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator nodes_iterator1;
      typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator nodes_iterator2;

      // We position an iterator at the beginning of the information associated with the a region owned by the calling process
      // and another iterator right at the end of the information
      nodes_iterator1 = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerNode<GlobalOrdinal> >(nodes_.begin(), nodes_.end(), unaryPredicate);
      nodes_iterator2 = std::find_if_not<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerNode<GlobalOrdinal> >(nodes_iterator1, nodes_.end(), unaryPredicate);

      // Coun the number of mesh nodes inside a region
      int num_region_nodes = nodes_iterator2 - nodes_iterator1;

      typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator nodes_iterator_aux;

      // The policy assumes that in the input file the indexBase for the node label is 1
      GlobalOrdinal region_node_label = 1;
      for (nodes_iterator_aux = nodes_iterator1; nodes_iterator_aux != nodes_iterator2; ++nodes_iterator_aux) {
        GlobalOrdinal node = std::get<0>(*nodes_iterator_aux);
        checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
        typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
        nodes_to_region_iterator           = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);
        Array<GlobalOrdinal> nodal_regions = std::get<1>(*nodes_to_region_iterator);

        // By default, I choose that a node is owned by the process associated with the region that shows up first in its list of beloning
        // This guarantees that each row of the composite stiffness matrix is owned only by a single process, as Trilinos requires
        if (*iter_array == nodal_regions[0])
          elements.push_back(node);

        // Nodes on the interface still belong to multiple regions, so
        // it is important to keep track of this for the row maps of region matrices
        region_elements.push_back(region_node_label);

        // If a process owns a region (or even a portion of it), we provide to it a map
        // from region indices to composite indices for all the nodes inside that region,
        // even if a specific node is not owned by the calling process
        regionToAll[*iter_array - 1].push_back(std::make_tuple(region_node_label, std::get<0>(*nodes_iterator_aux)));
        region_node_label++;
      }

      // C++ indexing starts from 0, so everything is shofted backward by one to make it consistent with programming language's policies
      for (typename Array<GlobalOrdinal>::iterator iter = region_elements.begin(); iter != region_elements.end(); ++iter)
        *iter = *iter - 1;

      elements_per_region[*iter_array - 1] = region_elements;
      TEUCHOS_TEST_FOR_EXCEPTION((num_region_nodes != regionToAll[*iter_array - 1].size()), Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Number of region nodes does not match with number of nodes stored in regionToAll \n"
                                                                                                                                     << "num_region_nodes= " << num_region_nodes << " whereas regionToAll[" << *iter_array - 1 << "].size()= " << regionToAll[*iter_array - 1].size() << "\n");
    }
    TEUCHOS_TEST_FOR_EXCEPTION((num_total_regions_ != regionToAll.size()), Exceptions::RuntimeError, "regionToAll size has been corrupted\n"
                                                                                                         << "num_total_regions_ = " << num_total_regions_ << " whereas regionToAll.size() = " << regionToAll.size() << "\n");
  } else {
    // The code enters the scope of these curly brackets if the number of processes exceeds the number of regions in the domain.
    // Therefore, each process owns only a subregion (or at most a single entire region).
    // In this situation, the calling process must identify the region it is associated with
    bool region_found      = false;
    GlobalOrdinal myRegion = -1;
    TEUCHOS_TEST_FOR_EXCEPTION(!(procs_per_region_.size() == num_total_regions_), Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - Number of total regions does not match with regionHandler structures \n");

    Array<GlobalOrdinal> region_procs;
    while (!region_found) {
      typename Array<GlobalOrdinal>::iterator iter_proc;
      for (GlobalOrdinal region_index = 1; region_index <= procs_per_region_.size(); ++region_index) {
        region_procs = std::get<1>(procs_per_region_[region_index - 1]);
        iter_proc    = std::find(region_procs.begin(), region_procs.end(), myPID);
        if (iter_proc != region_procs.end()) {
          myRegion     = region_index;
          region_found = true;
        }
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION((myRegion == -1 || !region_found), Exceptions::RuntimeError, ("Region containing PROC ID: " + std::to_string(myPID) + " NOT FOUND \n"));
    region_procs = std::get<1>(procs_per_region_[myRegion - 1]);

    checkerNode<GlobalOrdinal> unaryPredicate(myRegion);
    typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator nodes_iterator1;
    typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator nodes_iterator2;
    nodes_iterator1 = std::find_if<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerNode<GlobalOrdinal> >(nodes_.begin(), nodes_.end(), unaryPredicate);
    nodes_iterator2 = std::find_if_not<typename Array<std::tuple<GlobalOrdinal, GlobalOrdinal> >::iterator, checkerNode<GlobalOrdinal> >(nodes_iterator1, nodes_.end(), unaryPredicate);

    int num_region_nodes = nodes_iterator2 - nodes_iterator1;
    int num_region_procs = region_procs.size();

    if (num_region_nodes < num_region_procs) {
      Array<GlobalOrdinal> region_procs_reduced;
      region_procs_reduced.clear();
      for (int i = 0; i < num_region_nodes; ++i)
        region_procs_reduced.push_back(region_procs[i]);

      typename Array<GlobalOrdinal>::iterator proc_iterator;
      proc_iterator = std::find<typename Array<GlobalOrdinal>::iterator, GlobalOrdinal>(region_procs_reduced.begin(), region_procs_reduced.end(), myPID);

      if (proc_iterator != region_procs_reduced.end())  // This reasoning works because the PROC ID for each region has been previously sorted in ascending order
      {
        GlobalOrdinal node              = std::get<0>(*(nodes_iterator1 + (proc_iterator - region_procs_reduced.begin() + 1)));
        GlobalOrdinal region_node_label = proc_iterator - region_procs_reduced.begin() + 1;
        checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
        typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
        nodes_to_region_iterator           = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);
        Array<GlobalOrdinal> nodal_regions = std::get<1>(*nodes_to_region_iterator);

        // The follolwing if statement is necessary to guarantee uniqueness of the RowMap used for the composite matrix
        if (myRegion == nodal_regions[0])
          elements.push_back(node);

        // Although a process does not own a row in the composite matrix, it may still happen that it owns the row
        // from a region matrix perspective
        region_elements.push_back(region_node_label);
      }

      // If a process owns a region (or even a portion of it), we provide to it a map
      // from region indices to composite indices for all the nodes inside that region,
      // even if a specific node is not owned by the calling process
      // If a process owns something of a region, then the process has a global view of who owns what for that region
      // Although this may seem more information than what actually needed, it is important for the computation of the collapsing.
      // If the collapsing is not calculated, then this structure actually overestimates what a process needs to know.
      for (proc_iterator = region_procs_reduced.begin(); proc_iterator != region_procs_reduced.end(); ++proc_iterator) {
        GlobalOrdinal node              = std::get<0>(*(nodes_iterator1 + (proc_iterator - region_procs_reduced.begin() + 1)));
        GlobalOrdinal region_node_label = proc_iterator - region_procs_reduced.begin() + 1;
        checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
        typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
        nodes_to_region_iterator = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);
        regionToAll[myRegion - 1].push_back(std::make_tuple(region_node_label, node));
      }

    } else if (num_region_nodes == num_region_procs) {
      typename Array<GlobalOrdinal>::iterator proc_iterator;
      proc_iterator = std::find<typename Array<GlobalOrdinal>::iterator, GlobalOrdinal>(region_procs.begin(), region_procs.end(), myPID);

      if (proc_iterator != region_procs.end())  // This reasoning works because the PROC ID for each region has been previously sorted in ascending order
      {
        GlobalOrdinal node              = std::get<0>(*(nodes_iterator1 + (proc_iterator - region_procs.begin() + 1)));
        GlobalOrdinal region_node_label = proc_iterator - region_procs.begin() + 1;
        checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
        typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
        nodes_to_region_iterator           = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);
        Array<GlobalOrdinal> nodal_regions = std::get<1>(*nodes_to_region_iterator);
        if (myRegion == nodal_regions[0])
          elements.push_back(node);

        region_elements.push_back(region_node_label);
      }

      // If a process owns a region (or even a portion of it), we provide to it a map
      // from region indices to composite indices for all the nodes inside that region,
      // even if a specific node is not owned by the calling process
      for (proc_iterator = region_procs.begin(); proc_iterator != region_procs.end(); ++proc_iterator) {
        GlobalOrdinal node              = std::get<0>(*(nodes_iterator1 + (proc_iterator - region_procs.begin() + 1)));
        GlobalOrdinal region_node_label = proc_iterator - region_procs.begin() + 1;
        checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
        typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
        nodes_to_region_iterator = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);
        regionToAll[myRegion - 1].push_back(std::make_tuple(region_node_label, node));
      }
    } else {
      typename Array<GlobalOrdinal>::iterator proc_iterator;
      proc_iterator = std::find<typename Array<GlobalOrdinal>::iterator, GlobalOrdinal>(region_procs.begin(), region_procs.end(), myPID);

      int num_nodes_proc       = std::ceil(static_cast<double>(num_region_nodes) / static_cast<double>(num_region_procs));
      int num_procs_extra_node = num_region_nodes % num_region_procs;

      if (proc_iterator - region_procs.begin() + 1 <= num_procs_extra_node || num_procs_extra_node == 0) {
        int init_node = num_nodes_proc * (proc_iterator - region_procs.begin());
        for (int i = 0; i < num_nodes_proc; ++i) {
          GlobalOrdinal node              = std::get<0>(*(nodes_iterator1 + init_node + i));
          GlobalOrdinal region_node_label = init_node + i + 1;
          checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
          typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
          nodes_to_region_iterator           = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);
          Array<GlobalOrdinal> nodal_regions = std::get<1>(*nodes_to_region_iterator);
          if (myRegion == nodal_regions[0])
            elements.push_back(node);

          region_elements.push_back(region_node_label);
        }
      } else {
        int init_node = num_nodes_proc * num_procs_extra_node + (proc_iterator - region_procs.begin() - num_procs_extra_node) * (num_nodes_proc - 1);
        for (int i = 0; i < num_nodes_proc - 1; ++i) {
          GlobalOrdinal node              = std::get<0>(*(nodes_iterator1 + init_node + i));
          GlobalOrdinal region_node_label = init_node + i + 1;
          checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
          typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
          nodes_to_region_iterator           = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);
          Array<GlobalOrdinal> nodal_regions = std::get<1>(*nodes_to_region_iterator);
          if (myRegion == nodal_regions[0])
            elements.push_back(node);

          region_elements.push_back(region_node_label);
        }
      }

      // If a process owns a region (or even a portion of it), we provide to it a map
      // from region indices to composite indices for all the nodes inside that region,
      // even if a specific node is not owned by the calling process
      for (proc_iterator = region_procs.begin(); proc_iterator != region_procs.end(); ++proc_iterator) {
        if (proc_iterator - region_procs.begin() + 1 <= num_procs_extra_node || num_procs_extra_node == 0) {
          int init_node = num_nodes_proc * (proc_iterator - region_procs.begin());
          for (int i = 0; i < num_nodes_proc; ++i) {
            GlobalOrdinal node              = std::get<0>(*(nodes_iterator1 + init_node + i));
            GlobalOrdinal region_node_label = init_node + i + 1;
            checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
            typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
            nodes_to_region_iterator = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);
            regionToAll[myRegion - 1].push_back(std::make_tuple(region_node_label, node));
          }
        } else {
          int init_node = num_nodes_proc * num_procs_extra_node + (proc_iterator - region_procs.begin() - num_procs_extra_node) * (num_nodes_proc - 1);
          for (int i = 0; i < num_nodes_proc - 1; ++i) {
            GlobalOrdinal node              = std::get<0>(*(nodes_iterator1 + init_node + i));
            GlobalOrdinal region_node_label = init_node + i + 1;
            checkerNodesToRegion<GlobalOrdinal> unaryPredicateNode(node);
            typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator nodes_to_region_iterator;
            nodes_to_region_iterator = std::find_if<typename Array<std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > >::iterator, checkerNodesToRegion<GlobalOrdinal> >(nodesToRegion_.begin(), nodesToRegion_.end(), unaryPredicateNode);
            regionToAll[myRegion - 1].push_back(std::make_tuple(region_node_label, node));
          }
        }
      }
    }
    TEUCHOS_TEST_FOR_EXCEPTION((num_total_regions_ != regionToAll.size()), Exceptions::RuntimeError, "Process ID: " << comm_->getRank() << " - regionToAll size has been corrupted\n"
                                                                                                                    << "num_total_regions_ = " << num_total_regions_ << " whereas regionToAll.size()= " << regionToAll.size() << "\n");

    // C++ indexing starts from 0, so everything is shofted backward by one to make it consistent with programming language's policies
    for (typename Array<GlobalOrdinal>::iterator iter = region_elements.begin(); iter != region_elements.end(); ++iter)
      *iter = *iter - 1;

    elements_per_region[myRegion - 1] = region_elements;
  }

  // C++ indexing starts from 0, so everything is shofted backward by one to make it consistent with programming language's policies
  for (typename Array<GlobalOrdinal>::iterator iter = elements.begin(); iter != elements.end(); ++iter)
    *iter = *iter - 1;

  maps_.composite_map_ = elements;
  maps_.region_maps_   = elements_per_region;
  maps_.regionToAll_   = regionToAll;

  for (int i = 0; i < comm_->getSize(); ++i) {
    comm_->barrier();

    if (comm_->getRank() == i) {
      std::cout << "Proc " << i << std::endl;
      for (int i = 0; i < regionToAll.size(); ++i) {
        Teuchos::Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > currArray = regionToAll[i];
        for (int j = 0; j < currArray.size(); ++j) {
          std::cout << std::get<0>(currArray[j]) << "\t" << std::get<1>(currArray[j]) << std::endl;
        }
      }
    }

    comm_->barrier();
  }
}

// Get methods to allow a user to interface with private members of the RegionHandler class
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Array<GlobalOrdinal> RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetRegionRowMap(GlobalOrdinal region_index) const {
  TEUCHOS_TEST_FOR_EXCEPTION(region_index >= num_total_regions_, Exceptions::RuntimeError, "Value of region index exceeds total number of regions stored \n"
                                                                                               << "Trying to access informaiton about region " << region_index << " when the total number of regions stored is " << num_total_regions_ << "\n");
  return maps_.region_maps_[region_index];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetRegionToAll() const {
  return maps_.regionToAll_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetRegionToAll(GlobalOrdinal region_index) const {
  TEUCHOS_TEST_FOR_EXCEPTION(region_index >= num_total_regions_, Exceptions::RuntimeError, "Value of region index exceeds total number of regions stored \n"
                                                                                               << "Trying to access informaiton about region " << region_index << " when the total number of regions stored is " << num_total_regions_ << "\n");
  return maps_.regionToAll_[region_index];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>::printView() const {
  if (0 == comm_->getRank()) {
    std::cout << "Total number of mesh nodes: " << num_total_nodes_ << std::endl;
    std::cout << "Total number of mesh regions: " << num_total_regions_ << std::endl;
    std::cout << "Number of rows in nodes_ structure: " << nodes_.size() << std::endl;
    for (int i = 0; i < nodes_.size(); ++i) {
      std::cout << std::get<0>(nodes_[i]) << "\t" << std::get<1>(nodes_[i]) << std::endl;
    }
  }
}

// Print methods
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>::printNodesToRegion() const {
  if (0 == comm_->getRank()) {
    std::cout << "Total number of mesh nodes: " << num_total_nodes_ << std::endl;
    std::cout << "Total number of mesh regions: " << num_total_regions_ << std::endl;
    std::cout << "Number of rows in nodes_ structure: " << nodes_.size() << std::endl;
    for (int i = 0; i < nodesToRegion_.size(); ++i) {
      std::cout << "Node " << std::get<0>(nodesToRegion_[i]) << "\t belongs to regions: " << std::get<1>(nodesToRegion_[i]) << std::endl;
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void RegionHandler<Scalar, LocalOrdinal, GlobalOrdinal, Node>::printInactive() const {
  if (maps_.composite_map_.empty())
    std::cout << "INACTIVE PROC ID: " << comm_->getRank() << std::endl;
}

}  // namespace Xpetra

#endif
