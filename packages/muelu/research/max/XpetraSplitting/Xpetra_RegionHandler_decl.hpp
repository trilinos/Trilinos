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
 * Xpetra_RegionHandler_decl.hpp
 *
 * Created on: August 17, 2017
 * 	Author: Massimiliano Lupo Pasini (massimiliano.lupo.pasini@gmail.com)
 */
#ifndef XPETRA_REGIONHANDLER_DECL_HPP
#define XPETRA_REGIONHANDLER_DECL_HPP

#include "Xpetra_Map.hpp"
#include <string>
#include <fstream>
#include <vector>

namespace Xpetra {

template <class GlobalOrdinal>
bool compareRegions(const std::tuple<GlobalOrdinal, GlobalOrdinal> &, const std::tuple<GlobalOrdinal, GlobalOrdinal> &);

template <class GlobalOrdinal>
bool compareNodes(const std::tuple<GlobalOrdinal, GlobalOrdinal> &, const std::tuple<GlobalOrdinal, GlobalOrdinal> &);

// Definition of the predicate for the node_ structure.
// Given a tuple made of node index and a specific region it belongs to,
// this predicate returns true if the node belongs to the region specified in input to the predicate.
template <class GlobalOrdinal>
class checkerNode {
 public:
  // Constructor
  checkerNode(GlobalOrdinal region_index) { region_index_ = region_index; };

  // Unary Operator
  bool operator()(const std::tuple<GlobalOrdinal, GlobalOrdinal> &node) { return (std::get<1>(node) == region_index_); }

 private:
  GlobalOrdinal region_index_;
};

// Definition of the predicate for the nodesToRegion_ sitructure
// Given a tuple made of node index and a vector with labels of regions it belongs to,
// this predicate returns true if the node coincides with the node specified in input to the predicate.
template <class GlobalOrdinal>
class checkerNodesToRegion {
 public:
  // Constructor
  checkerNodesToRegion(GlobalOrdinal node_index) { node_index_ = node_index; };

  // Unary Operator
  bool operator()(const std::tuple<GlobalOrdinal, Array<GlobalOrdinal> > &node) { return (std::get<0>(node) == node_index_); }

 private:
  GlobalOrdinal node_index_;
};

// This is an auxiliary class to store row maps for the composite matrix, region matrices and
// a regionToAll map to link region node indices with the composite ones
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class Splitting_MapsInfo {
 public:
  Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > regionToAll_;  // used as a map for a RegionToAll node index
  Array<GlobalOrdinal> composite_map_;                                    // used as RowMap for composite matrices
  Array<Array<GlobalOrdinal> > region_maps_;                              // used as RowMap for region matrices
};

// This is the actual class that defines the regionHandler
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class RegionHandler {
 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor specifying the file name containing region information.
  RegionHandler(const std::string &, RCP<const Teuchos::Comm<int> >);

  //}
  //! @Interface methods
  //@{
  GlobalOrdinal GetNumGlobalElements() const { return num_total_nodes_; };
  GlobalOrdinal GetNumTotalRegions() const { return num_total_regions_; };
  GlobalOrdinal GetNumRegionNodes(GlobalOrdinal region_idx) const { return num_region_nodes_[region_idx]; };
  Array<GlobalOrdinal> GetGlobalRowMap() const { return maps_.composite_map_; };
  Array<GlobalOrdinal> GetRegionRowMap(GlobalOrdinal region_index) const;
  Array<Array<GlobalOrdinal> > GetRegionRowMaps() const { return maps_.region_maps_; };
  Array<Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > > GetRegionToAll() const;       // used as a map for a RegionToAll node index
  Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > GetRegionToAll(GlobalOrdinal) const;  // used as a map for a RegionToAll node index
  Array<std::tuple<int, Array<GlobalOrdinal> > > GetInterfaceNodes() const { return interfaceNodes_; };
  //}
  //! @Printout methods
  void printView() const;
  void printNodesToRegion() const;
  void printInactive() const;
  //}

 private:
  //! @Private variables
  //@{

  RCP<const Teuchos::Comm<int> > comm_;
  bool nodes_sorted_by_regions_ = false;

  // Global information
  GlobalOrdinal num_total_nodes_   = 0;
  GlobalOrdinal num_total_regions_ = 0;
  Array<std::tuple<GlobalOrdinal, GlobalOrdinal> > nodes_;  // basic structure that imports the information from the input file

  // the following two Array are used to handle the situation where either the number of processes exceeds the number of regions or viceversa
  Array<GlobalOrdinal> regions_per_proc_;                            // if num_proc > num_regions, then it says how many regions are owned by a single process, empty otherwise
  Array<std::tuple<int, Array<GlobalOrdinal> > > procs_per_region_;  // lists of processes instantiated for each region

  Array<std::tuple<int, Array<GlobalOrdinal> > > nodesToRegion_;   // for each node it lists the regions it belongs to
  Array<std::tuple<int, Array<GlobalOrdinal> > > interfaceNodes_;  // for each node on the interface it lists the regions it belongs to
  // vector which contains the number of region nodes for each domain region
  Array<GlobalOrdinal> num_region_nodes_;

  // Maps used for composite and region operators
  Splitting_MapsInfo<Scalar, LocalOrdinal, GlobalOrdinal, Node> maps_;
  //@}

  //! @Private Methods
  //@{
  void ReadFileInfo(const std::string &);
  void ComputeProcRegions();
  void NodesToRegion();
  void CreateRowMaps();
  //@}

};  // class RegionHandler

// This compare class is used to run the sorting algorithm on the list of nodes with associated regions they belong to.
// First, nodes are sorted in ascending order for region labels. Then, the sorting shuffles the nodes in ascending node index for
// each given region
template <class GlobalOrdinal>
bool compareRegions(const std::tuple<GlobalOrdinal, GlobalOrdinal> &lhs, const std::tuple<GlobalOrdinal, GlobalOrdinal> &rhs) {
  // First we prioritize the sorting according to the region label
  // If the region is the same, then the sorting looks at the composite node index
  if (std::get<1>(lhs) < std::get<1>(rhs))
    return true;
  else if (std::get<1>(lhs) == std::get<1>(rhs))
    return std::get<0>(lhs) < std::get<0>(rhs);
  else
    return false;
}

// This compare is sed to run the sorting algorithm where the nodes are ordered in ascendin order for thei node indes, regardless of the
// associated region index
template <class GlobalOrdinal>
bool compareNodes(const std::tuple<GlobalOrdinal, GlobalOrdinal> &lhs, const std::tuple<GlobalOrdinal, GlobalOrdinal> &rhs) {
  return std::get<0>(lhs) < std::get<0>(rhs);
}

}  // namespace Xpetra

#endif
