// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef MUELU_CLASSICALMAPFACTORY_DEF_HPP_
#define MUELU_CLASSICALMAPFACTORY_DEF_HPP_

#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Xpetra_Vector.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_ClassicalMapFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_GraphBase.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_LWGraph.hpp"


// NOTE: We should be checking for KokkosKernels here, but
// MueLu doesn't have a macro for that
#ifdef HAVE_MUELU_KOKKOSCORE
#include "MueLu_LWGraph_kokkos.hpp"
#include <KokkosGraph_Distance1ColorHandle.hpp>
#include <KokkosGraph_Distance1Color.hpp>
#endif

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const
  {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregation: deterministic");
    SET_VALID_ENTRY("aggregation: coloring algorithm");
#undef SET_VALID_ENTRY
    validParamList->set< RCP<const FactoryBase> >("A",              Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
    validParamList->set< RCP<const FactoryBase> >("Graph",       null, "Generating factory of the graph");
    validParamList->set< RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");
 
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const
  {
    Input(currentLevel, "A");
    Input(currentLevel, "UnAmalgamationInfo");
    Input(currentLevel, "Graph");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const
  {
    FactoryMonitor m(*this, "Build", currentLevel);
    RCP<const GraphBase> graph = Get<RCP<GraphBase> >(currentLevel,"Graph");
    RCP<const Matrix> A = Get<RCP<Matrix> >(currentLevel,"A");
    const ParameterList& pL = GetParameterList();
    /* ============================================================= */
    /* Phase 1 : Compute an initial MIS                              */
    /* ============================================================= */
    ArrayRCP<LO> myColors;
    LO numColors=0;
    // FIXME:  This is not going to respect coloring at or near processor
    // boundaries, so that'll need to get cleaned up later
    RCP<LocalOrdinalVector> fc_splitting;
    std::string coloringAlgo = pL.get<std::string>("aggregation: coloring algorithm");
    if(coloringAlgo == "file") {
      // Read the CF splitting from disk
      // NOTE: For interoperability reasons, this is dependent on the point_type enum not changing
      std::string map_file   = std::string("map_fcsplitting_") + std::to_string(currentLevel.GetLevelID()) + std::string(".m");
      std::string color_file = std::string("fcsplitting_")     + std::to_string(currentLevel.GetLevelID()) + std::string(".m");
        
      FILE * mapfile = fopen(map_file.c_str(),"r");
      using real_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;
      using RealValuedMultiVector = typename Xpetra::MultiVector<real_type,LO,GO,NO>;
      RCP<RealValuedMultiVector> mv;

      if(mapfile) {
        fclose(mapfile);
        RCP<const Map> colorMap = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadMap(map_file, A->getRowMap()->lib(), A->getRowMap()->getComm());
        TEUCHOS_TEST_FOR_EXCEPTION(!colorMap->isCompatible(*A->getRowMap()),std::invalid_argument,"Coloring on disk has incompatible map with A");

        mv = Xpetra::IO<real_type, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(color_file,colorMap);
      }
      else {
        // Use A's rowmap and hope it matches
        mv = Xpetra::IO<real_type, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVector(color_file,A->getRowMap());
      }
      TEUCHOS_TEST_FOR_EXCEPTION(mv.is_null(),std::invalid_argument,"Coloring on disk cannot be read");      
      fc_splitting = LocalOrdinalVectorFactory::Build(A->getRowMap());
      TEUCHOS_TEST_FOR_EXCEPTION(mv->getLocalLength() != fc_splitting->getLocalLength(),std::invalid_argument,"Coloring map mismatch");
      ArrayRCP<const real_type> mv_data= mv->getData(0);
      ArrayRCP<LO> fc_data= fc_splitting->getDataNonConst(0);
      
      for(LO i=0; i<(LO)fc_data.size(); i++)
        fc_data[i] = (LO) mv_data[i];

    }
#ifdef HAVE_MUELU_KOKKOSCORE  
    else if(coloringAlgo != "MIS" && graph->GetDomainMap()->lib() == Xpetra::UseTpetra) {
      SubFactoryMonitor sfm(*this,"GraphColoring",currentLevel);
      DoGraphColoring(*graph,myColors,numColors);
    }
#endif
    else if(coloringAlgo == "MIS") { 
      SubFactoryMonitor sfm(*this,"MIS",currentLevel);
      DoMISNaive(*graph,myColors,numColors);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Unrecognized distance 1 coloring algorithm");
    }

    // FIXME: This coloring will either need to be done MPI parallel, or
    // there needs to be a cleanup phase to fix mistakes

    /* ============================================================= */
    /* Phase 2 : Mark the C-Points                                   */
    /* ============================================================= */
    LO num_c_points = 0, num_d_points=0, num_f_points = 0;
    if(fc_splitting.is_null()) {
      // We just have a coloring, so we need to generate a splitting
      auto boundaryNodes = graph->GetBoundaryNodeMap();
      fc_splitting = LocalOrdinalVectorFactory::Build(A->getRowMap());
      ArrayRCP<LO> myPointType = fc_splitting->getDataNonConst(0);
      for(LO i=0; i<(LO)myColors.size(); i++) {
        if(boundaryNodes[i]) {
          myPointType[i] = DIRICHLET_PT;
          num_d_points++;
        }
        else if ((LO)myColors[i] == 1) {
        myPointType[i] = C_PT;
        num_c_points++;
        }
        else
          myPointType[i] = F_PT;
      }
      num_f_points = (LO)myColors.size() - num_d_points - num_c_points;
    }
    else {
      // If we read the splitting off disk, we just need to count
      ArrayRCP<LO> myPointType = fc_splitting->getDataNonConst(0);

      for(LO i=0; i<(LO)myPointType.size(); i++) {
        if(myPointType[i] == DIRICHLET_PT)
          num_d_points++;
        else if (myPointType[i] == C_PT)
          num_c_points++;
      }
      num_f_points = (LO)myPointType.size() - num_d_points - num_c_points;
    }

    // FIXME: This array will need to be ghosted so we can get the point_types
    // of the neighbors

    // FIXME:  These stats will need to be reduced
    GetOStream(Statistics1) << "ClassicalMapFactory: C/F/D = "<<num_c_points<<"/"<<num_f_points<<"/"<<num_d_points<<std::endl;

    /* Generate the Coarse map */
    RCP<const Map> coarseMap;
    {
      SubFactoryMonitor sfm(*this,"Coarse Map",currentLevel);
      GenerateCoarseMap(*A->getRowMap(),num_c_points,coarseMap);
    }
   
    Set(currentLevel, "FC Splitting",fc_splitting);
    Set(currentLevel, "CoarseMap", coarseMap);
   
  }

/* ************************************************************************* */
template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
GenerateCoarseMap(const Map & fineMap, LO num_c_points, RCP<const Map> & coarseMap) const {

  // FIXME: Assumes scalar PDE
  std::vector<size_t> stridingInfo_(1);
  stridingInfo_[0]=1;
  GO domainGIDOffset = 0;

  coarseMap = StridedMapFactory::Build(fineMap.lib(),
                                       Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                       num_c_points,
                                       fineMap.getIndexBase(),
                                       stridingInfo_,
                                       fineMap.getComm(),
                                       domainGIDOffset); 
}


/* ************************************************************************* */
template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
DoGraphColoring(const GraphBase & graph, ArrayRCP<LO> & myColors_out, LO & numColors) const {
  const ParameterList& pL = GetParameterList();
#ifdef HAVE_MUELU_KOKKOSCORE  
  using graph_t = typename LWGraph_kokkos::local_graph_type;
  using KernelHandle = KokkosKernels::Experimental::
    KokkosKernelsHandle<typename graph_t::row_map_type::value_type,
                        typename graph_t::entries_type::value_type,
                        typename graph_t::entries_type::value_type,
                        typename graph_t::device_type::execution_space,
                        typename graph_t::device_type::memory_space,
                        typename graph_t::device_type::memory_space>;
  KernelHandle kh;

  // Leave gc algorithm choice as the default
  kh.create_graph_coloring_handle();
  
  // Get the distance-1 graph coloring handle
  auto coloringHandle = kh.get_graph_coloring_handle();      
  
  // Set the distance-1 coloring algorithm to use
  if(pL.get<bool>("aggregation: deterministic") == true) {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_SERIAL );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: serial" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "serial") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_SERIAL );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: serial" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "vertex based") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_VB );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "vertex based bit array") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_VBBIT );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based bit array" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "vertex based color set") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_VBCS );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based color set" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "vertex based deterministic") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_VBD );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based deterministic" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "vertex based deterministic bit array") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_VBDBIT );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: vertex based deterministic bit array" << std::endl;
  } else if(pL.get<std::string>("aggregation: coloring algorithm") == "edge based") {
    coloringHandle->set_algorithm( KokkosGraph::COLORING_EB );
    if(IsPrint(Statistics1)) GetOStream(Statistics1) << "  algorithm: edge based" << std::endl;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Unrecognized distance 1 coloring algorithm");
  }
  
  // Create device views for graph rowptrs/colinds
  size_t numRows = graph.GetNodeNumVertices();
  auto graphLWK = dynamic_cast<const LWGraph_kokkos*>(&graph);
  auto graphLW  = dynamic_cast<const LWGraph*>(&graph);
  auto graphG   = dynamic_cast<const Graph*>(&graph);
  TEUCHOS_TEST_FOR_EXCEPTION(!graphLW && !graphLWK && !graphG,std::invalid_argument,"Graph is not a LWGraph or LWGraph_kokkos object");
    // Run d1 graph coloring
    // Assume that the graph is symmetric so row map/entries and col map/entries are the same

  if(graphLWK) {
    KokkosGraph::Experimental::graph_color(&kh, 
                                           numRows, 
                                           numRows, // FIXME: This should be the number of columns
                                           graphLWK->getRowPtrs(),
                                           graphLWK->getEntries(),
                                         true);
  }
  else if(graphLW) {
    auto rowptrs = graphLW->getRowPtrs();
    auto entries = graphLW->getEntries();
    // Copy rowptrs to a size_t, because kokkos-kernels doesn't like rowptrs as LO's
    Teuchos::Array<size_t> rowptrs_s(rowptrs.size());
    std::copy(rowptrs.begin(),rowptrs.end(),rowptrs_s.begin());
    Kokkos::View<const size_t*,Kokkos::LayoutLeft,Kokkos::HostSpace> rowptrs_v(rowptrs_s.data(),(size_t)rowptrs.size());
    Kokkos::View<const LO*,Kokkos::LayoutLeft,Kokkos::HostSpace> entries_v(entries.getRawPtr(),(size_t)entries.size());
    KokkosGraph::Experimental::graph_color(&kh, 
                                           numRows, 
                                           numRows, // FIXME: This should be the number of columns
                                           rowptrs_v,
                                           entries_v,
                                           true);    
  }
  else if(graphG) {  
    // FIXME:  This is a terrible, terrible hack, based on 0-based local indexing.
    RCP<const CrsGraph> graphC = graphG->GetGraph();
    size_t numEntries = graphC->getNodeNumEntries();
    ArrayView<const LO> indices;
    graphC->getLocalRowView(0,indices);
    Kokkos::View<size_t*,Kokkos::LayoutLeft,Kokkos::HostSpace> rowptrs_v("rowptrs_v",graphC->getNodeNumRows()+1);
    rowptrs_v[0]=0;
    for(LO i=0; i<(LO)graphC->getNodeNumRows()+1; i++) 
      rowptrs_v[i+1] = rowptrs_v[i] + graphC->getNumEntriesInLocalRow(i);
    Kokkos::View<const LO*,Kokkos::LayoutLeft,Kokkos::HostSpace> entries_v(&indices[0],numEntries);    
    KokkosGraph::Experimental::graph_color(&kh, 
                                           numRows, 
                                           numRows, // FIXME: This should be the number of columns
                                           rowptrs_v,
                                           entries_v,
                                           true);       
  }

  
  // Extract the colors and store them in the aggregates
  auto myColors_d = coloringHandle->get_vertex_colors();
  numColors = static_cast<LO>(coloringHandle->get_num_colors());

  // Copy back to host
  auto myColors_h = Kokkos::create_mirror_view(myColors_d);
  myColors_out.resize(myColors_h.size());
  Kokkos::View<LO*,Kokkos::LayoutLeft,Kokkos::HostSpace> myColors_v(&myColors_out[0],myColors_h.size());
  Kokkos::deep_copy(myColors_v,myColors_h);
  
  //clean up coloring handle
  kh.destroy_graph_coloring_handle();
#else
  TEUCHOS_TEST_FOR_EXCEPTION(1, Exceptions::RuntimeError,"ClassicalMapFactory: Requires KokkosKernels");
#endif
  
}// end DoGraphColoring
    

/* ************************************************************************* */
template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalMapFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
DoMISNaive(const GraphBase & graph, ArrayRCP<LO> & myColors, LO & numColors) const {
  // This is a fall-back routine for when we don't have Kokkos or when it isn't initialized
  // We just do greedy MIS because this is easy to write.

  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  LO MIS = Teuchos::ScalarTraits<LO>::one();

  //FIXME: Not efficient
  myColors.resize(0);
  myColors.resize(graph.GetNodeNumVertices(),LO_INVALID);
  auto boundaryNodes = graph.GetBoundaryNodeMap();
  LO Nrows = (LO)graph.GetNodeNumVertices();

  
  for(LO row=0; row < Nrows; row++) {
    if(boundaryNodes[row]) 
      continue;
    ArrayView<const LO> indices = graph.getNeighborVertices(row);
    bool has_colored_neighbor=false;
    for(LO j=0; !has_colored_neighbor && j<(LO)indices.size(); j++) {
      // FIXME: This does not handle ghosting correctly
      if(myColors[indices[j]] == MIS) 
        has_colored_neighbor=true;
    }
    if(!has_colored_neighbor)
      myColors[row] = MIS;   
  } 
  numColors=1;
}


} //namespace MueLu

#endif /* MUELU_CLASSICALMAPFACTORY_DEF_HPP_ */
