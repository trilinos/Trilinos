// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TrilinosCouplings_Statistics_hpp
#define __TrilinosCouplings_Statistics_hpp

// Intrepid includes
#include <Intrepid_FunctionSpaceTools.hpp>
#include <Intrepid_CellTools.hpp>
#include <Intrepid_ArrayTools.hpp>
#include <Intrepid_RealSpaceTools.hpp>
#include <Intrepid_Utils.hpp>

// Xpetra
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_CrsGraph.hpp>

#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraCrsGraph.hpp>
#endif

#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_EpetraCrsGraph.hpp>
#endif


// Teuchos
#include <Teuchos_Comm.hpp>


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class MachineLearningStatistics_Hex3D {
  using ST = Scalar;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using crsgraph_type    = Xpetra::CrsGraph<LO, GO, Node>;
  using multivector_type = Xpetra::MultiVector<ST, LO, GO, Node>;
  using vector_type      = Xpetra::Vector<ST, LO, GO, Node>;
  
  public:

  MachineLearningStatistics_Hex3D(long long numElemsGlobal_):numElemsGlobal(numElemsGlobal_) {
    local_stat_max.resize(NUM_STATISTICS);
    local_stat_min.resize(NUM_STATISTICS);
    local_stat_sum.resize(NUM_STATISTICS);
    global_stat_max.resize(NUM_STATISTICS);
    global_stat_min.resize(NUM_STATISTICS);
    global_stat_sum.resize(NUM_STATISTICS);

    for(int i=0; i<NUM_STATISTICS; i++) {
      local_stat_max[i] = 0.0;
      local_stat_min[i] = std::numeric_limits<double>::max();
      local_stat_sum[i] = 0.0;
    }
  }

  double distance2(Intrepid::FieldContainer<Scalar> & coord, int n1, int n2) {
    double dist = 0.0;
    for(int i=0; i<coord.dimension(1); i++)
      dist += (coord(n2,i) -coord(n1,i)) * (coord(n2,i) -coord(n1,i));
    return sqrt(dist);
  }

  double myDistance2(const Xpetra::MultiVector<ST, LO, GO, Node> &v, int i0, int i1) {
    const size_t numVectors = v.getNumVectors();
    double distance = 0.0;
    for (size_t j=0; j<numVectors; j++) {
      distance += (v.getData(j)[i0]-v.getData(j)[i1])*(v.getData(j)[i0]-v.getData(j)[i1]);
    }
    return distance;
  }


  
  /**********************************************************************************/
  /****************************** STATISTICS (Part I) *******************************/
  /**********************************************************************************/
  // Statistics: Compute max / min of sigma parameter, mesh information  
  void Phase1(Intrepid::FieldContainer<int>&elemToNode,Intrepid::FieldContainer<int>&elemToEdge,Intrepid::FieldContainer<int>&edgeToNode,Intrepid::FieldContainer<double>&nodeCoord,Intrepid::FieldContainer<double>&sigmaVal) {
    double maxmin_ratio = 0;
    double stretch = 0;
    double diag_ratio = 0;
    double diag_length_max = 0;
    double diag_length_min = 0;
    double edge_length_max = 0;
    double edge_length_min = 0;
    double principle_axis_1[] = {0, 0, 0};
    double principle_axis_2[] = {0, 0, 0};
    double principle_axis_3[] = {0, 0, 0};
    double hex_0[] = {0, 0, 0};
    double hex_1[] = {0, 0, 0};
    double hex_2[] = {0, 0, 0};
    double hex_3[] = {0, 0, 0};
    double hex_4[] = {0, 0, 0};
    double hex_5[] = {0, 0, 0};
    double hex_6[] = {0, 0, 0};
    double hex_7[] = {0, 0, 0};
    
    double dist = 0;
    int diag_nodes1[] = {0, 1, 2, 3};
    int diag_nodes2[] = {6, 7, 4, 5};
    int edge_nodes_1[] = {0, 0, 0, 1, 1, 2, 2, 3, 4, 4, 5, 6};
    int edge_nodes_2[] = {1, 3, 4, 2, 5, 3, 6, 7, 5, 7, 6, 7};
    int edge_opposites_1[] = { 0,  1, 2, 3, 4, 5};
    int edge_opposites_2[] = {11, 10, 6, 9, 7, 8};
    double edge_length_1 = 0;
    double edge_length_2 = 0;
  
    int edge = elemToEdge(0, 0);
    int node1 = edgeToNode(edge, 0);
    int node2 = edgeToNode(edge, 1);
    int node3 = edgeToNode(edge, 0);
    int node4 = edgeToNode(edge, 1);

    int numElems = elemToNode.dimension(0);
    for(int i=0; i<numElems; i++) {
      // Set up hex nodes
      int hexnode0 = elemToNode(i, 0);
      int hexnode1 = elemToNode(i, 1);
      int hexnode2 = elemToNode(i, 2);
      int hexnode3 = elemToNode(i, 3);
      int hexnode4 = elemToNode(i, 4);
      int hexnode5 = elemToNode(i, 5);
      int hexnode6 = elemToNode(i, 6);
      int hexnode7 = elemToNode(i, 7);
      
      // Now I have to get the node coordinates
      for(int j=0; j<dim; j++) {
        hex_0[j] = nodeCoord(hexnode0, j);
        hex_1[j] = nodeCoord(hexnode1, j);
        hex_2[j] = nodeCoord(hexnode2, j);
        hex_3[j] = nodeCoord(hexnode3, j);
        hex_4[j] = nodeCoord(hexnode4, j);
        hex_5[j] = nodeCoord(hexnode5, j);
        hex_6[j] = nodeCoord(hexnode6, j);
        hex_7[j] = nodeCoord(hexnode7, j);
      }
      
      double pr1_norm = 0;
      double pr2_norm = 0;
      double pr3_norm = 0;
      for(int j=0; j<dim; j++){
        principle_axis_1[j] = hex_1[j] - hex_0[j] + hex_2[j] - hex_3[j] + hex_5[j] - hex_4[j] + hex_6[j] - hex_7[j];
        pr1_norm += principle_axis_1[j] * principle_axis_1[j];
        principle_axis_2[j] = hex_3[j] - hex_0[j] + hex_2[j] - hex_1[j] + hex_7[j] - hex_4[j] + hex_6[j] - hex_5[j];
        pr2_norm += principle_axis_2[j] * principle_axis_2[j];
        principle_axis_3[j] = hex_4[j] - hex_0[j] + hex_5[j] - hex_1[j] + hex_6[j] - hex_2[j] + hex_7[j] - hex_3[j];
        pr3_norm += principle_axis_3[j] * principle_axis_3[j];
      }
      
      pr1_norm = sqrt(pr1_norm);
      pr2_norm = sqrt(pr2_norm);
      pr3_norm = sqrt(pr3_norm);
      
      for(int j=0; j<dim; j++) {
        principle_axis_1[j] = principle_axis_1[j] / pr1_norm;
        principle_axis_2[j] = principle_axis_2[j] / pr2_norm;
        principle_axis_3[j] = principle_axis_3[j] / pr3_norm;
      }
      
      
      // 0 - Material property
      local_stat_max[0] = std::max(local_stat_max[0],sigmaVal(i));
      local_stat_min[0] = std::min(local_stat_min[0],sigmaVal(i));
      local_stat_sum[0] += sigmaVal(i);
      
      edge = elemToEdge(i, 0);
      node1 = edgeToNode(edge, 0);
      node2 = edgeToNode(edge, 1);
      
      // 1 - Max/min edge - ratio of max to min edge length
      edge_length_max = distance2(nodeCoord, node1, node2);
      edge_length_min = edge_length_max;
      int numEdgesPerElem = elemToEdge.dimension(1);
      for (int j=0; j<numEdgesPerElem; j++) {
        edge = elemToEdge(i,j);
        node1 = edgeToNode(edge,0);
        node2 = edgeToNode(edge,1);
        dist = distance2(nodeCoord,node1,node2);
        edge_length_max = std::max(edge_length_max,dist);
        edge_length_min = std::min(edge_length_min,dist);
      }
      maxmin_ratio = edge_length_max / edge_length_min;
      local_stat_max[1] = std::max(local_stat_max[1],maxmin_ratio);
      local_stat_min[1] = std::min(local_stat_min[1],maxmin_ratio);
      local_stat_sum[1] += maxmin_ratio;
      
      // 2 - det of cell Jacobian (later)
      
      // 3 - Stretch
      diag_length_max = distance2(nodeCoord, elemToNode(i, diag_nodes1[0]),
                                  elemToNode(i, diag_nodes2[0]));
      diag_length_min = distance2(nodeCoord, elemToNode(i, diag_nodes1[0]),
                                  elemToNode(i, diag_nodes2[0]));
      for (int j=0; j<NUM_NODE_PAIRS; j++) {
        node1 = elemToNode(i, diag_nodes1[j]);
        node2 = elemToNode(i, diag_nodes2[j]);
        dist = distance2(nodeCoord, node1, node2);
        diag_length_max = std::max(diag_length_max, dist);
        diag_length_min = std::min(diag_length_min, dist);
      }
      stretch = sqrt(3) * edge_length_min / diag_length_max;
      diag_ratio = diag_length_min / diag_length_max;
      local_stat_max[3] = std::max(local_stat_max[3], stretch);
      local_stat_min[3] = std::min(local_stat_min[3], stretch);
      local_stat_sum[3] += stretch;
      
      // 4 - Diagonal Ratio
      local_stat_max[4] = std::max(local_stat_max[4], diag_ratio);
      local_stat_min[4] = std::min(local_stat_min[4], diag_ratio);
      local_stat_sum[4] += diag_ratio;
      
      // 5 - Inverse Taper
      node1 = elemToNode(i, diag_nodes1[0]);
      node2 = elemToNode(i, diag_nodes1[1]);
      node3 = elemToNode(i, diag_nodes2[0]);
      node4 = elemToNode(i, diag_nodes2[1]);
      edge_length_1 = distance2(nodeCoord, node1, node2);
      edge_length_2 = distance2(nodeCoord, node3, node4);
      double ratio = edge_length_1 / edge_length_2;
      ratio = std::min(ratio, 1/ratio);
      for (int j=0; j<NUM_EDGE_PAIRS; j++) {
        node1 = elemToNode(i, edge_nodes_1[edge_opposites_1[j]]);
        node2 = elemToNode(i, edge_nodes_2[edge_opposites_1[j]]);
        node3 = elemToNode(i, edge_nodes_1[edge_opposites_2[j]]);
        node4 = elemToNode(i, edge_nodes_2[edge_opposites_2[j]]);
        edge_length_1 = distance2(nodeCoord, node1, node2);
        edge_length_2 = distance2(nodeCoord, node3, node4);
        double my_ratio = edge_length_1 / edge_length_2;
        my_ratio = std::min(my_ratio, 1/my_ratio);
        ratio = std::min(ratio, my_ratio);
      }
      local_stat_max[5] = std::max(local_stat_max[5], ratio);
      local_stat_min[5] = std::min(local_stat_min[5], ratio);
      local_stat_sum[5] += ratio;
      
      // 6 - Skew
      double skew = 0;
      double skew1 = 0;
      double skew2 = 0;
      double skew3 = 0;
      
      for(int j=0; j<dim; j++) {
        skew1 += principle_axis_1[j] * principle_axis_2[j];
        skew2 += principle_axis_1[j] * principle_axis_3[j];
        skew3 += principle_axis_1[j] * principle_axis_3[j];
      }
      skew1 = std::abs(skew1);
      skew2 = std::abs(skew2);
      skew3 = std::abs(skew3);
      skew = std::max(skew1, skew2);
      skew = std::max(skew, skew3);
      local_stat_max[6] = std::max(local_stat_max[6], skew);
      local_stat_min[6] = std::min(local_stat_min[6], skew);
      local_stat_sum[6] += skew;
    }
     
  }

  /**********************************************************************************/
  /***************************** STATISTICS (Part IIa) ******************************/
  /**********************************************************************************/
  void Phase2a(Intrepid::FieldContainer<Scalar> &worksetJacobDet,Intrepid::FieldContainer<Scalar> &worksetCubWeights) {
    int worksetSize  = worksetJacobDet.dimension(0);
    int numCubPoints = worksetJacobDet.dimension(1);

    bool weightsWorkset = (worksetCubWeights.rank()==2)?true:false;

    for(int i=0; i<worksetSize; i++) {
      // 0 - Material property
      // 1 - Max/min edge - ratio of max to min edge length
      // 2 - det of cell Jacobian (later)
      double elementdetJ = 0.0, elementWeight=0.0;
      for(int j=0; j<numCubPoints; j++) {
        double weight = weightsWorkset ? worksetCubWeights(i,j) : worksetCubWeights(j);
        elementdetJ   += worksetJacobDet(i,j) * weight;
        elementWeight += weight;
      }
      double detJ = elementdetJ / elementWeight;
      local_stat_max[2] = std::max(local_stat_max[2],detJ);
      local_stat_min[2] = std::min(local_stat_min[2],detJ);
      local_stat_sum[2] += detJ;
    }
  }


  /**********************************************************************************/
  /***************************** STATISTICS (Part IIb) ******************************/
  /**********************************************************************************/

  void Phase2b(Teuchos::RCP<const Xpetra::CrsGraph<LO, GO, Node> > gl_StiffGraph, Teuchos::RCP<Xpetra::MultiVector<ST, LO, GO,Node> > coords) {
    using multivector_factory = Xpetra::MultiVectorFactory<ST,LO,GO,Node>;
    using vector_factory      = Xpetra::VectorFactory<ST,LO,GO,Node>;

    Teuchos::RCP<multivector_type> coordsOwnedPlusShared;
    comm = gl_StiffGraph->getRowMap()->getComm();

    if (!(gl_StiffGraph->getImporter().is_null())) {
      coordsOwnedPlusShared = multivector_factory::Build(gl_StiffGraph->getColMap(), 3, true);
      coordsOwnedPlusShared->doImport(*coords, *gl_StiffGraph->getImporter(), Xpetra::CombineMode::ADD);
    }
    else {
      coordsOwnedPlusShared = coords;
    }

    Teuchos::RCP<const Xpetra::Map<LO, GO, Node> > rowMap = gl_StiffGraph->getRowMap();
    Teuchos::RCP<vector_type > laplDiagOwned = vector_factory::Build(rowMap, true);
    Teuchos::ArrayView<const LO> indices;
    size_t numOwnedRows = rowMap->getLocalNumElements();
    for (size_t row=0; row<numOwnedRows; row++) {
      gl_StiffGraph->getLocalRowView(row, indices);
      size_t numIndices = indices.size();
      for (size_t j=0; j<numIndices; j++) {
        size_t col = indices[j];
        if (row == col) continue;
        laplDiagOwned->sumIntoLocalValue(row, 1/myDistance2(*coordsOwnedPlusShared, row, col));
      }
    }
    Teuchos::RCP<vector_type> laplDiagOwnedPlusShared;
    if (!gl_StiffGraph->getImporter().is_null()) {
      laplDiagOwnedPlusShared = vector_factory::Build(gl_StiffGraph->getColMap(), true);
      laplDiagOwnedPlusShared->doImport(*laplDiagOwned, *gl_StiffGraph->getImporter(), Xpetra::CombineMode::ADD);
    }
    else {
      laplDiagOwnedPlusShared = laplDiagOwned;
    }

    for (size_t row=0; row<numOwnedRows; row++) {
      gl_StiffGraph->getLocalRowView(row, indices);
      size_t numIndices = indices.size();
      for(size_t j=0; j<numIndices; j++) {
        size_t col = indices[j];
        if (row==col) continue;
        double laplVal = 1.0 / myDistance2(*coordsOwnedPlusShared, row, col);
        double aiiajj = std::abs(laplDiagOwnedPlusShared->getData(0)[row]*laplDiagOwnedPlusShared->getData(0)[col]);
        double aij = laplVal * laplVal;
        double ratio = sqrt(aij / aiiajj);
        local_stat_max[7] = std::max(local_stat_max[7], ratio);
        local_stat_min[7] = std::min(local_stat_min[7], ratio);
        local_stat_sum[7] += ratio;
      }
    }
    globalNumMatrixEntries = gl_StiffGraph->getGlobalNumEntries();
  }

#ifdef HAVE_XPETRA_TPETRA
  void Phase2b(Teuchos::RCP<const Tpetra::CrsGraph<LO, GO, Node> > gl_StiffGraph, Teuchos::RCP<Tpetra::MultiVector<ST, LO, GO,Node> > coords) {
    Teuchos::RCP<multivector_type> coords_X = Teuchos::rcp(new Xpetra::TpetraMultiVector<ST,LO,GO,Node>(coords));
    Teuchos::RCP<const crsgraph_type> graph_X = Teuchos::rcp(new Xpetra::TpetraCrsGraph<LO,GO,Node>(Teuchos::rcp_const_cast<Tpetra::CrsGraph<LO,GO,Node> >(gl_StiffGraph)));

    Phase2b(graph_X, coords_X);
  }
#endif

#ifdef HAVE_XPETRA_EPETRA
  void Phase2b(Teuchos::RCP<const Epetra_CrsGraph> gl_StiffGraph, Teuchos::RCP<Epetra_MultiVector> coords) {
    Teuchos::RCP<multivector_type> coords_X = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GO,Node>(coords));
    Teuchos::RCP<const crsgraph_type> graph_X = Teuchos::rcp(new Xpetra::EpetraCrsGraphT<GO,Node>(Teuchos::rcp_const_cast<Epetra_CrsGraph>(gl_StiffGraph)));

    Phase2b(graph_X, coords_X);
  }
#endif


  /**********************************************************************************/
  /***************************** STATISTICS (Part III) ******************************/
  /**********************************************************************************/
  void Phase3() {
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MIN,NUM_STATISTICS,local_stat_min.data(),global_stat_min.data());
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,NUM_STATISTICS,local_stat_max.data(),global_stat_max.data());
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM,NUM_STATISTICS,local_stat_sum.data(),global_stat_sum.data());
    // NOTE: All output properties should be unitless if we want to compare across problems.
    // NOTE: Should the mean be weighted by cell volume?  That is not currently done.
    
    // 0 - Material property
    problemStatistics.set("sigma: min/mean",global_stat_min[0]/global_stat_sum[0]*numElemsGlobal);
    problemStatistics.set("sigma: max/mean",global_stat_max[0]/global_stat_sum[0]*numElemsGlobal);
    
    // 1 - Max/min edge ratio
    problemStatistics.set("element edge ratio: min",global_stat_min[1]);
    problemStatistics.set("element edge ratio: max",global_stat_max[1]);
    problemStatistics.set("element edge ratio: mean",global_stat_sum[1] / numElemsGlobal);
    
    // 2 - det of cell Jacobian (later)
    problemStatistics.set("element det jacobian: min/mean",global_stat_min[2]/global_stat_sum[2]*numElemsGlobal);
    problemStatistics.set("element det jacobian: max/mean",global_stat_max[2]/global_stat_sum[2]*numElemsGlobal);
    
    // 3 - Stretch
    problemStatistics.set("Stretch max", global_stat_max[3]);
    problemStatistics.set("Stretch min", global_stat_min[3]);
    problemStatistics.set("Stretch mean", global_stat_sum[3] / numElemsGlobal);
    
    // 4 - Diagonal Ratio
    problemStatistics.set("Diagonal Ratio max", global_stat_max[4]);
    problemStatistics.set("Diagonal Ratio min", global_stat_min[4]);
    problemStatistics.set("Diagonal Ratio mean", global_stat_sum[4]/numElemsGlobal);
    
    // 5 - Inverse Taper
    problemStatistics.set("Inverse Taper max", global_stat_max[5]);
    problemStatistics.set("Inverse Taper min", global_stat_min[5]);
    problemStatistics.set("Inverse Taper mean", global_stat_sum[5] / numElemsGlobal);
    
    // 6 - Skew
    problemStatistics.set("Skew max", global_stat_max[6]);
    problemStatistics.set("Skew min", global_stat_min[6]);
    problemStatistics.set("Skew mean", global_stat_sum[6] / numElemsGlobal);
    
    // 7 - Lapl Diag
    problemStatistics.set("Lapl Diag max", global_stat_max[7]);
    problemStatistics.set("Lapl Diag min", global_stat_min[7]);
    problemStatistics.set("Lapl Diag mean", global_stat_sum[7] / globalNumMatrixEntries);
  }  
  

  Teuchos::ParameterList GetStatistics() {return problemStatistics;}

private:

  // Internal
  long long numElemsGlobal;
  long long globalNumMatrixEntries;
  const int NUM_STATISTICS = 8;
  std::vector<double> local_stat_max;
  std::vector<double> local_stat_min;
  std::vector<double> local_stat_sum;
  std::vector<double> global_stat_max;
  std::vector<double> global_stat_min;
  std::vector<double> global_stat_sum;

  const int NUM_NODE_PAIRS = 4;
  const int NUM_EDGE_PAIRS = 6;
  const int dim = 3;

  
  Teuchos::RCP<Tpetra::CrsGraph<LO, GO, Node> > crs_graph;
  Teuchos::RCP<const Teuchos::Comm<int> > comm;
  Teuchos::ParameterList problemStatistics;

};

#endif
