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
/*
 * MueLu_AggregationExportFactory_def.hpp
 *
 *  Created on: Feb 10, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_
#define MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include "MueLu_AggregationExportFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Monitor.hpp"
#include <vector>
#include <list>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <cstdio>
#include <cmath>

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

  std::string output_msg =
      "Output filename template (%TIMESTEP is replaced by \'Output file: time step\' variable,"
      "%ITER is replaced by \'Output file: iter\' variable, %LEVELID is replaced level id, %PROCID is replaced by processor id)";
  std::string output_def = "aggs_level%LEVELID_proc%PROCID.out";

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Factory for A.");
  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory for Coordinates.");
  validParamList->set<RCP<const FactoryBase> >("Graph", Teuchos::null, "Factory for Graph.");
  validParamList->set<RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Factory for Aggregates.");
  validParamList->set<RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Factory for UnAmalgamationInfo.");
  validParamList->set<RCP<const FactoryBase> >("DofsPerNode", Teuchos::null, "Factory for DofsPerNode.");
  // CMS/BMK: Old style factory-only options.  Deprecate me.
  validParamList->set<std::string>("Output filename", output_def, output_msg);
  validParamList->set<int>("Output file: time step", 0, "time step variable for output file name");
  validParamList->set<int>("Output file: iter", 0, "nonlinear iteration variable for output file name");

  // New-style master list options (here are same defaults as in masterList.xml)
  validParamList->set<std::string>("aggregation: output filename", "", "filename for VTK-style visualization output");
  validParamList->set<int>("aggregation: output file: time step", 0, "time step variable for output file name");       // Remove me?
  validParamList->set<int>("aggregation: output file: iter", 0, "nonlinear iteration variable for output file name");  // Remove me?
  validParamList->set<std::string>("aggregation: output file: agg style", "Point Cloud", "style of aggregate visualization for VTK output");
  validParamList->set<bool>("aggregation: output file: fine graph edges", false, "Whether to draw all fine node connections along with the aggregates.");
  validParamList->set<bool>("aggregation: output file: coarse graph edges", false, "Whether to draw all coarse node connections along with the aggregates.");
  validParamList->set<bool>("aggregation: output file: build colormap", false, "Whether to output a random colormap for ParaView in a separate XML file.");
  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "Aggregates");          //< factory which created aggregates
  Input(fineLevel, "DofsPerNode");         //< CoalesceAndDropFactory (needed for DofsPerNode variable)
  Input(fineLevel, "UnAmalgamationInfo");  //< AmalgamationFactory (needed for UnAmalgamationInfo variable)

  const ParameterList& pL = GetParameterList();
  // Only pull in coordinates if the user explicitly requests direct VTK output, so as not to break uses of old code
  if (pL.isParameter("aggregation: output filename") && pL.get<std::string>("aggregation: output filename").length()) {
    Input(fineLevel, "Coordinates");
    Input(fineLevel, "A");
    Input(fineLevel, "Graph");
    if (pL.get<bool>("aggregation: output file: coarse graph edges")) {
      Input(coarseLevel, "Coordinates");
      Input(coarseLevel, "A");
      Input(coarseLevel, "Graph");
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  using namespace std;
  // Decide which build function to follow, based on input params
  const ParameterList& pL = GetParameterList();
  FactoryMonitor m(*this, "AggregationExportFactory", coarseLevel);
  Teuchos::RCP<Aggregates> aggregates          = Get<Teuchos::RCP<Aggregates> >(fineLevel, "Aggregates");
  Teuchos::RCP<const Teuchos::Comm<int> > comm = aggregates->GetMap()->getComm();
  int numProcs                                 = comm->getSize();
  int myRank                                   = comm->getRank();
  string masterFilename                        = pL.get<std::string>("aggregation: output filename");  // filename parameter from master list
  string pvtuFilename                          = "";                                                   // only root processor will set this
  string localFilename                         = pL.get<std::string>("Output filename");
  string filenameToWrite;
  bool useVTK         = false;
  doCoarseGraphEdges_ = pL.get<bool>("aggregation: output file: coarse graph edges");
  doFineGraphEdges_   = pL.get<bool>("aggregation: output file: fine graph edges");
  if (masterFilename.length()) {
    useVTK          = true;
    filenameToWrite = masterFilename;
    if (filenameToWrite.rfind(".vtu") == string::npos)  // Must have the file extension in the name
      filenameToWrite.append(".vtu");
    if (numProcs > 1 && filenameToWrite.rfind("%PROCID") == string::npos)  // filename can't be identical between processsors in parallel problem
      filenameToWrite.insert(filenameToWrite.rfind(".vtu"), "-proc%PROCID");
  } else
    filenameToWrite = localFilename;
  LocalOrdinal DofsPerNode                 = Get<LocalOrdinal>(fineLevel, "DofsPerNode");
  Teuchos::RCP<AmalgamationInfo> amalgInfo = Get<RCP<AmalgamationInfo> >(fineLevel, "UnAmalgamationInfo");
  Teuchos::RCP<Matrix> Amat                = Get<RCP<Matrix> >(fineLevel, "A");
  Teuchos::RCP<Matrix> Ac;
  if (doCoarseGraphEdges_)
    Ac = Get<RCP<Matrix> >(coarseLevel, "A");
  Teuchos::RCP<CoordinateMultiVector> coords       = Teuchos::null;
  Teuchos::RCP<CoordinateMultiVector> coordsCoarse = Teuchos::null;
  Teuchos::RCP<GraphBase> fineGraph                = Teuchos::null;
  Teuchos::RCP<GraphBase> coarseGraph              = Teuchos::null;
  if (doFineGraphEdges_)
    fineGraph = Get<RCP<GraphBase> >(fineLevel, "Graph");
  if (doCoarseGraphEdges_)
    coarseGraph = Get<RCP<GraphBase> >(coarseLevel, "Graph");
  if (useVTK)  // otherwise leave null, will not be accessed by non-vtk code
  {
    coords  = Get<RCP<CoordinateMultiVector> >(fineLevel, "Coordinates");
    coords_ = coords;
    if (doCoarseGraphEdges_)
      coordsCoarse = Get<RCP<CoordinateMultiVector> >(coarseLevel, "Coordinates");
    dims_ = coords->getNumVectors();  // 2D or 3D?
    if (numProcs > 1) {
      if (aggregates->AggregatesCrossProcessors()) {  // Do we want to use the map from aggregates here instead of the map from A? Using the map from A seems to be problematic with multiple dofs per node
        RCP<Import> coordImporter                = Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(coords->getMap(), Amat->getColMap());
        RCP<CoordinateMultiVector> ghostedCoords = Xpetra::MultiVectorFactory<coordinate_type, LocalOrdinal, GlobalOrdinal, Node>::Build(Amat->getColMap(), dims_);
        ghostedCoords->doImport(*coords, *coordImporter, Xpetra::INSERT);
        coords  = ghostedCoords;
        coords_ = ghostedCoords;
      }
      if (doCoarseGraphEdges_) {
        RCP<Import> coordImporter                = Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(coordsCoarse->getMap(), Ac->getColMap());
        RCP<CoordinateMultiVector> ghostedCoords = Xpetra::MultiVectorFactory<coordinate_type, LocalOrdinal, GlobalOrdinal, Node>::Build(Ac->getColMap(), dims_);
        ghostedCoords->doImport(*coordsCoarse, *coordImporter, Xpetra::INSERT);
        coordsCoarse  = ghostedCoords;
        coordsCoarse_ = ghostedCoords;
      }
    }
  }
  GetOStream(Runtime0) << "AggregationExportFactory: DofsPerNode: " << DofsPerNode << std::endl;
  Teuchos::RCP<LocalOrdinalMultiVector> vertex2AggId_vector = aggregates->GetVertex2AggId();
  Teuchos::RCP<LocalOrdinalVector> procWinner_vector        = aggregates->GetProcWinner();
  Teuchos::ArrayRCP<LocalOrdinal> vertex2AggId              = aggregates->GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LocalOrdinal> procWinner                = aggregates->GetProcWinner()->getDataNonConst(0);

  vertex2AggId_ = vertex2AggId;

  // prepare for calculating global aggregate ids
  std::vector<GlobalOrdinal> numAggsGlobal(numProcs, 0);
  std::vector<GlobalOrdinal> numAggsLocal(numProcs, 0);
  std::vector<GlobalOrdinal> minGlobalAggId(numProcs, 0);

  numAggsLocal[myRank] = aggregates->GetNumAggregates();
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, numProcs, &numAggsLocal[0], &numAggsGlobal[0]);
  for (int i = 1; i < Teuchos::as<int>(numAggsGlobal.size()); ++i) {
    numAggsGlobal[i] += numAggsGlobal[i - 1];
    minGlobalAggId[i] = numAggsGlobal[i - 1];
  }
  if (numProcs == 0)
    aggsOffset_ = 0;
  else
    aggsOffset_ = minGlobalAggId[myRank];
  ArrayRCP<LO> aggStart;
  ArrayRCP<GlobalOrdinal> aggToRowMap;
  amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMap);
  int timeStep    = pL.get<int>("Output file: time step");
  int iter        = pL.get<int>("Output file: iter");
  filenameToWrite = this->replaceAll(filenameToWrite, "%LEVELID", toString(fineLevel.GetLevelID()));
  filenameToWrite = this->replaceAll(filenameToWrite, "%TIMESTEP", toString(timeStep));
  filenameToWrite = this->replaceAll(filenameToWrite, "%ITER", toString(iter));
  // Proc id MUST be included in vtu filenames to distinguish them (if multiple procs)
  // In all other cases (else), including processor # in filename is optional
  string masterStem = "";
  if (useVTK) {
    masterStem = filenameToWrite.substr(0, filenameToWrite.rfind(".vtu"));
    masterStem = this->replaceAll(masterStem, "%PROCID", "");
  }
  pvtuFilename     = masterStem + "-master.pvtu";
  string baseFname = filenameToWrite;  // get a version of the filename string with the %PROCID token, but without substituting myRank (needed for pvtu output)
  filenameToWrite  = this->replaceAll(filenameToWrite, "%PROCID", toString(myRank));
  GetOStream(Runtime0) << "AggregationExportFactory: outputfile \"" << filenameToWrite << "\"" << std::endl;
  ofstream fout(filenameToWrite.c_str());
  GO numAggs = aggregates->GetNumAggregates();
  if (!useVTK) {
    GO indexBase = aggregates->GetMap()->getIndexBase();  // extract indexBase from overlapping map within aggregates structure. The indexBase is constant throughout the whole simulation (either 0 = C++ or 1 = Fortran)
    GO offset    = amalgInfo->GlobalOffset();             // extract offset for global dof ids
    vector<GlobalOrdinal> nodeIds;
    for (int i = 0; i < numAggs; ++i) {
      fout << "Agg " << minGlobalAggId[myRank] + i << " Proc " << myRank << ":";

      // TODO: Use k+=DofsPerNode instead of ++k and get rid of std::unique call afterwards
      for (int k = aggStart[i]; k < aggStart[i + 1]; ++k) {
        nodeIds.push_back((aggToRowMap[k] - offset - indexBase) / DofsPerNode + indexBase);
      }

      // remove duplicate entries from nodeids
      std::sort(nodeIds.begin(), nodeIds.end());
      typename std::vector<GlobalOrdinal>::iterator endLocation = std::unique(nodeIds.begin(), nodeIds.end());
      nodeIds.erase(endLocation, nodeIds.end());

      // print out nodeids
      for (typename std::vector<GlobalOrdinal>::iterator printIt = nodeIds.begin(); printIt != nodeIds.end(); printIt++)
        fout << " " << *printIt;
      nodeIds.clear();
      fout << std::endl;
    }
  } else {
    using namespace std;
    // Make sure we have coordinates
    TEUCHOS_TEST_FOR_EXCEPTION(coords.is_null(), Exceptions::RuntimeError, "AggExportFactory could not get coordinates, but they are required for VTK output.");
    numAggs_  = numAggs;
    numNodes_ = coords->getLocalLength();
    // Get the sizes of the aggregates to speed up grabbing node IDs
    aggSizes_       = aggregates->ComputeAggregateSizesArrayRCP();
    myRank_         = myRank;
    string aggStyle = "Point Cloud";
    try {
      aggStyle = pL.get<string>("aggregation: output file: agg style");  // Let "Point Cloud" be the default style
    } catch (std::exception& e) {
    }
    vector<int> vertices;
    vector<int> geomSizes;
    string indent = "";
    nodeMap_      = Amat->getMap();
    for (LocalOrdinal i = 0; i < numNodes_; i++) {
      isRoot_.push_back(aggregates->IsRoot(i));
    }
    // If problem is serial and not outputting fine nor coarse graph edges, don't make pvtu file
    // Otherwise create it
    if (myRank == 0 && (numProcs != 1 || doCoarseGraphEdges_ || doFineGraphEdges_)) {
      ofstream pvtu(pvtuFilename.c_str());
      writePVTU_(pvtu, baseFname, numProcs);
      pvtu.close();
    }
    if (aggStyle == "Point Cloud")
      this->doPointCloud(vertices, geomSizes, numAggs_, numNodes_);
    else if (aggStyle == "Jacks")
      this->doJacks(vertices, geomSizes, numAggs_, numNodes_, isRoot_, vertex2AggId_);
    else if (aggStyle == "Jacks++")  // Not actually implemented
      doJacksPlus_(vertices, geomSizes);
    else if (aggStyle == "Convex Hulls")
      doConvexHulls(vertices, geomSizes);
    else {
      GetOStream(Warnings0) << "   Unrecognized agg style.\n   Possible values are Point Cloud, Jacks, Jacks++, and Convex Hulls.\n   Defaulting to Point Cloud." << std::endl;
      aggStyle = "Point Cloud";
      this->doPointCloud(vertices, geomSizes, numAggs_, numNodes_);
    }
    writeFile_(fout, aggStyle, vertices, geomSizes);
    if (doCoarseGraphEdges_) {
      string fname     = filenameToWrite;
      string cEdgeFile = fname.insert(fname.rfind(".vtu"), "-coarsegraph");
      std::ofstream edgeStream(cEdgeFile.c_str());
      doGraphEdges_(edgeStream, Ac, coarseGraph, false, DofsPerNode);
    }
    if (doFineGraphEdges_) {
      string fname     = filenameToWrite;
      string fEdgeFile = fname.insert(fname.rfind(".vtu"), "-finegraph");
      std::ofstream edgeStream(fEdgeFile.c_str());
      doGraphEdges_(edgeStream, Amat, fineGraph, true, DofsPerNode);
    }
    if (myRank == 0 && pL.get<bool>("aggregation: output file: build colormap")) {
      buildColormap_();
    }
  }
  fout.close();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doJacksPlus_(std::vector<int>& /* vertices */, std::vector<int>& /* geomSizes */) const {
  // TODO
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doConvexHulls(std::vector<int>& vertices, std::vector<int>& geomSizes) const {
  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> xCoords = coords_->getData(0);
  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> yCoords = coords_->getData(1);
  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> zCoords = Teuchos::null;

  if (dims_ == 2) {
    this->doConvexHulls2D(vertices, geomSizes, numAggs_, numNodes_, isRoot_, vertex2AggId_, xCoords, yCoords);
  } else {
    zCoords = coords_->getData(2);
    this->doConvexHulls3D(vertices, geomSizes, numAggs_, numNodes_, isRoot_, vertex2AggId_, xCoords, yCoords, zCoords);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doGraphEdges_(std::ofstream& fout, Teuchos::RCP<Matrix>& A, Teuchos::RCP<GraphBase>& G, bool fine, int dofs) const {
  using namespace std;
  ArrayView<const Scalar> values;
  ArrayView<const LocalOrdinal> neighbors;
  // Allow two different colors of connections (by setting "aggregates" scalar to CONTRAST_1 or CONTRAST_2)
  vector<pair<int, int> > vert1;  // vertices (node indices)
  vector<pair<int, int> > vert2;  // size of every cell is assumed to be 2 vertices, since all edges are drawn as lines

  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> xCoords = coords_->getData(0);
  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> yCoords = coords_->getData(1);
  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> zCoords = Teuchos::null;
  if (dims_ == 3)
    zCoords = coords_->getData(2);

  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> cx = coordsCoarse_->getData(0);
  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> cy = coordsCoarse_->getData(1);
  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> cz = Teuchos::null;
  if (dims_ == 3)
    cz = coordsCoarse_->getData(2);

  if (A->isGloballyIndexed()) {
    ArrayView<const GlobalOrdinal> indices;
    for (GlobalOrdinal globRow = 0; globRow < GlobalOrdinal(A->getGlobalNumRows()); globRow++) {
      if (dofs == 1)
        A->getGlobalRowView(globRow, indices, values);
      neighbors = G->getNeighborVertices((LocalOrdinal)globRow);
      int gEdge = 0;
      int aEdge = 0;
      while (gEdge != int(neighbors.size())) {
        if (dofs == 1) {
          if (neighbors[gEdge] == indices[aEdge]) {
            // graph and matrix both have this edge, wasn't filtered, show as color 1
            vert1.push_back(pair<int, int>(int(globRow), neighbors[gEdge]));
            gEdge++;
            aEdge++;
          } else {
            // graph contains an edge at gEdge which was filtered from A, show as color 2
            vert2.push_back(pair<int, int>(int(globRow), neighbors[gEdge]));
            gEdge++;
          }
        } else  // for multiple DOF problems, don't try to detect filtered edges and ignore A
        {
          vert1.push_back(pair<int, int>(int(globRow), neighbors[gEdge]));
          gEdge++;
        }
      }
    }
  } else {
    ArrayView<const LocalOrdinal> indices;
    for (LocalOrdinal locRow = 0; locRow < LocalOrdinal(A->getLocalNumRows()); locRow++) {
      if (dofs == 1)
        A->getLocalRowView(locRow, indices, values);
      neighbors = G->getNeighborVertices(locRow);
      // Add those local indices (columns) to the list of connections (which are pairs of the form (localM, localN))
      int gEdge = 0;
      int aEdge = 0;
      while (gEdge != int(neighbors.size())) {
        if (dofs == 1) {
          if (neighbors[gEdge] == indices[aEdge]) {
            vert1.push_back(pair<int, int>(locRow, neighbors[gEdge]));
            gEdge++;
            aEdge++;
          } else {
            vert2.push_back(pair<int, int>(locRow, neighbors[gEdge]));
            gEdge++;
          }
        } else {
          vert1.push_back(pair<int, int>(locRow, neighbors[gEdge]));
          gEdge++;
        }
      }
    }
  }
  for (size_t i = 0; i < vert1.size(); i++) {
    if (vert1[i].first > vert1[i].second) {
      int temp        = vert1[i].first;
      vert1[i].first  = vert1[i].second;
      vert1[i].second = temp;
    }
  }
  for (size_t i = 0; i < vert2.size(); i++) {
    if (vert2[i].first > vert2[i].second) {
      int temp        = vert2[i].first;
      vert2[i].first  = vert2[i].second;
      vert2[i].second = temp;
    }
  }
  sort(vert1.begin(), vert1.end());
  vector<pair<int, int> >::iterator newEnd = unique(vert1.begin(), vert1.end());  // remove duplicate edges
  vert1.erase(newEnd, vert1.end());
  sort(vert2.begin(), vert2.end());
  newEnd = unique(vert2.begin(), vert2.end());
  vert2.erase(newEnd, vert2.end());
  vector<int> points1;
  points1.reserve(2 * vert1.size());
  for (size_t i = 0; i < vert1.size(); i++) {
    points1.push_back(vert1[i].first);
    points1.push_back(vert1[i].second);
  }
  vector<int> points2;
  points2.reserve(2 * vert2.size());
  for (size_t i = 0; i < vert2.size(); i++) {
    points2.push_back(vert2[i].first);
    points2.push_back(vert2[i].second);
  }
  vector<int> unique1 = this->makeUnique(points1);
  vector<int> unique2 = this->makeUnique(points2);
  fout << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << endl;
  fout << "  <UnstructuredGrid>" << endl;
  fout << "    <Piece NumberOfPoints=\"" << unique1.size() + unique2.size() << "\" NumberOfCells=\"" << vert1.size() + vert2.size() << "\">" << endl;
  fout << "      <PointData Scalars=\"Node Aggregate Processor\">" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"Node\" format=\"ascii\">" << endl;  // node and aggregate will be set to CONTRAST_1|2, but processor will have its actual value
  string indent = "          ";
  fout << indent;
  for (size_t i = 0; i < unique1.size(); i++) {
    fout << CONTRAST_1_ << " ";
    if (i % 25 == 24)
      fout << endl
           << indent;
  }
  for (size_t i = 0; i < unique2.size(); i++) {
    fout << CONTRAST_2_ << " ";
    if ((i + 2 * vert1.size()) % 25 == 24)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"Aggregate\" format=\"ascii\">" << endl;
  fout << indent;
  for (size_t i = 0; i < unique1.size(); i++) {
    fout << CONTRAST_1_ << " ";
    if (i % 25 == 24)
      fout << endl
           << indent;
  }
  for (size_t i = 0; i < unique2.size(); i++) {
    fout << CONTRAST_2_ << " ";
    if ((i + 2 * vert2.size()) % 25 == 24)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << endl;
  fout << indent;
  for (size_t i = 0; i < unique1.size() + unique2.size(); i++) {
    fout << myRank_ << " ";
    if (i % 25 == 24)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "      </PointData>" << endl;
  fout << "      <Points>" << endl;
  fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
  fout << indent;
  for (size_t i = 0; i < unique1.size(); i++) {
    if (fine) {
      fout << xCoords[unique1[i]] << " " << yCoords[unique1[i]] << " ";
      if (dims_ == 3)
        fout << zCoords[unique1[i]] << " ";
      else
        fout << "0 ";
      if (i % 2)
        fout << endl
             << indent;
    } else {
      fout << cx[unique1[i]] << " " << cy[unique1[i]] << " ";
      if (dims_ == 3)
        fout << cz[unique1[i]] << " ";
      else
        fout << "0 ";
      if (i % 2)
        fout << endl
             << indent;
    }
  }
  for (size_t i = 0; i < unique2.size(); i++) {
    if (fine) {
      fout << xCoords[unique2[i]] << " " << yCoords[unique2[i]] << " ";
      if (dims_ == 3)
        fout << zCoords[unique2[i]] << " ";
      else
        fout << "0 ";
      if (i % 2)
        fout << endl
             << indent;
    } else {
      fout << cx[unique2[i]] << " " << cy[unique2[i]] << " ";
      if (dims_ == 3)
        fout << cz[unique2[i]] << " ";
      else
        fout << "0 ";
      if ((i + unique1.size()) % 2)
        fout << endl
             << indent;
    }
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "      </Points>" << endl;
  fout << "      <Cells>" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
  fout << indent;
  for (size_t i = 0; i < points1.size(); i++) {
    fout << points1[i] << " ";
    if (i % 10 == 9)
      fout << endl
           << indent;
  }
  for (size_t i = 0; i < points2.size(); i++) {
    fout << points2[i] + unique1.size() << " ";
    if ((i + points1.size()) % 10 == 9)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
  fout << indent;
  int offset = 0;
  for (size_t i = 0; i < vert1.size() + vert2.size(); i++) {
    offset += 2;
    fout << offset << " ";
    if (i % 10 == 9)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
  fout << indent;
  for (size_t i = 0; i < vert1.size() + vert2.size(); i++) {
    fout << "3 ";
    if (i % 25 == 24)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "      </Cells>" << endl;
  fout << "    </Piece>" << endl;
  fout << "  </UnstructuredGrid>" << endl;
  fout << "</VTKFile>" << endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeFile_(std::ofstream& fout, std::string styleName, std::vector<int>& vertices, std::vector<int>& geomSizes) const {
  using namespace std;

  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> xCoords = coords_->getData(0);
  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> yCoords = coords_->getData(1);
  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> zCoords = Teuchos::null;
  if (dims_ == 3)
    zCoords = coords_->getData(2);

  vector<int> uniqueFine = this->makeUnique(vertices);
  string indent          = "      ";
  fout << "<!--" << styleName << " Aggregates Visualization-->" << endl;
  fout << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << endl;
  fout << "  <UnstructuredGrid>" << endl;
  fout << "    <Piece NumberOfPoints=\"" << uniqueFine.size() << "\" NumberOfCells=\"" << geomSizes.size() << "\">" << endl;
  fout << "      <PointData Scalars=\"Node Aggregate Processor\">" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"Node\" format=\"ascii\">" << endl;
  indent = "          ";
  fout << indent;
  bool localIsGlobal = GlobalOrdinal(nodeMap_->getGlobalNumElements()) == GlobalOrdinal(nodeMap_->getLocalNumElements());
  for (size_t i = 0; i < uniqueFine.size(); i++) {
    if (localIsGlobal) {
      fout << uniqueFine[i] << " ";  // if all nodes are on this processor, do not map from local to global
    } else
      fout << nodeMap_->getGlobalElement(uniqueFine[i]) << " ";
    if (i % 10 == 9)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"Aggregate\" format=\"ascii\">" << endl;
  fout << indent;
  for (size_t i = 0; i < uniqueFine.size(); i++) {
    if (vertex2AggId_[uniqueFine[i]] == -1)
      fout << vertex2AggId_[uniqueFine[i]] << " ";
    else
      fout << aggsOffset_ + vertex2AggId_[uniqueFine[i]] << " ";
    if (i % 10 == 9)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << endl;
  fout << indent;
  for (size_t i = 0; i < uniqueFine.size(); i++) {
    fout << myRank_ << " ";
    if (i % 20 == 19)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "      </PointData>" << endl;
  fout << "      <Points>" << endl;
  fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
  fout << indent;
  for (size_t i = 0; i < uniqueFine.size(); i++) {
    fout << xCoords[uniqueFine[i]] << " " << yCoords[uniqueFine[i]] << " ";
    if (dims_ == 2)
      fout << "0 ";
    else
      fout << zCoords[uniqueFine[i]] << " ";
    if (i % 3 == 2)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "      </Points>" << endl;
  fout << "      <Cells>" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
  fout << indent;
  for (size_t i = 0; i < vertices.size(); i++) {
    fout << vertices[i] << " ";
    if (i % 10 == 9)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
  fout << indent;
  int accum = 0;
  for (size_t i = 0; i < geomSizes.size(); i++) {
    accum += geomSizes[i];
    fout << accum << " ";
    if (i % 10 == 9)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
  fout << indent;
  for (size_t i = 0; i < geomSizes.size(); i++) {
    switch (geomSizes[i]) {
      case 1:
        fout << "1 ";  // Point
        break;
      case 2:
        fout << "3 ";  // Line
        break;
      case 3:
        fout << "5 ";  // Triangle
        break;
      default:
        fout << "7 ";  // Polygon
    }
    if (i % 30 == 29)
      fout << endl
           << indent;
  }
  fout << endl;
  fout << "        </DataArray>" << endl;
  fout << "      </Cells>" << endl;
  fout << "    </Piece>" << endl;
  fout << "  </UnstructuredGrid>" << endl;
  fout << "</VTKFile>" << endl;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::buildColormap_() const {
  using namespace std;
  try {
    ofstream color("random-colormap.xml");
    color << "<ColorMap name=\"MueLu-Random\" space=\"RGB\">" << endl;
    // Give -1, -2, -3 distinctive colors (so that the style functions can have constrasted geometry types)
    // Do red, orange, yellow to constrast with cool color spectrum for other types
    color << "  <Point x=\"" << CONTRAST_1_ << "\" o=\"1\" r=\"1\" g=\"0\" b=\"0\"/>" << endl;
    color << "  <Point x=\"" << CONTRAST_2_ << "\" o=\"1\" r=\"1\" g=\"0.6\" b=\"0\"/>" << endl;
    color << "  <Point x=\"" << CONTRAST_3_ << "\" o=\"1\" r=\"1\" g=\"1\" b=\"0\"/>" << endl;
    srand(time(NULL));
    for (int i = 0; i < 5000; i += 4) {
      color << "  <Point x=\"" << i << "\" o=\"1\" r=\"" << (rand() % 50) / 256.0 << "\" g=\"" << (rand() % 256) / 256.0 << "\" b=\"" << (rand() % 256) / 256.0 << "\"/>" << endl;
    }
    color << "</ColorMap>" << endl;
    color.close();
  } catch (std::exception& e) {
    GetOStream(Warnings0) << "   Error while building colormap file: " << e.what() << endl;
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writePVTU_(std::ofstream& pvtu, std::string baseFname, int numProcs) const {
  using namespace std;
  // If using vtk, filenameToWrite now contains final, correct ***.vtu filename (for the current proc)
  // So the root proc will need to use its own filenameToWrite to make a list of the filenames of all other procs to put in
  // pvtu file.
  pvtu << "<VTKFile type=\"PUnstructuredGrid\" byte_order=\"LittleEndian\">" << endl;
  pvtu << "  <PUnstructuredGrid GhostLevel=\"0\">" << endl;
  pvtu << "    <PPointData Scalars=\"Node Aggregate Processor\">" << endl;
  pvtu << "      <PDataArray type=\"Int32\" Name=\"Node\"/>" << endl;
  pvtu << "      <PDataArray type=\"Int32\" Name=\"Aggregate\"/>" << endl;
  pvtu << "      <PDataArray type=\"Int32\" Name=\"Processor\"/>" << endl;
  pvtu << "    </PPointData>" << endl;
  pvtu << "    <PPoints>" << endl;
  pvtu << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>" << endl;
  pvtu << "    </PPoints>" << endl;
  for (int i = 0; i < numProcs; i++) {
    // specify the piece for each proc (the replaceAll expression matches what the filenames will be on other procs)
    pvtu << "    <Piece Source=\"" << this->replaceAll(baseFname, "%PROCID", toString(i)) << "\"/>" << endl;
  }
  // reference files with graph pieces, if applicable
  if (doFineGraphEdges_) {
    for (int i = 0; i < numProcs; i++) {
      string fn = this->replaceAll(baseFname, "%PROCID", toString(i));
      pvtu << "    <Piece Source=\"" << fn.insert(fn.rfind(".vtu"), "-finegraph") << "\"/>" << endl;
    }
  }
  if (doCoarseGraphEdges_) {
    for (int i = 0; i < numProcs; i++) {
      string fn = this->replaceAll(baseFname, "%PROCID", toString(i));
      pvtu << "    <Piece Source=\"" << fn.insert(fn.rfind(".vtu"), "-coarsegraph") << "\"/>" << endl;
    }
  }
  pvtu << "  </PUnstructuredGrid>" << endl;
  pvtu << "</VTKFile>" << endl;
  pvtu.close();
}
}  // namespace MueLu

#endif /* MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_ */
