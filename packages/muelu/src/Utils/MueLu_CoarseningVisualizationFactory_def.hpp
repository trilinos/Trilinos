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

#ifndef MUELU_COARSENINGVISUALIZATIONFACTORY_DEF_HPP_
#define MUELU_COARSENINGVISUALIZATIONFACTORY_DEF_HPP_

#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include "MueLu_CoarseningVisualizationFactory_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_Graph.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"
#include <vector>
#include <list>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <cstdio>
#include <cmath>

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    std::string output_msg = "Output filename template (%TIMESTEP is replaced by \'Output file: time step\' variable,"
        "%ITER is replaced by \'Output file: iter\' variable, %LEVELID is replaced level id, %PROCID is replaced by processor id)";
    std::string output_def = "aggs_level%LEVELID_proc%PROCID.out";

    //validParamList->set< RCP<const FactoryBase> >("A", Teuchos::null, "Factory for A.");
    validParamList->set< RCP<const FactoryBase> >("P",              Teuchos::null, "Tentative prolongator factory");
    validParamList->set< RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory for Coordinates.");
    //validParamList->set< RCP<const FactoryBase> >("Graph", Teuchos::null, "Factory for Graph.");
    //validParamList->set< RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Factory for Aggregates.");
    //validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Factory for UnAmalgamationInfo.");
    //validParamList->set< RCP<const FactoryBase> >("DofsPerNode", Teuchos::null, "Factory for DofsPerNode.");

    // New-style master list options (here are same defaults as in masterList.xml)
    validParamList->set< std::string >           ("visualization: output filename",                    "viz%LEVELID",                    "filename for VTK-style visualization output");
    //validParamList->set< int >                   ("aggregation: output file: time step",             0,                     "time step variable for output file name");// Remove me?
    //validParamList->set< int >                   ("aggregation: output file: iter",                  0,                     "nonlinear iteration variable for output file name");//Remove me?
    //validParamList->set<std::string>             ("aggregation: output file: agg style",             "Point Cloud",         "style of aggregate visualization for VTK output");
    //validParamList->set<bool>                    ("aggregation: output file: fine graph edges",      false,                 "Whether to draw all fine node connections along with the aggregates.");
    //validParamList->set<bool>                    ("aggregation: output file: coarse graph edges",    false,                 "Whether to draw all coarse node connections along with the aggregates.");
    //validParamList->set<bool>                    ("aggregation: output file: build colormap",        false,                 "Whether to output a random colormap for ParaView in a separate XML file.");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    //Input(fineLevel, "Aggregates");         //< factory which created aggregates
    //Input(fineLevel, "DofsPerNode");        //< CoalesceAndDropFactory (needed for DofsPerNode variable)
    //Input(fineLevel, "UnAmalgamationInfo"); //< AmalgamationFactory (needed for UnAmalgamationInfo variable)

    //const ParameterList & pL = GetParameterList();
    Input(fineLevel, "Coordinates");
    //Input(fineLevel, "A");
    //Input(fineLevel, "Graph");
    Input(coarseLevel, "P");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {

    RCP<Matrix>           P             = Get< RCP<Matrix> >          (coarseLevel, "P");

    RCP<const Teuchos::Comm<int> > comm = P->getRowMap()->getComm();

    LocalOrdinal dofsPerNode = 1;
    LocalOrdinal stridedRowOffset = 0;
    RCP<const StridedMap> strRowMap    = Teuchos::null;
    if (P->IsView("stridedMaps") && Teuchos::rcp_dynamic_cast<const StridedMap>(P->getRowMap("stridedMaps")) != Teuchos::null) {
      strRowMap = Teuchos::rcp_dynamic_cast<const StridedMap>(P->getRowMap("stridedMaps"));
      LocalOrdinal blockid       = strRowMap->getStridedBlockId();
      if (blockid > -1) {
        std::vector<size_t> stridingInfo = strRowMap->getStridingData();
        for (size_t j = 0; j < Teuchos::as<size_t>(blockid); j++)
          stridedRowOffset += stridingInfo[j];
        dofsPerNode = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);
      } else {
        dofsPerNode = strRowMap->getFixedBlockSize();
      }
      GetOStream(Runtime1) << "CoarseningVisualizationFactory::Build():" << " #dofs per node = " << dofsPerNode << std::endl;
    }

    LocalOrdinal columnsPerNode = dofsPerNode;
    LocalOrdinal stridedColumnOffset = 0;
    RCP<const StridedMap> strDomainMap = Teuchos::null;
    if (P->IsView("stridedMaps") && Teuchos::rcp_dynamic_cast<const StridedMap>(P->getRowMap("stridedMaps")) != Teuchos::null) {
      strDomainMap = Teuchos::rcp_dynamic_cast<const StridedMap>(P->getColMap("stridedMaps"));
      LocalOrdinal blockid = strDomainMap->getStridedBlockId();

      if (blockid > -1) {
        std::vector<size_t> stridingInfo = strDomainMap->getStridingData();
        for (size_t j = 0; j < Teuchos::as<size_t>(blockid); j++)
          stridedColumnOffset += stridingInfo[j];
        columnsPerNode = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);
      } else {
        columnsPerNode = strDomainMap->getFixedBlockSize();
      }
      GetOStream(Runtime1) << "CoarseningVisualizationFactory::Build():" << " #columns per node = " << columnsPerNode << std::endl;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(strDomainMap->getNodeNumElements() != P->getColMap()->getNodeNumElements(), Exceptions::RuntimeError,
                                               "CoarseningVisualization only supports non-overlapping transfers");

    // number of local "aggregates"
    LocalOrdinal numLocalAggs = strDomainMap->getNodeNumElements() / columnsPerNode;
    std::vector< std::set<LocalOrdinal> > localAggs(numLocalAggs);

    // do loop over all local rows of prolongator and extract column information
    for (LO row = 0; row < Teuchos::as<LO>(P->getRowMap()->getNodeNumElements()); ++row) {
      //size_t nnz = P->getNumEntriesInLocalRow(row);
      ArrayView<const LO> indices;
      ArrayView<const SC> vals;
      P->getLocalRowView(row, indices, vals);

      for(typename ArrayView<const LO>::iterator c = indices.begin(); c != indices.end(); ++c) {
        localAggs[(*c)/columnsPerNode].insert(row/dofsPerNode);
      }
    }

    // determine number of "aggs" per proc and calculate local "agg" offset...
    std::vector<int> myLocalAggsPerProc(comm->getSize(),0);
    std::vector<int> numLocalAggsPerProc(comm->getSize(),0);
    myLocalAggsPerProc[comm->getRank()] = numLocalAggs;
    Teuchos::reduceAll(*comm,Teuchos::REDUCE_MAX,comm->getSize(),&myLocalAggsPerProc[0],&numLocalAggsPerProc[0]);

    LocalOrdinal myAggOffset = 0;
    for(int i = 0; i < comm->getRank(); ++i) {
      myAggOffset += numLocalAggsPerProc[i];
    }

    /*for (LocalOrdinal i = 0; i < numLocalAggs; ++i) {

      std::cout << "PROC: " << comm->getRank() << " Local aggregate: " << i + myAggOffset << " with nodes: ";
      for( typename std::set<LocalOrdinal>::iterator it = localAggs[i].begin(); it != localAggs[i].end(); ++it) {
        std::cout << *it << ", ";
      }
      std::cout << std::endl;
    }*/

    // get fine level coordinates
    Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > coords = Get<RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > >(fineLevel, "Coordinates");

    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<LO>(P->getRowMap()->getNodeNumElements()) / dofsPerNode != Teuchos::as<LocalOrdinal>(coords->getLocalLength()), Exceptions::RuntimeError,
                                           "Number of fine level nodes in coordinates is inconsistent with dof based information");

    Teuchos::RCP<const Map> nodeMap = coords->getMap();

    Teuchos::ArrayRCP<const double> xCoords = Teuchos::arcp_reinterpret_cast<const double>(coords->getData(0));
    Teuchos::ArrayRCP<const double> yCoords = Teuchos::arcp_reinterpret_cast<const double>(coords->getData(1));
    Teuchos::ArrayRCP<const double> zCoords = Teuchos::null;
    if(coords->getNumVectors() == 3) {
      zCoords = Teuchos::arcp_reinterpret_cast<const double>(coords->getData(2));
    }


    LocalOrdinal numFineNodes = Teuchos::as<LocalOrdinal>(coords->getLocalLength());

    std::vector<LocalOrdinal> vertex2AggId(numFineNodes, -1);
    for (LocalOrdinal i = 0; i < numLocalAggs; ++i) {
      // TODO: check if entry = -1
      for( typename std::set<LocalOrdinal>::iterator it = localAggs[i].begin(); it != localAggs[i].end(); ++it) {
        vertex2AggId[*it] = i;
      }
    }

    std::vector<bool> isRoot(numFineNodes, false);

    // we have no information which node is root and which is not
    // we could either look at the entries in P again or build some new root nodes
    // assuming that root nodes are usually at the center of the aggregate
    for (LocalOrdinal i = 0; i < numLocalAggs; ++i) {

      double xCenter = 0.0;
      double yCenter = 0.0;
      double zCenter = 0.0;

      // loop over all nodes in aggregate i and determine center coordinates of aggregate
      for( typename std::set<LocalOrdinal>::iterator it = localAggs[i].begin(); it != localAggs[i].end(); ++it) {
        xCenter += xCoords[*it];
        yCenter += yCoords[*it];
        if(coords->getNumVectors() == 3) zCenter += zCoords[*it];
      }
      xCenter /= Teuchos::as<LocalOrdinal>(localAggs[i].size());
      yCenter /= Teuchos::as<LocalOrdinal>(localAggs[i].size());
      zCenter /= Teuchos::as<LocalOrdinal>(localAggs[i].size());

      // loop over all nodes in aggregate i and find node which is closest to aggregate center
      LocalOrdinal rootCandidate = -1;
      double minDistance = Teuchos::ScalarTraits<double>::one() / Teuchos::ScalarTraits<double>::sfmin();
      for( typename std::set<LocalOrdinal>::iterator it = localAggs[i].begin(); it != localAggs[i].end(); ++it) {
        double tempx = xCenter - xCoords[*it];
        double tempy = yCenter - yCoords[*it];
        double tempz = 0.0;
        if(coords->getNumVectors() == 3) tempz = zCenter - zCoords[*it];
        double distance = 0.0;
        distance += tempx*tempx;
        distance += tempy*tempy;
        distance += tempz*tempz;
        distance = sqrt(distance);
        if (distance <= minDistance) {
          minDistance = distance;
          rootCandidate = *it;
        }
      }

      isRoot[rootCandidate] = true;
    }


    std::vector<LocalOrdinal> vertices;
    std::vector<LocalOrdinal> geomSize;

    // TODO replace this...
#if 0
    vertices.reserve(numFineNodes);
    geomSize.reserve(numFineNodes);
    for(LocalOrdinal i = 0; i < numFineNodes; i++)
    {
      vertices.push_back(i);
      geomSize.push_back(1);
    }
#endif

#if 1
    //For each aggregate, find the root node then connect it with the other nodes in the aggregate
    //Doesn't matter the order, as long as all edges are found.
    vertices.reserve(vertices.size() + 3 * (numFineNodes - numLocalAggs));
    geomSize.reserve(vertices.size() + 2 * (numFineNodes - numLocalAggs));
    int root = 0;
    for(int i = 0; i < numLocalAggs; i++) //TODO: Replace this O(n^2) with a better way
    {
      while(!isRoot[root])
        root++;
      int numInAggFound = 0;
      for(int j = 0; j < numFineNodes; j++)
      {
        if(j == root) //don't make a connection from the root to itself
        {
          numInAggFound++;
          continue;
        }
        if(vertex2AggId[root] == vertex2AggId[j])
        {
          vertices.push_back(root);
          vertices.push_back(j);
          geomSize.push_back(2);
          //Also draw the free endpoint explicitly for the current line
          vertices.push_back(j);
          geomSize.push_back(1);
          numInAggFound++;
          //if(numInAggFound == aggSizes_[vertex2AggId_[root]]) //don't spend more time looking if done with that root
          //  break;
        }
      }
      root++; //get set up to look for the next root
    }
#endif
    // end replace this

    const ParameterList & pL = GetParameterList();
    std::string masterFilename = pL.get<std::string>("visualization: output filename"); //filename parameter from master list
    std::string filenameToWrite;
    filenameToWrite = masterFilename;
    if(filenameToWrite.rfind(".vtu") == std::string::npos) //Must have the file extension in the name
      filenameToWrite.append(".vtu");
    if(comm->getSize() > 1 && filenameToWrite.rfind("%PROCID") == std::string::npos) //filename can't be identical between processsors in parallel problem
      filenameToWrite.insert(filenameToWrite.rfind(".vtu"), "-proc%PROCID");
    filenameToWrite = replaceAll(filenameToWrite, "%PROCID", toString(comm->getRank()));
    filenameToWrite = replaceAll(filenameToWrite, "%LEVELID",  toString(fineLevel.GetLevelID()));

    std::string styleName = "PointCloud";

    std::ofstream fout(filenameToWrite.c_str());
    std::vector<int> uniqueFine = makeUnique(vertices);
    std::string indent = "      ";
    fout << "<!--" << styleName << " Aggregates Visualization-->" << std::endl;
    fout << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << std::endl;
    fout << "  <UnstructuredGrid>" << std::endl;
    fout << "    <Piece NumberOfPoints=\"" << uniqueFine.size() << "\" NumberOfCells=\"" << geomSize.size() << "\">" << std::endl;
    fout << "      <PointData Scalars=\"Node Aggregate Processor\">" << std::endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Node\" format=\"ascii\">" << std::endl;
    indent = "          ";
    fout << indent;
    for(int i = 0; i < int(uniqueFine.size()); i++)
    {
      fout << nodeMap->getGlobalElement(uniqueFine[i]) << " ";
      if(i % 10 == 9)
        fout << std::endl << indent;
    }
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Aggregate\" format=\"ascii\">" << std::endl;
    fout << indent;
    for(int i = 0; i < int(uniqueFine.size()); i++)
    {
      fout <<  myAggOffset + vertex2AggId[uniqueFine[i]] << " ";
      if(i % 10 == 9)
        fout << std::endl << indent;
    }
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << std::endl;
    fout << indent;
    for(int i = 0; i < int(uniqueFine.size()); i++)
    {
      fout << comm->getRank() << " ";
      if(i % 20 == 19)
        fout << std::endl << indent;
    }
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    fout << "      </PointData>" << std::endl;
    fout << "      <Points>" << std::endl;
    fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    fout << indent;
    for(int i = 0; i < int(uniqueFine.size()); i++)
    {
      fout << xCoords[uniqueFine[i]] << " " << yCoords[uniqueFine[i]] << " ";
      if(coords->getNumVectors() == 2)
        fout << "0 ";
      else
        fout << zCoords[uniqueFine[i]] << " ";
      if(i % 3 == 2)
        fout << std::endl << indent;
    }
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    fout << "      </Points>" << std::endl;
    fout << "      <Cells>" << std::endl;
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    fout << indent;
    for(int i = 0; i < int(vertices.size()); i++)
    {
      fout << vertices[i] << " ";
      if(i % 10 == 9)
        fout << std::endl << indent;
    }
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    fout << indent;
    int accum = 0;
    for(int i = 0; i < int(geomSize.size()); i++)
    {
      accum += geomSize[i];
      fout << accum << " ";
      if(i % 10 == 9)
        fout << std::endl << indent;
    }
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    fout << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << std::endl;
    fout << indent;
    for(int i = 0; i < int(geomSize.size()); i++)
    {
      switch(geomSize[i])
      {
        case 1:
          fout << "1 "; //Point
          break;
        case 2:
          fout << "3 "; //Line
          break;
        case 3:
          fout << "5 "; //Triangle
          break;
        default:
          fout << "7 "; //Polygon
      }
      if(i % 30 == 29)
        fout << std::endl << indent;
    }
    fout << std::endl;
    fout << "        </DataArray>" << std::endl;
    fout << "      </Cells>" << std::endl;
    fout << "    </Piece>" << std::endl;
    fout << "  </UnstructuredGrid>" << std::endl;
    fout << "</VTKFile>" << std::endl;

    /*Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > coords = Teuchos::null;
    Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > coordsCoarse = Teuchos::null;
    Teuchos::RCP<GraphBase> fineGraph = Teuchos::null;
    Teuchos::RCP<GraphBase> coarseGraph = Teuchos::null;
    if(doFineGraphEdges_)
      fineGraph = Get<RCP<GraphBase> >(fineLevel, "Graph");
    if(doCoarseGraphEdges_)
      coarseGraph = Get<RCP<GraphBase> >(coarseLevel, "Graph");
    coords = Get<RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > >(fineLevel, "Coordinates");
    if(doCoarseGraphEdges_)
      coordsCoarse = Get<RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > >(coarseLevel, "Coordinates");
    dims_ = coords->getNumVectors();  //2D or 3D?*/


#if 0
    //Decide which build function to follow, based on input params
    const ParameterList& pL = GetParameterList();
    FactoryMonitor m(*this, "AggregationExportFactory", coarseLevel);
    Teuchos::RCP<Aggregates> aggregates      = Get< Teuchos::RCP<Aggregates> >(fineLevel,"Aggregates");
    Teuchos::RCP<const Teuchos::Comm<int> > comm = aggregates->GetMap()->getComm();
    int numProcs = comm->getSize();
    int myRank   = comm->getRank();
    string masterFilename = pL.get<std::string>("aggregation: output filename"); //filename parameter from master list
    string pvtuFilename = ""; //only root processor will set this
    string localFilename = pL.get<std::string>("Output filename");
    string filenameToWrite;
    bool useVTK = false;
    doCoarseGraphEdges_ = pL.get<bool>("aggregation: output file: coarse graph edges");
    doFineGraphEdges_ = pL.get<bool>("aggregation: output file: fine graph edges");
    if(masterFilename.length())
    {
      useVTK = true;
      filenameToWrite = masterFilename;
      if(filenameToWrite.rfind(".vtu") == string::npos) //Must have the file extension in the name
        filenameToWrite.append(".vtu");
      if(numProcs > 1 && filenameToWrite.rfind("%PROCID") == string::npos) //filename can't be identical between processsors in parallel problem
        filenameToWrite.insert(filenameToWrite.rfind(".vtu"), "-proc%PROCID");
    }
    else
      filenameToWrite = localFilename;
    LocalOrdinal DofsPerNode                 = Get< LocalOrdinal >            (fineLevel,"DofsPerNode");
    Teuchos::RCP<AmalgamationInfo> amalgInfo = Get< RCP<AmalgamationInfo> >   (fineLevel,"UnAmalgamationInfo");
    Teuchos::RCP<Matrix> Amat = Get<RCP<Matrix> >(fineLevel, "A");
    Teuchos::RCP<Matrix> Ac;
    if(doCoarseGraphEdges_)
      Ac = Get<RCP<Matrix> >(coarseLevel, "A");
    Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > coords = Teuchos::null;
    Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > coordsCoarse = Teuchos::null;
    Teuchos::RCP<GraphBase> fineGraph = Teuchos::null;
    Teuchos::RCP<GraphBase> coarseGraph = Teuchos::null;
    if(doFineGraphEdges_)
      fineGraph = Get<RCP<GraphBase> >(fineLevel, "Graph");
    if(doCoarseGraphEdges_)
      coarseGraph = Get<RCP<GraphBase> >(coarseLevel, "Graph");
    if(useVTK) //otherwise leave null, will not be accessed by non-vtk code
    {
      coords = Get<RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > >(fineLevel, "Coordinates");
      if(doCoarseGraphEdges_)
        coordsCoarse = Get<RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > >(coarseLevel, "Coordinates");
      dims_ = coords->getNumVectors();  //2D or 3D?
      if(numProcs > 1)
      {
        if(doFineGraphEdges_)
        {
          RCP<Import> coordImporter = Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(coords->getMap(), Amat->getColMap());
          RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > ghostedCoords = Xpetra::MultiVectorFactory<double, LocalOrdinal, GlobalOrdinal, Node>::Build(Amat->getColMap(), dims_);
          ghostedCoords->doImport(*coords, *coordImporter, Xpetra::INSERT);
          coords = ghostedCoords;
        }
        if(doCoarseGraphEdges_)
        {
          RCP<Import> coordImporter = Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(coordsCoarse->getMap(), Ac->getColMap());
          RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > ghostedCoords = Xpetra::MultiVectorFactory<double, LocalOrdinal, GlobalOrdinal, Node>::Build(Ac->getColMap(), dims_);
          ghostedCoords->doImport(*coordsCoarse, *coordImporter, Xpetra::INSERT);
          coordsCoarse = ghostedCoords;
        }
      }
    }
    GetOStream(Runtime0) << "AggregationExportFactory: DofsPerNode: " << DofsPerNode << std::endl;
    Teuchos::RCP<LocalOrdinalVector> vertex2AggId_vector = aggregates->GetVertex2AggId();
    Teuchos::RCP<LocalOrdinalVector> procWinner_vector   = aggregates->GetProcWinner();
    Teuchos::ArrayRCP<LocalOrdinal>  vertex2AggId        = aggregates->GetVertex2AggId()->getDataNonConst(0);
    Teuchos::ArrayRCP<LocalOrdinal>  procWinner          = aggregates->GetProcWinner()->getDataNonConst(0);

    vertex2AggId_ = vertex2AggId;

    // prepare for calculating global aggregate ids
    std::vector<GlobalOrdinal> numAggsGlobal (numProcs, 0);
    std::vector<GlobalOrdinal> numAggsLocal  (numProcs, 0);
    std::vector<GlobalOrdinal> minGlobalAggId(numProcs, 0);

    numAggsLocal[myRank] = aggregates->GetNumAggregates();
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, numProcs, &numAggsLocal[0], &numAggsGlobal[0]);
    for (int i = 1; i < Teuchos::as<int>(numAggsGlobal.size()); ++i)
    {
      numAggsGlobal [i] += numAggsGlobal[i-1];
      minGlobalAggId[i]  = numAggsGlobal[i-1];
    }
    if(myRank == 0)
      aggsOffset_ = 0;
    else
      aggsOffset_ = minGlobalAggId[myRank];
    ArrayRCP<LO>            aggStart;
    ArrayRCP<GlobalOrdinal> aggToRowMap;
    amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMap);
    int timeStep = pL.get< int > ("Output file: time step");
    int iter = pL.get< int > ("Output file: iter");
    filenameToWrite = replaceAll(filenameToWrite, "%LEVELID",  toString(fineLevel.GetLevelID()));
    filenameToWrite = replaceAll(filenameToWrite, "%TIMESTEP", toString(timeStep));
    filenameToWrite = replaceAll(filenameToWrite, "%ITER",     toString(iter));
    //Proc id MUST be included in vtu filenames to distinguish them (if multiple procs)
    //In all other cases (else), including processor # in filename is optional
    string masterStem = "";
    if(useVTK)
    {
      masterStem = filenameToWrite.substr(0, filenameToWrite.rfind(".vtu"));
      masterStem = replaceAll(masterStem, "%PROCID", "");
    }
    pvtuFilename = masterStem + "-master.pvtu";
    string baseFname = filenameToWrite;  //get a version of the filename string with the %PROCID token, but without substituting myRank (needed for pvtu output)
    filenameToWrite = replaceAll(filenameToWrite, "%PROCID", toString(myRank));
    GetOStream(Runtime0) << "AggregationExportFactory: outputfile \"" << filenameToWrite << "\"" << std::endl;
    ofstream fout(filenameToWrite.c_str());
    GO numAggs = aggregates->GetNumAggregates();
    if(!useVTK)
    {
      GO indexBase = aggregates->GetMap()->getIndexBase(); // extract indexBase from overlapping map within aggregates structure. The indexBase is constant throughout the whole simulation (either 0 = C++ or 1 = Fortran)
      GO offset    = amalgInfo->GlobalOffset();            // extract offset for global dof ids
      vector<GlobalOrdinal> nodeIds;
      for (int i = 0; i < numAggs; ++i) {
        fout << "Agg " << minGlobalAggId[myRank] + i << " Proc " << myRank << ":";

        // TODO: Use k+=DofsPerNode instead of ++k and get rid of std::unique call afterwards
        for (int k = aggStart[i]; k < aggStart[i+1]; ++k) {
          nodeIds.push_back((aggToRowMap[k] - offset - indexBase) / DofsPerNode + indexBase);
        }

        // remove duplicate entries from nodeids
        std::sort(nodeIds.begin(), nodeIds.end());
        typename std::vector<GlobalOrdinal>::iterator endLocation = std::unique(nodeIds.begin(), nodeIds.end());
        nodeIds.erase(endLocation, nodeIds.end());

        // print out nodeids 
        for(typename std::vector<GlobalOrdinal>::iterator printIt = nodeIds.begin(); printIt != nodeIds.end(); printIt++)
          fout << " " << *printIt;
        nodeIds.clear();
        fout << std::endl;
      }
    }
    else
    {
      using namespace std;
      if(sizeof(Scalar) != sizeof(double)) //Don't try to process complex, from now on will assume scalar is double.
        throw runtime_error("Only double scalars not supported in VTK aggregate visualization.");
      //Make sure we have coordinates
      if(coords.is_null())
        throw runtime_error("AggExportFactory could not get coordinates, but they are required for VTK output.");
      numAggs_ = numAggs;
      numNodes_ = coords->getLocalLength();
      //get access to the coord data
      xCoords_ = Teuchos::arcp_reinterpret_cast<const double>(coords->getData(0));
      yCoords_ = Teuchos::arcp_reinterpret_cast<const double>(coords->getData(1));
      zCoords_ = Teuchos::null;
      if(doCoarseGraphEdges_)
      {
        cx_ = Teuchos::arcp_reinterpret_cast<const double>(coordsCoarse->getData(0));
        cy_ = Teuchos::arcp_reinterpret_cast<const double>(coordsCoarse->getData(1));
        cz_ = Teuchos::null;
      }
      if(dims_ == 3)
      {
        zCoords_ = Teuchos::arcp_reinterpret_cast<const double>(coords->getData(2));
        if(doCoarseGraphEdges_)
          cz_ = Teuchos::arcp_reinterpret_cast<const double>(coordsCoarse->getData(2));
      }
      //Get the sizes of the aggregates to speed up grabbing node IDs
      aggSizes_ = aggregates->ComputeAggregateSizes();
      myRank_ = myRank;
      string aggStyle = "Point Cloud";
      try
      {
        aggStyle = pL.get<string>("aggregation: output file: agg style"); //Let "Point Cloud" be the default style
      }
      catch(exception& e) {}
      vector<int> vertices;
      vector<int> geomSizes;
      vector<int> verticesCoarse;
      vector<int> geomSizesCoarse;
      string indent = "";
      nodeMap_ = Amat->getMap();
      for(LocalOrdinal i = 0; i < numNodes_; i++)
      {
        isRoot_.push_back(aggregates->IsRoot(i));
      }
      //If the output will be parallel VTU and this is root processor, write out the master PVTU file now
      //Only do this if the problem has multiple processors, since the .vtu from the root proc would work on its own in a serial problem
      if(myRank == 0 && (numProcs != 1 || doCoarseGraphEdges_ || doFineGraphEdges_))
      {
        ofstream pvtu(pvtuFilename.c_str());
        writePVTU_(pvtu, baseFname, numProcs);
        pvtu.close();
      }
      if(aggStyle == "Point Cloud")
        doPointCloud_(vertices, geomSizes);
      else if(aggStyle == "Jacks")
        doJacks_(vertices, geomSizes);
      else if(aggStyle == "Jacks++") //Not actually implemented
        doJacksPlus_(vertices, geomSizes);
      else if(aggStyle == "Convex Hulls")
        doConvexHulls_(vertices, geomSizes);
      else if(aggStyle == "Alpha Hulls") //Not implemented
        doAlphaHulls_(vertices, geomSizes);
      else
      {
        std::cout << "   Warning: Unrecognized agg style.\nPossible values are Point Cloud, Jacks, Jacks++, Convex Hulls and Alpha Hulls.\nDefaulting to Point Cloud." << std::endl;
        aggStyle = "Point Cloud";
        doPointCloud_(vertices, geomSizes);
      }
      writeFile_(fout, aggStyle, vertices, geomSizes, verticesCoarse, geomSizesCoarse);
      if(doCoarseGraphEdges_)
      {
        string fname = filenameToWrite;
        string cEdgeFile = fname.insert(fname.rfind(".vtu"), "-coarsegraph");
        std::ofstream edgeStream(cEdgeFile.c_str());
        doGraphEdges_(edgeStream, Ac, coarseGraph, false);
      }
      if(doFineGraphEdges_)
      {
        string fname = filenameToWrite;
        string fEdgeFile = fname.insert(fname.rfind(".vtu"), "-finegraph");
        std::ofstream edgeStream(fEdgeFile.c_str());
        doGraphEdges_(edgeStream, Amat, fineGraph, true);
      }
      if(myRank == 0 && pL.get<bool>("aggregation: output file: build colormap"))
      {
        buildColormap_();
      }
    }
    fout.close();
#endif
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doPointCloud_(std::vector<int>& vertices, std::vector<int>& geomSize) const
  {
    /*vertices.reserve(vertices.size() + numNodes_);
    geomSize.reserve(geomSize.size() + numNodes_);
    for(int i = 0; i < numNodes_; i++)
    {
      vertices.push_back(i);
      geomSize.push_back(1);
    }*/
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doJacks_(std::vector<int>& vertices, std::vector<int>& geomSizes) const
  {
    //For each aggregate, find the root node then connect it with the other nodes in the aggregate
    //Doesn't matter the order, as long as all edges are found.
    /*vertices.reserve(vertices.size() + 3 * (numNodes_ - numAggs_));
    geomSizes.reserve(vertices.size() + 2 * (numNodes_ - numAggs_));
    int root = 0;
    for(int i = 0; i < numAggs_; i++) //TODO: Replace this O(n^2) with a better way
    {
      while(!isRoot_[root])
        root++;
      int numInAggFound = 0;
      for(int j = 0; j < numNodes_; j++)
      {
        if(j == root) //don't make a connection from the root to itself
        {
          numInAggFound++;
          continue;
        }
        if(vertex2AggId_[root] == vertex2AggId_[j])
        {
          vertices.push_back(root);
          vertices.push_back(j);
          geomSizes.push_back(2);
          //Also draw the free endpoint explicitly for the current line
          vertices.push_back(j);
          geomSizes.push_back(1);
          numInAggFound++;
          if(numInAggFound == aggSizes_[vertex2AggId_[root]]) //don't spend more time looking if done with that root
            break;
        }
      }
      root++; //get set up to look for the next root
    }*/
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeFile_(std::ofstream& fout, std::string styleName, std::vector<int>& vertices, std::vector<int>& geomSizes, std::vector<int>& verticesCoarse, std::vector<int>& geomSizesCoarse) const
  {
#if 0
    using namespace std;
    vector<int> uniqueFine = makeUnique_(vertices);
    vector<int> uniqueCoarse = makeUnique_(verticesCoarse);
    string indent = "      ";
    fout << "<!--" << styleName << " Aggregates Visualization-->" << endl;
    fout << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << endl;
    fout << "  <UnstructuredGrid>" << endl;
    fout << "    <Piece NumberOfPoints=\"" << uniqueFine.size() + uniqueCoarse.size() << "\" NumberOfCells=\"" << geomSizes.size() + geomSizesCoarse.size() << "\">" << endl;
    fout << "      <PointData Scalars=\"Node Aggregate Processor\">" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Node\" format=\"ascii\">" << endl;
    indent = "          ";
    fout << indent;
    for(int i = 0; i < int(uniqueFine.size()); i++)
    {
      fout << nodeMap_->getGlobalElement(uniqueFine[i]) << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    if(doCoarseGraphEdges_)
    {
      for(int i = 0; i < int(uniqueCoarse.size()); i++)
      {
        fout << nodeMapCoarse_->getGlobalElement(uniqueCoarse[i]) << " ";
        if((i + uniqueFine.size()) % 10 == 9)  //pretty formatting, pick up in the line where other list left off
          fout << endl << indent;
      }
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Aggregate\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(uniqueFine.size()); i++)
    {
      fout << aggsOffset_ + vertex2AggId_[uniqueFine[i]] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    if(doCoarseGraphEdges_)
    {
      for(int i = 0; i < int(uniqueCoarse.size()); i++)
      {
        fout << "-1 ";
        if((i + uniqueFine.size()) % 10 == 9)
          fout << endl << indent;
      }
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(uniqueFine.size() + uniqueCoarse.size()); i++)
    {
      fout << myRank_ << " ";
      if(i % 20 == 19)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </PointData>" << endl;
    fout << "      <Points>" << endl;
    fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(uniqueFine.size()); i++)
    {
      fout << xCoords_[uniqueFine[i]] << " " << yCoords_[uniqueFine[i]] << " ";
      if(dims_ == 2)
        fout << "0 ";
      else
        fout << zCoords_[uniqueFine[i]] << " ";
      if(i % 3 == 2)
        fout << endl << indent;
    }
    for(int i = 0; i < int(uniqueCoarse.size()); i++)
    {
      fout << cx_[uniqueCoarse[i]] << " " << cy_[uniqueCoarse[i]] << " ";
      if(dims_ == 2)
        fout << "0 ";
      else
        fout << cz_[uniqueCoarse[i]] << " ";
      if((i + uniqueFine.size()) % 3 == 2)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Points>" << endl;
    fout << "      <Cells>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(vertices.size()); i++)
    {
      fout << vertices[i] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    for(int i = 0; i < int(verticesCoarse.size()); i++)
    {
      fout << verticesCoarse[i] + uniqueFine.size() << " ";
      if((i + vertices.size()) % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    fout << indent;
    int accum = 0;
    for(int i = 0; i < int(geomSizes.size()); i++)
    {
      accum += geomSizes[i];
      fout << accum << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    for(int i = 0; i < int(geomSizesCoarse.size()); i++)
    {
      accum += geomSizesCoarse[i];
      fout << accum << " ";
      if((i + geomSizes.size()) % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(geomSizes.size()); i++)
    {
      switch(geomSizes[i])
      {
        case 1:
          fout << "1 "; //Point
          break;
        case 2:
          fout << "3 "; //Line
          break;
        case 3:
          fout << "5 "; //Triangle
          break;
        default:
          fout << "7 "; //Polygon
      }
      if(i % 30 == 29)
        fout << endl << indent;
    }
    for(int i = 0; i < int(geomSizesCoarse.size()); i++)
    {
      switch(geomSizesCoarse[i])
      {
        case 1:
          fout << "1 "; //Point
          break;
        case 2:
          fout << "3 "; //Line
          break;
        case 3:
          fout << "5 "; //Triangle
          break;
        default:
          fout << "7 "; //Polygon
      }
      if((i + geomSizes.size()) % 30 == 29)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Cells>" << endl;
    fout << "    </Piece>" << endl;
    fout << "  </UnstructuredGrid>" << endl;
    fout << "</VTKFile>" << endl;
#endif
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::vector<int> CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::makeUnique(std::vector<int>& vertices) const
  {
    using namespace std;
    vector<int> uniqueNodes = vertices;
    sort(uniqueNodes.begin(), uniqueNodes.end());
    vector<int>::iterator newUniqueFineEnd = unique(uniqueNodes.begin(), uniqueNodes.end());
    uniqueNodes.erase(newUniqueFineEnd, uniqueNodes.end());
    //uniqueNodes is now a sorted list of the nodes whose info actually goes in file
    //Now replace values in vertices with locations of the old values in uniqueFine
    for(int i = 0; i < int(vertices.size()); i++)
    {
      int lo = 0;
      int hi = uniqueNodes.size() - 1;
      int mid = 0;
      int search = vertices[i];
      while(lo <= hi)
      {
        mid = lo + (hi - lo) / 2;
        if(uniqueNodes[mid] == search)
          break;
        else if(uniqueNodes[mid] > search)
          hi = mid - 1;
        else
          lo = mid + 1;
      }
      if(uniqueNodes[mid] != search)
        throw runtime_error("Issue in makeUnique_() - a point wasn't found in list.");
      vertices[i] = mid;
    }
    return uniqueNodes;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceAll(std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) const {
    while(1) {
      const int pos = result.find(replaceWhat);
      if (pos == -1)
        break;
      result.replace(pos, replaceWhat.size(), replaceWithWhat);
    }
    return result;
  }
} // namespace MueLu

#endif /* MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_ */
