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
    validParamList->set< RCP<const FactoryBase> >("P",           Teuchos::null, "Prolongator factory. The user has to declare either P or Ptent but not both at the same time.");
    validParamList->set< RCP<const FactoryBase> >("Ptent",       Teuchos::null, "Tentative prolongator factory. The user has to declare either P or Ptent as input but not both at the same time");
    validParamList->set< RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory for Coordinates.");
    //validParamList->set< RCP<const FactoryBase> >("Graph", Teuchos::null, "Factory for Graph.");
    //validParamList->set< RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Factory for Aggregates.");
    //validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Factory for UnAmalgamationInfo.");
    //validParamList->set< RCP<const FactoryBase> >("DofsPerNode", Teuchos::null, "Factory for DofsPerNode.");

    // New-style master list options (here are same defaults as in masterList.xml)
    validParamList->set< std::string >           ("visualization: output filename",                    "viz%LEVELID",                    "filename for VTK-style visualization output");
    validParamList->set<std::string>             ("visualization: style", "Point Cloud", "style of aggregate visualization for VTK output. Can be 'Point Cloud', 'Jacks', 'Convex Hulls'");
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



    //Input(fineLevel, "A");
    //Input(fineLevel, "Graph");

    Input(fineLevel, "Coordinates");

    const ParameterList & pL = GetParameterList();
    TEUCHOS_TEST_FOR_EXCEPTION(pL.isParameter("P") && GetFactory("P") != Teuchos::null &&
                               pL.isParameter("Ptent") && GetFactory("Ptent") != Teuchos::null, Exceptions::RuntimeError,
                "You must not declare both P and Ptent. Use only once for visualization.");
    TEUCHOS_TEST_FOR_EXCEPTION(GetFactory("P") == Teuchos::null && GetFactory("Ptent") == Teuchos::null, Exceptions::RuntimeError,
                "You have to either declare P or Ptent for visualization, but not both.");

    if (GetFactory("P") != Teuchos::null && GetFactory("Ptent") == Teuchos::null)
      Input(coarseLevel, "P");
    else if (GetFactory("Ptent") != Teuchos::null && GetFactory("P") == Teuchos::null)
      Input(coarseLevel, "Ptent");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {

    RCP<Matrix> P = Teuchos::null;
    const ParameterList & pL = GetParameterList();
    if (GetFactory("P") != Teuchos::null && GetFactory("Ptent") == Teuchos::null)
      P = Get< RCP<Matrix> >(coarseLevel, "P");
    if (GetFactory("Ptent") != Teuchos::null && GetFactory("P") == Teuchos::null)
      P = Get< RCP<Matrix> >(coarseLevel, "Ptent");

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

    // TODO add support for overlapping aggregates
    TEUCHOS_TEST_FOR_EXCEPTION(strDomainMap->getNodeNumElements() != P->getColMap()->getNodeNumElements(), Exceptions::RuntimeError,
                                               "CoarseningVisualization only supports non-overlapping transfers");

    // number of local "aggregates"
    LocalOrdinal numLocalAggs = strDomainMap->getNodeNumElements() / columnsPerNode;
    std::vector< std::set<LocalOrdinal> > localAggs(numLocalAggs);

    // do loop over all local rows of prolongator and extract column information
    for (LO row = 0; row < Teuchos::as<LO>(P->getRowMap()->getNodeNumElements()); ++row) {
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
        double mydistance = 0.0;
        mydistance += tempx*tempx;
        mydistance += tempy*tempy;
        mydistance += tempz*tempz;
        mydistance = sqrt(mydistance);
        if (mydistance <= minDistance) {
          minDistance = mydistance;
          rootCandidate = *it;
        }
      }

      isRoot[rootCandidate] = true;
    }


    std::vector<LocalOrdinal> vertices;
    std::vector<LocalOrdinal> geomSize;

    std::string aggStyle = pL.get<std::string>("aggregation: output file: agg style");
    if(aggStyle == "Point Cloud")
      doPointCloud(vertices, geomSize, numLocalAggs, numFineNodes);
    else if(aggStyle == "Jacks")
      doJacks(vertices, geomSize, numLocalAggs, numFineNodes, isRoot, vertex2AggId);
    else if(aggStyle == "Convex Hulls") {
      if(coords->getNumVectors() == 3)
        doConvexHulls3D(vertices, geomSize, numLocalAggs, numFineNodes, isRoot, vertex2AggId, xCoords, yCoords, zCoords);
      else if(coords->getNumVectors() == 2)
        doConvexHulls2D(vertices, geomSize, numLocalAggs, numFineNodes, isRoot, vertex2AggId, xCoords, yCoords, zCoords);
    }
    else
    {
      GetOStream(Warnings0) << "   Warning: Unrecognized agg style.\nPossible values are Point Cloud, Jacks, Convex Hulls.\nDefaulting to Point Cloud." << std::endl;
      aggStyle = "Point Cloud";
      doPointCloud(vertices, geomSize, numLocalAggs, numFineNodes);
    }


    // end replace this

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
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doPointCloud(std::vector<int>& vertices, std::vector<int>& geomSizes, LO numLocalAggs, LO numFineNodes) {
    vertices.reserve(numFineNodes);
    geomSizes.reserve(numFineNodes);
    for(LocalOrdinal i = 0; i < numFineNodes; i++)
    {
      vertices.push_back(i);
      geomSizes.push_back(1);
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doJacks(std::vector<int>& vertices, std::vector<int>& geomSizes, LO numLocalAggs, LO numFineNodes, const std::vector<bool>& isRoot, const std::vector<LO>& vertex2AggId) {
    //For each aggregate, find the root node then connect it with the other nodes in the aggregate
    //Doesn't matter the order, as long as all edges are found.
    vertices.reserve(vertices.size() + 3 * (numFineNodes - numLocalAggs));
    geomSizes.reserve(vertices.size() + 2 * (numFineNodes - numLocalAggs));
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
          geomSizes.push_back(2);
          //Also draw the free endpoint explicitly for the current line
          vertices.push_back(j);
          geomSizes.push_back(1);
          numInAggFound++;
          //if(numInAggFound == aggSizes_[vertex2AggId_[root]]) //don't spend more time looking if done with that root
          //  break;
        }
      }
      root++; //get set up to look for the next root
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doConvexHulls2D(std::vector<int>& vertices, std::vector<int>& geomSizes, LO numLocalAggs, LO numFineNodes, const std::vector<bool>& isRoot, const std::vector<LO>& vertex2AggId, const Teuchos::ArrayRCP<const double>& xCoords, const Teuchos::ArrayRCP<const double>& yCoords, const Teuchos::ArrayRCP<const double>& zCoords) {
    using namespace std;
    for(int agg = 0; agg < numLocalAggs; agg++)
    {
      list<int> aggNodes;
      for(int i = 0; i < numFineNodes; i++)
      {
        if(vertex2AggId[i] == agg)
          aggNodes.push_back(i);
      }
      //have a list of nodes in the aggregate
      TEUCHOS_TEST_FOR_EXCEPTION(aggNodes.size() == 0, Exceptions::RuntimeError,
               "CoarseningVisualization::doConvexHulls2D: aggregate contains zero nodes!");
      if(aggNodes.size() == 1)
      {
        vertices.push_back(aggNodes.front());
        geomSizes.push_back(1);
        continue;
      }
      if(aggNodes.size() == 2)
      {
        vertices.push_back(aggNodes.front());
        vertices.push_back(aggNodes.back());
        continue;
      }
      //check if all points are collinear, need to explicitly draw a line in that case.
      bool collinear = true; //assume true at first, if a segment not parallel to others then clear
      {
        list<int>::iterator it = aggNodes.begin();
        myVec3 firstPoint(xCoords[*it], yCoords[*it], 0);
        it++;
        myVec3 secondPoint(xCoords[*it], yCoords[*it], 0);
        it++;  //it now points to third node in the aggregate
        myVec3 norm1(-(secondPoint.y - firstPoint.y), secondPoint.x - firstPoint.x, 0);
        do
        {
          myVec3 thisNorm(yCoords[*it] - firstPoint.y, firstPoint.x - xCoords[*it], 0);
          //rotate one of the vectors by 90 degrees so that dot product is 0 if the two are parallel
          double temp = thisNorm.x;
          thisNorm.x = thisNorm.y;
          thisNorm.y = temp;
          double comp = dotProduct(norm1, thisNorm);
          if(-1e-8 > comp || comp > 1e-8)
          {
            collinear = false;
            break;
          }
          it++;
        }
        while(it != aggNodes.end());
      }
      if(collinear) {
        //find the most distant two points in the plane and use as endpoints of line representing agg
        list<int>::iterator min = aggNodes.begin();    //min X then min Y where x is a tie
        list<int>::iterator max = aggNodes.begin(); //max X then max Y where x is a tie
        for(list<int>::iterator it = ++aggNodes.begin(); it != aggNodes.end(); it++) {
          if(xCoords[*it] < xCoords[*min])
            min = it;
          else if(xCoords[*it] == xCoords[*min]) {
            if(yCoords[*it] < yCoords[*min]) min = it;
          }
          if(xCoords[*it] > xCoords[*max]) max = it;
          else if(xCoords[*it] == xCoords[*max]) {
            if(yCoords[*it] > yCoords[*max]) max = it;
          }
        }
        //Just set up a line between nodes *min and *max
        vertices.push_back(*min);
        vertices.push_back(*max);
        geomSizes.push_back(2);
        continue; //jump to next aggregate in loop
      }
      list<int>::iterator min = aggNodes.begin();
      for(list<int>::iterator it = ++aggNodes.begin(); it != aggNodes.end(); it++) {
        if(xCoords[*it] < xCoords[*min]) min = it;
        else if(xCoords[*it] == xCoords[*min]) {
          if(yCoords[*it] < yCoords[*min]) min = it;
        }
      }
      //this is the most common case: at least 3 nodes making up a polygon with positive area
      //do Jarvis march on these points to compute convex hull
      //start with "min" point
      int thisHullSize = 1;
      vertices.push_back(*min);
      //repeatedly sweep through aggNodes (with a pivot at thisHull.back()) counterclockwise. Connect up with last node found this way ("left-most"). If there is a tie for the left-most node, take the most distant one from thisHull.back().
      bool includeMin = false;
      while(1) {
        list<int>::iterator leftMost = aggNodes.begin();
        if(!includeMin && leftMost == min) {
          leftMost++;
        }
        list<int>::iterator it = leftMost;
        it++;
        while(it != aggNodes.end()) {
          if(it == min && !includeMin) { //don't compare to min on very first sweep
            it++;
            continue;
          }
          //see if it is in front of line containing nodes thisHull.back() and leftMost
          //first get the left normal of leftMost - thisHull.back() (<dy, -dx>)
          myVec3 testNorm(yCoords[*leftMost] - yCoords[vertices.back()], -(xCoords[*leftMost] - xCoords[vertices.back()]), 0);
          //now dot testNorm with *it - leftMost. If dot is positive, leftMost becomes it. If dot is zero, take one further from thisHull.back().
          myVec3 itVec(xCoords[*it] - xCoords[*leftMost], yCoords[*it] - yCoords[*leftMost], 0);
          double dotProd = dotProduct(testNorm, itVec);
          if(-1e-10 < dotProd && dotProd < 1e-10) {
            //thisHull.back(), it and leftMost are collinear.
            //Just sum the differences in x and differences in y for each and compare to get further one, don't need distance formula
            double itDist = fabs(xCoords[*it] - xCoords[vertices.back()]) + fabs(yCoords[*it] - yCoords[vertices.back()]);
            double leftMostDist = fabs(xCoords[*leftMost] - xCoords[vertices.back()]) + fabs(yCoords[*leftMost] - yCoords[vertices.back()]);
            if(itDist > leftMostDist) leftMost = it;
          }
          else if(dotProd > 0)
            leftMost = it;
          it++;
        }
        //if leftMost is min, then the loop is complete.
        if(leftMost == min) {
          geomSizes.push_back(thisHullSize);
          break; //this goes to the next aggregate
        }
        //add leftMost to thisHull, and remove it from aggNodes so later searches are faster
        vertices.push_back(*leftMost);
        thisHullSize++;
        aggNodes.erase(leftMost);
        includeMin = true; //have found second point (the one after min) so now include min in the searches
      }
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doConvexHulls3D(std::vector<int>& vertices, std::vector<int>& geomSizes, LO numLocalAggs, LO numFineNodes, const std::vector<bool>& isRoot, const std::vector<LO>& vertex2AggId, const Teuchos::ArrayRCP<const double>& xCoords, const Teuchos::ArrayRCP<const double>& yCoords, const Teuchos::ArrayRCP<const double>& zCoords) {
    //Use 3D quickhull algo.
    //Vector of node indices representing triangle vertices
    //Note: Calculate the hulls first since will only include point data for points in the hulls
    //Effectively the size() of vertIndices after each hull is added to it
    typedef std::list<int>::iterator Iter;
    for(int agg = 0; agg < numLocalAggs; agg++)
    {
      std::list<int> aggNodes; //At first, list of all nodes in the aggregate. As nodes are enclosed or included by/in hull, remove them
      for(int i = 0; i < numFineNodes; i++)
      {
        if(vertex2AggId[i] == agg)
          aggNodes.push_back(i);
      }
      //First, check anomalous cases
      TEUCHOS_TEST_FOR_EXCEPTION(aggNodes.size() == 0, Exceptions::RuntimeError,
               "CoarseningVisualization::doConvexHulls3D: aggregate contains zero nodes!");
      if(aggNodes.size() == 1)
      {
        vertices.push_back(aggNodes.front());
        geomSizes.push_back(1);
        continue;
      }
      else if(aggNodes.size() == 2)
      {
        vertices.push_back(aggNodes.front());
        vertices.push_back(aggNodes.back());
        geomSizes.push_back(2);
        continue;
      }
      //check for collinearity
      bool areCollinear = true;
      {
        Iter colCheck = aggNodes.begin();
        int firstNode = *colCheck;
        colCheck++;
        int secondNode = *colCheck;
        colCheck++;
        double colCheckDX = xCoords[firstNode] - xCoords[secondNode];
        double colCheckDY = yCoords[firstNode] - yCoords[secondNode];
        double colCheckDZ = zCoords[firstNode] - zCoords[secondNode];
        myVec3 collinearVec(colCheckDX, colCheckDY, colCheckDZ);
        //normalize collinearVec
        double collinearVecMag = mymagnitude(collinearVec);
        collinearVec.x /= collinearVecMag;
        collinearVec.y /= collinearVecMag;
        collinearVec.z /= collinearVecMag;
        myVec3 firstPoint(xCoords[aggNodes.front()], yCoords[aggNodes.front()], zCoords[aggNodes.front()]);
        for(Iter it = colCheck; it != aggNodes.end(); it++)
        {
          myVec3 thisVec(xCoords[*it], yCoords[*it], zCoords[*it]);
          myVec3 vecDiff = vecSubtract(thisVec, firstPoint);
          //normalize vecDiff so that it can be directly compared to collinearVec
          double vecDiffMag = mymagnitude(vecDiff);
          vecDiff.x /= vecDiffMag;
          vecDiff.y /= vecDiffMag;
          vecDiff.z /= vecDiffMag;
          //compare x, y and z separately and give some slack for rounding error
          myVec3 compare = vecSubtract(vecDiff, collinearVec);
          if(compare.x < -1e-10 || compare.x > 1e-10 || compare.y < -1e-10 || compare.y > 1e-10 || compare.z < -1e-10 || compare.z > 1e-10)
          {
            areCollinear = false;
            break;
          }
        }
      }
      if(areCollinear) {
        //find the endpoints of segment describing all the points
        //compare x, if tie compare y, if tie compare z
        Iter min = aggNodes.begin();
        Iter max = aggNodes.begin();
        Iter it = ++aggNodes.begin();
        for(; it != aggNodes.end(); it++) {
          if(xCoords[*it] < xCoords[*min]) min = it;
          else if(xCoords[*it] == xCoords[*min]) {
            if(yCoords[*it] < yCoords[*min]) min = it;
            else if(yCoords[*it] == yCoords[*min]) {
              if(zCoords[*it] < zCoords[*min]) min = it;
            }
          }
          if(xCoords[*it] > xCoords[*max]) max = it;
          else if(xCoords[*it] == xCoords[*max]) {
            if(yCoords[*it] > yCoords[*max]) max = it;
            else if(yCoords[*it] == yCoords[*max]) {
              if(zCoords[*it] > zCoords[*max]) max = it;
            }
          }
        }
        vertices.push_back(*min);
        vertices.push_back(*max);
        geomSizes.push_back(2);
        continue;
      }
      Iter exIt = aggNodes.begin(); //iterator to be used for searching for min/max x/y/z
      int extremeSix[] = {*exIt, *exIt, *exIt, *exIt, *exIt, *exIt}; //nodes with minimumX, maxX, minY ...
      exIt++;
      for(; exIt != aggNodes.end(); exIt++) {
        Iter it = exIt;
        if(xCoords[*it] < xCoords[extremeSix[0]] ||
          (xCoords[*it] == xCoords[extremeSix[0]] && yCoords[*it] < yCoords[extremeSix[0]]) ||
          (xCoords[*it] == xCoords[extremeSix[0]] && yCoords[*it] == yCoords[extremeSix[0]] && zCoords[*it] < zCoords[extremeSix[0]]))
            extremeSix[0] = *it;
        if(xCoords[*it] > xCoords[extremeSix[1]] ||
          (xCoords[*it] == xCoords[extremeSix[1]] && yCoords[*it] > yCoords[extremeSix[1]]) ||
          (xCoords[*it] == xCoords[extremeSix[1]] && yCoords[*it] == yCoords[extremeSix[1]] && zCoords[*it] > zCoords[extremeSix[1]]))
            extremeSix[1] = *it;
        if(yCoords[*it] < yCoords[extremeSix[2]] ||
          (yCoords[*it] == yCoords[extremeSix[2]] && zCoords[*it] < zCoords[extremeSix[2]]) ||
          (yCoords[*it] == yCoords[extremeSix[2]] && zCoords[*it] == zCoords[extremeSix[2]] && xCoords[*it] < xCoords[extremeSix[2]]))
            extremeSix[2] = *it;
        if(yCoords[*it] > yCoords[extremeSix[3]] ||
          (yCoords[*it] == yCoords[extremeSix[3]] && zCoords[*it] > zCoords[extremeSix[3]]) ||
          (yCoords[*it] == yCoords[extremeSix[3]] && zCoords[*it] == zCoords[extremeSix[3]] && xCoords[*it] > xCoords[extremeSix[3]]))
            extremeSix[3] = *it;
        if(zCoords[*it] < zCoords[extremeSix[4]] ||
          (zCoords[*it] == zCoords[extremeSix[4]] && xCoords[*it] < xCoords[extremeSix[4]]) ||
          (zCoords[*it] == zCoords[extremeSix[4]] && xCoords[*it] == xCoords[extremeSix[4]] && yCoords[*it] < yCoords[extremeSix[4]]))
            extremeSix[4] = *it;
        if(zCoords[*it] > zCoords[extremeSix[5]] ||
          (zCoords[*it] == zCoords[extremeSix[5]] && xCoords[*it] > xCoords[extremeSix[5]]) ||
          (zCoords[*it] == zCoords[extremeSix[5]] && xCoords[*it] == xCoords[extremeSix[5]] && yCoords[*it] > zCoords[extremeSix[5]]))
            extremeSix[5] = *it;
      }
      myVec3 extremeVectors[6];
      for(int i = 0; i < 6; i++) {
        myVec3 thisExtremeVec(xCoords[extremeSix[i]], yCoords[extremeSix[i]], zCoords[extremeSix[i]]);
        extremeVectors[i] = thisExtremeVec;
      }
      double maxDist = 0;
      int base1 = 0; //ints from 0-5: which pair out of the 6 extreme points are the most distant? (indices in extremeSix and extremeVectors)
      int base2 = 0;
      for(int i = 0; i < 5; i++) {
        for(int j = i + 1; j < 6; j++) {
          double thisDist = distance(extremeVectors[i], extremeVectors[j]);
          if(thisDist > maxDist) {
            maxDist = thisDist;
            base1 = i;
            base2 = j;
          }
        }
      }
      std::list<myTriangle> hullBuilding;    //each Triangle is a triplet of nodes (int IDs) that form a triangle
      //remove base1 and base2 iters from aggNodes, they are known to be in the hull
      aggNodes.remove(extremeSix[base1]);
      aggNodes.remove(extremeSix[base2]);
      //extremeSix[base1] and [base2] still have the vec3_ representation
      myTriangle tri;
      tri.v1 = extremeSix[base1];
      tri.v2 = extremeSix[base2];
      //Now find the point that is furthest away from the line between base1 and base2
      maxDist = 0;
      //need the vectors to do "quadruple product" formula
      myVec3 b1 = extremeVectors[base1];
      myVec3 b2 = extremeVectors[base2];
      Iter thirdNode;
      for(Iter node = aggNodes.begin(); node != aggNodes.end(); node++) {
        myVec3 nodePos(xCoords[*node], yCoords[*node], zCoords[*node]);
        double dist = mymagnitude(crossProduct(vecSubtract(nodePos, b1), vecSubtract(nodePos, b2))) / mymagnitude(vecSubtract(b2, b1));
        if(dist > maxDist) {
          maxDist = dist;
          thirdNode = node;
        }
      }
      //Now know the last node in the first triangle
      tri.v3 = *thirdNode;
      hullBuilding.push_back(tri);
      myVec3 b3(xCoords[*thirdNode], yCoords[*thirdNode], zCoords[*thirdNode]);
      aggNodes.erase(thirdNode);
      //Find the fourth node (most distant from triangle) to form tetrahedron
      maxDist = 0;
      int fourthVertex = -1;
      for(Iter node = aggNodes.begin(); node != aggNodes.end(); node++) {
        myVec3 thisNode(xCoords[*node], yCoords[*node], zCoords[*node]);
        double nodeDist = pointDistFromTri(thisNode, b1, b2, b3);
        if(nodeDist > maxDist) {
          maxDist = nodeDist;
          fourthVertex = *node;
        }
      }
      aggNodes.remove(fourthVertex);
      myVec3 b4(xCoords[fourthVertex], yCoords[fourthVertex], zCoords[fourthVertex]);
      //Add three new triangles to hullBuilding to form the first tetrahedron
      //use tri to hold the triangle info temporarily before being added to list
      tri = hullBuilding.front();
      tri.v1 = fourthVertex;
      hullBuilding.push_back(tri);
      tri = hullBuilding.front();
      tri.v2 = fourthVertex;
      hullBuilding.push_back(tri);
      tri = hullBuilding.front();
      tri.v3 = fourthVertex;
      hullBuilding.push_back(tri);
      //now orient all four triangles so that the vertices are oriented clockwise (so getNorm_ points outward for each)
      myVec3 barycenter((b1.x + b2.x + b3.x + b4.x) / 4.0, (b1.y + b2.y + b3.y + b4.y) / 4.0, (b1.z + b2.z + b3.z + b4.z) / 4.0);
      for(std::list<myTriangle>::iterator tetTri = hullBuilding.begin(); tetTri != hullBuilding.end(); tetTri++) {
        myVec3 triVert1(xCoords[tetTri->v1], yCoords[tetTri->v1], zCoords[tetTri->v1]);
        myVec3 triVert2(xCoords[tetTri->v2], yCoords[tetTri->v2], zCoords[tetTri->v2]);
        myVec3 triVert3(xCoords[tetTri->v3], yCoords[tetTri->v3], zCoords[tetTri->v3]);
        myVec3 trinorm = getNorm(triVert1, triVert2, triVert3);
        myVec3 ptInPlane = tetTri == hullBuilding.begin() ? b1 : b4; //first triangle definitely has b1 in it, other three definitely have b4
        if(isInFront(barycenter, ptInPlane, trinorm)) {
          //don't want the faces of the tetrahedron to face inwards (towards barycenter) so reverse orientation
          //by swapping two vertices
          int temp = tetTri->v1;
          tetTri->v1 = tetTri->v2;
          tetTri->v2 = temp;
        }
      }
      //now, have starting polyhedron in hullBuilding (all faces are facing outwards according to getNorm_) and remaining nodes to process are in aggNodes
      //recursively, for each triangle, make a list of the points that are 'in front' of the triangle. Find the point with the maximum distance from the triangle.
      //Add three new triangles, each sharing one edge with the original triangle but now with the most distant point as a vertex. Remove the most distant point from
      //the list of all points that need to be processed. Also from that list remove all points that are in front of the original triangle but not in front of all three
      //new triangles, since they are now enclosed in the hull.
      //Construct point lists for each face of the tetrahedron individually.
      myVec3 trinorms[4]; //normals to the four tetrahedron faces, now oriented outwards
      int index = 0;
      for(std::list<myTriangle>::iterator tetTri = hullBuilding.begin(); tetTri != hullBuilding.end(); tetTri++) {
        myVec3 triVert1(xCoords[tetTri->v1], yCoords[tetTri->v1], zCoords[tetTri->v1]);
        myVec3 triVert2(xCoords[tetTri->v2], yCoords[tetTri->v2], zCoords[tetTri->v2]);
        myVec3 triVert3(xCoords[tetTri->v3], yCoords[tetTri->v3], zCoords[tetTri->v3]);
        trinorms[index] = getNorm(triVert1, triVert2, triVert3);
        index++;
      }
      std::list<int> startPoints1;
      std::list<int> startPoints2;
      std::list<int> startPoints3;
      std::list<int> startPoints4;
      //scope this so that 'point' is not in function scope
      {
        Iter point = aggNodes.begin();
        while(!aggNodes.empty()) { //this removes points one at a time as they are put in startPointsN or are already done
          myVec3 pointVec(xCoords[*point], yCoords[*point], zCoords[*point]);
          //Note: Because of the way the tetrahedron faces are constructed above,
          //face 1 definitely contains b1 and faces 2-4 definitely contain b4.
          if(isInFront(pointVec, b1, trinorms[0])) {
            startPoints1.push_back(*point);
            point = aggNodes.erase(point);
          } else if(isInFront(pointVec, b4, trinorms[1])) {
            startPoints2.push_back(*point);
            point = aggNodes.erase(point);
          } else if(isInFront(pointVec, b4, trinorms[2])) {
            startPoints3.push_back(*point);
            point = aggNodes.erase(point);
          } else if(isInFront(pointVec, b4, trinorms[3])) {
            startPoints4.push_back(*point);
            point = aggNodes.erase(point);
          } else {
            point = aggNodes.erase(point); //points here are already inside tetrahedron.
          }
        }
        //Call processTriangle_ for each triangle in the initial tetrahedron, one at a time.
      }
      typedef std::list<myTriangle>::iterator TriIter;
      TriIter firstTri = hullBuilding.begin();
      myTriangle start1 = *firstTri;
      firstTri++;
      myTriangle start2 = *firstTri;
      firstTri++;
      myTriangle start3 = *firstTri;
      firstTri++;
      myTriangle start4 = *firstTri;
      //kick off depth-first recursive filling of hullBuilding list with all triangles in the convex hull
      if(!startPoints1.empty())
        processTriangle(hullBuilding, start1, startPoints1, barycenter, xCoords, yCoords, zCoords);
      if(!startPoints2.empty())
        processTriangle(hullBuilding, start2, startPoints2, barycenter, xCoords, yCoords, zCoords);
      if(!startPoints3.empty())
        processTriangle(hullBuilding, start3, startPoints3, barycenter, xCoords, yCoords, zCoords);
      if(!startPoints4.empty())
        processTriangle(hullBuilding, start4, startPoints4, barycenter, xCoords, yCoords, zCoords);
      //hullBuilding now has all triangles that make up this hull.
      //Dump hullBuilding info into the list of all triangles for the scene.
      vertices.reserve(vertices.size() + 3 * hullBuilding.size());
      for(TriIter hullTri = hullBuilding.begin(); hullTri != hullBuilding.end(); hullTri++) {
        vertices.push_back(hullTri->v1);
        vertices.push_back(hullTri->v2);
        vertices.push_back(hullTri->v3);
        geomSizes.push_back(3);
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  myVec3 CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::crossProduct(myVec3 v1, myVec3 v2)
  {
    return myVec3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dotProduct(myVec3 v1, myVec3 v2)
  {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isInFront(myVec3 point, myVec3 inPlane, myVec3 n)
  {
    myVec3 rel(point.x - inPlane.x, point.y - inPlane.y, point.z - inPlane.z); //position of the point relative to the plane
    return dotProduct(rel, n) > 1e-12 ? true : false;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mymagnitude(myVec3 vec)
  {
    return sqrt(dotProduct(vec, vec));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::distance(myVec3 p1, myVec3 p2)
  {
    return mymagnitude(vecSubtract(p1, p2));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  myVec3 CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::vecSubtract(myVec3 v1, myVec3 v2)
  {
    return myVec3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  myVec3 CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNorm(myVec3 v1, myVec3 v2, myVec3 v3) //normal to face of triangle (will be outward rel. to polyhedron) (v1, v2, v3 are in CCW order when normal is toward viewpoint)
  {
    return crossProduct(vecSubtract(v2, v1), vecSubtract(v3, v1));
  }

  //get minimum distance from 'point' to plane containing v1, v2, v3 (or the triangle with v1, v2, v3 as vertices)
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::pointDistFromTri(myVec3 point, myVec3 v1, myVec3 v2, myVec3 v3)
  {
    using namespace std;
    myVec3 norm = getNorm(v1, v2, v3);
    //must normalize the normal vector
    double normScl = mymagnitude(norm);
    norm.x /= normScl;
    norm.y /= normScl;
    norm.z /= normScl;
    double rv = fabs(dotProduct(norm, vecSubtract(point, v1)));
    return rv;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::vector<myTriangle> CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::processTriangle(std::list<myTriangle>& tris, myTriangle tri, std::list<int>& pointsInFront, myVec3& barycenter, const Teuchos::ArrayRCP<const double>& xCoords, const Teuchos::ArrayRCP<const double>& yCoords, const Teuchos::ArrayRCP<const double>& zCoords) {
    //*tri is in the tris list, and is the triangle to process here. tris is a complete list of all triangles in the hull so far. pointsInFront is only a list of the nodes in front of tri. Need coords also.
    //precondition: each triangle is already oriented so that getNorm_(v1, v2, v3) points outward (away from interior of hull)
    //First find the point furthest from triangle.
    using namespace std;
    typedef std::list<int>::iterator Iter;
    typedef std::list<myTriangle>::iterator TriIter;
    typedef list<pair<int, int> >::iterator EdgeIter;
    double maxDist = 0;
    //Need vector representations of triangle's vertices
    myVec3 v1(xCoords[tri.v1], yCoords[tri.v1], zCoords[tri.v1]);
    myVec3 v2(xCoords[tri.v2], yCoords[tri.v2], zCoords[tri.v2]);
    myVec3 v3(xCoords[tri.v3], yCoords[tri.v3], zCoords[tri.v3]);
    myVec3 farPointVec; //useful to have both the point's coordinates and it's position in the list
    Iter farPoint = pointsInFront.begin();
    for(Iter point = pointsInFront.begin(); point != pointsInFront.end(); point++)
    {
      myVec3 pointVec(xCoords[*point], yCoords[*point], zCoords[*point]);
      double dist = pointDistFromTri(pointVec, v1, v2, v3);
      if(dist > maxDist)
      {
        dist = maxDist;
        farPointVec = pointVec;
        farPoint = point;
      }
    }
    //Find all the triangles that the point is in front of (can be more than 1)
    //At the same time, remove them from tris, as every one will be replaced later
    vector<myTriangle> visible; //use a list of iterators so that the underlying object is still in tris
    for(TriIter it = tris.begin(); it != tris.end();)
    {
      myVec3 vec1(xCoords[it->v1], yCoords[it->v1], zCoords[it->v1]);
      myVec3 vec2(xCoords[it->v2], yCoords[it->v2], zCoords[it->v2]);
      myVec3 vec3(xCoords[it->v3], yCoords[it->v3], zCoords[it->v3]);
      myVec3 norm = getNorm(vec1, vec2, vec3);
      if(isInFront(farPointVec, vec1, norm))
      {
        visible.push_back(*it);
        it = tris.erase(it);
      }
      else
        it++;
    }
    //Figure out what triangles need to be destroyed/created
    //First create a list of edges (as std::pair<int, int>, where the two ints are the node endpoints)
    list<pair<int, int> > horizon;
    //For each triangle, add edges to the list iff the edge only appears once in the set of all
    //Have members of horizon have the lower node # first, and the higher one second
    for(vector<myTriangle>::iterator it = visible.begin(); it != visible.end(); it++)
    {
      pair<int, int> e1(it->v1, it->v2);
      pair<int, int> e2(it->v2, it->v3);
      pair<int, int> e3(it->v1, it->v3);
      //"sort" the pair values
      if(e1.first > e1.second)
      {
        int temp = e1.first;
        e1.first = e1.second;
        e1.second = temp;
      }
      if(e2.first > e2.second)
      {
        int temp = e2.first;
        e2.first = e2.second;
        e2.second = temp;
      }
      if(e3.first > e3.second)
      {
        int temp = e3.first;
        e3.first = e3.second;
        e3.second = temp;
      }
      horizon.push_back(e1);
      horizon.push_back(e2);
      horizon.push_back(e3);
    }
    //sort based on lower node first, then higher node (lexicographical ordering built in to pair)
    horizon.sort();
    //Remove all edges from horizon, except those that appear exactly once
    {
      EdgeIter it = horizon.begin();
      while(it != horizon.end())
      {
        int occur = count(horizon.begin(), horizon.end(), *it);
        if(occur > 1)
        {
          pair<int, int> removeVal = *it;
          while(removeVal == *it && !(it == horizon.end()))
            it = horizon.erase(it);
        }
        else
          it++;
      }
    }
    //Now make a list of new triangles being created, each of which take 2 vertices from an edge and one from farPoint
    list<myTriangle> newTris;
    for(EdgeIter it = horizon.begin(); it != horizon.end(); it++)
    {
      myTriangle t(it->first, it->second, *farPoint);
      newTris.push_back(t);
    }
    //Ensure every new triangle is oriented outwards, using the barycenter of the initial tetrahedron
    vector<myTriangle> trisToProcess;
    vector<list<int> > newFrontPoints;
    for(TriIter it = newTris.begin(); it != newTris.end(); it++)
    {
      myVec3 t1(xCoords[it->v1], yCoords[it->v1], zCoords[it->v1]);
      myVec3 t2(xCoords[it->v2], yCoords[it->v2], zCoords[it->v2]);
      myVec3 t3(xCoords[it->v3], yCoords[it->v3], zCoords[it->v3]);
      if(isInFront(barycenter, t1, getNorm(t1, t2, t3)))
      {
        //need to swap two vertices to flip orientation of triangle
        int temp = it->v1;
        myVec3 tempVec = t1;
        it->v1 = it->v2;
        t1 = t2;
        it->v2 = temp;
        t2 = tempVec;
      }
      myVec3 outwardNorm = getNorm(t1, t2, t3); //now definitely points outwards
      //Add the triangle to tris
      tris.push_back(*it);
      trisToProcess.push_back(tris.back());
      //Make a list of the points that are in front of nextToProcess, to be passed in for processing
      list<int> newInFront;
      for(Iter point = pointsInFront.begin(); point != pointsInFront.end();)
      {
        myVec3 pointVec(xCoords[*point], yCoords[*point], zCoords[*point]);
        if(isInFront(pointVec, t1, outwardNorm))
        {
          newInFront.push_back(*point);
          point = pointsInFront.erase(point);
        }
        else
          point++;
      }
      newFrontPoints.push_back(newInFront);
    }
    vector<myTriangle> allRemoved; //list of all invalid iterators that were erased by calls to processmyTriangle below
    for(int i = 0; i < int(trisToProcess.size()); i++)
    {
      if(!newFrontPoints[i].empty())
      {
        //Comparing the 'triangle to process' to the one for this call prevents infinite recursion/stack overflow.
        //TODO: Why was it doing that? Rounding error? Make more robust fix. But this does work for the time being.
        if(find(allRemoved.begin(), allRemoved.end(), trisToProcess[i]) == allRemoved.end() && !(trisToProcess[i] == tri))
        {
          vector<myTriangle> removedList = processTriangle(tris, trisToProcess[i], newFrontPoints[i], barycenter, xCoords, yCoords, zCoords);
          for(int j = 0; j < int(removedList.size()); j++)
            allRemoved.push_back(removedList[j]);
        }
      }
    }
    return visible;
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
