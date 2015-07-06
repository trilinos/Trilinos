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

#include "MueLu_AggregationExportFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"
#include <vector>
#include <list>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <cmath>

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    std::string output_msg = "Output filename template (%TIMESTEP is replaced by \'Output file: time step\' variable,"
        "%ITER is replaced by \'Output file: iter\' variable, %LEVELID is replaced level id, %PROCID is replaced by processor id)";
    std::string output_def = "aggs_level%LEVELID_proc%PROCID.out";

    validParamList->set< RCP<const FactoryBase> >("Aggregates",             Teuchos::null, "Generating factory for aggregates");
    validParamList->set< RCP<const FactoryBase> >("DofsPerNode",            Teuchos::null, "Generating factory for number of dofs per node");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo",     Teuchos::null, "Generating factory for amalgamation");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",            Teuchos::null, "Transfer factory for coordinates");

    // CMS/BMK: Old style factory-only options.  Deprecate me.
    validParamList->set< std::string >           ("Output filename",           output_def, output_msg);
    validParamList->set< int >                   ("Output file: time step",             0, "time step variable for output file name");
    validParamList->set< int >                   ("Output file: iter",                  0, "nonlinear iteration variable for output file name");

    // New-style master list options (here are same defaults as in masterList.xml)
    validParamList->set< std::string >           ("aggregation: output filename",                    "",         "filename for VTK-style visualization output");
    validParamList->set< int >                   ("aggregation: output file: time step",              0,          "time step variable for output file name");// Remove me?
    validParamList->set< int >                   ("aggregation: output file: iter",                  0,          "nonlinear iteration variable for output file name");//Remove me?
    validParamList->set<std::string>             ("aggregation: output file: agg style",             "Point Cloud",         "style of aggregate visualization for VTK output");
    validParamList->set<bool>                    ("aggregation: output file: graph edges",           false,                 "Whether to draw all node connections along with the aggregates.");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, "Aggregates");         //< factory which created aggregates
    Input(fineLevel, "DofsPerNode");        //< CoalesceAndDropFactory (needed for DofsPerNode variable)
    Input(fineLevel, "UnAmalgamationInfo"); //< AmalgamationFactory (needed for UnAmalgamationInfo variable)

    const ParameterList & pL = GetParameterList();
    //Only pull in coordinates if the user explicitly requests direct VTK output
    if(pL.isParameter("aggregation: output filename") && pL.get<std::string>("aggregation: output filename").length())
    {
      Input(fineLevel, "Coordinates");
      if(pL.isParameter("aggregation: output file: graph edges") && pL.get<bool>("aggregation: output file: graph edges"))
        Input(fineLevel, "A");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
    using namespace std;
    //Decide which build function to follow, based on input params
    const ParameterList& pL = GetParameterList();
    FactoryMonitor m(*this, "AggregationExportFactory", coarseLevel);
    std::string masterFilename = pL.get<std::string>("aggregation: output filename"); //filename parameter from master list
    std::string localFilename = pL.get<std::string>("Output filename");
    std::string filenameToWrite;
    bool useVTK = false;
    cout << "Master filename is \"" << masterFilename << "\"." << endl;
    if(masterFilename.length())
    {
      useVTK = true;
      filenameToWrite = masterFilename;
    }
    else
      filenameToWrite = localFilename;
    if(useVTK)
      cout << "Using VTK." << endl;
    else
      cout << "Not using vtk." << endl;
    Teuchos::RCP<Aggregates> aggregates      = Get< Teuchos::RCP<Aggregates> >(fineLevel,"Aggregates");
    LocalOrdinal DofsPerNode                 = Get< LocalOrdinal >            (fineLevel,"DofsPerNode");
    Teuchos::RCP<AmalgamationInfo> amalgInfo = Get< RCP<AmalgamationInfo> >   (fineLevel,"UnAmalgamationInfo");
    Teuchos::RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > coords = Teuchos::null;
    if(useVTK) //otherwise leave null, will not be accessed by non-vtk code
      coords = Get<RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> > >(fineLevel, "Coordinates");
    GetOStream(Runtime0) << "AggregationExportFactory: DofsPerNode: " << DofsPerNode << std::endl;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = aggregates->GetMap()->getComm();
    int numProcs = comm->getSize();
    int myRank   = comm->getRank();

    Teuchos::RCP<LocalOrdinalVector> vertex2AggId_vector = aggregates->GetVertex2AggId();
    Teuchos::RCP<LocalOrdinalVector> procWinner_vector   = aggregates->GetProcWinner();
    Teuchos::ArrayRCP<LocalOrdinal>  vertex2AggId        = aggregates->GetVertex2AggId()->getDataNonConst(0);
    Teuchos::ArrayRCP<LocalOrdinal>  procWinner          = aggregates->GetProcWinner()->getDataNonConst(0);

    // prepare for calculating global aggregate ids
    std::vector<GlobalOrdinal> numAggsGlobal (numProcs, 0);
    std::vector<GlobalOrdinal> numAggsLocal  (numProcs, 0);
    std::vector<GlobalOrdinal> minGlobalAggId(numProcs, 0);

    numAggsLocal[myRank] = aggregates->GetNumAggregates();
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, numProcs, &numAggsLocal[0], &numAggsGlobal[0]);
    for (int i = 1; i < Teuchos::as<int>(numAggsGlobal.size()); ++i) {
      numAggsGlobal [i] += numAggsGlobal[i-1];
      minGlobalAggId[i]  = numAggsGlobal[i-1];
    }

    ArrayRCP<LO>            aggStart;
    ArrayRCP<GlobalOrdinal> aggToRowMap;
    amalgInfo->UnamalgamateAggregates(*aggregates, aggStart, aggToRowMap);

    int timeStep = pL.get< int > ("Output file: time step");
    int iter = pL.get< int > ("Output file: iter");
    filenameToWrite = replaceAll(filenameToWrite, "%LEVELID",  toString(fineLevel.GetLevelID()));
    filenameToWrite = replaceAll(filenameToWrite, "%PROCID",   toString(myRank));
    filenameToWrite = replaceAll(filenameToWrite, "%TIMESTEP", toString(timeStep));
    filenameToWrite = replaceAll(filenameToWrite, "%ITER",     toString(iter));

    GetOStream(Runtime0) << "AggregationExportFactory: outputfile \"" << filenameToWrite << "\"" << std::endl;
    //does the user want a widely compatible .vtp file (xml formatted data) to visualize aggregates in ParaView?
    //If filename ends in .vtp (which it will have to be for an unstructured 'PolyData' VTK file), do that
    //This is the filename that VTK users will set, so check that for a valid filename

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
      //Only the main process will output anything
      if(myRank != 0)
        return;
      //Note: For now, this will only work with real scalars.
      using namespace std;
      if(sizeof(Scalar) != sizeof(double))
        throw runtime_error("Complex scalars not supported in VTK aggregate visualization.");
      //Make sure we have coordinates
      if(coords.is_null())
        throw runtime_error("AggExportFactory could not get coordinates, but they are required for VTK output.");
      int numNodes = coords->getGlobalLength();
      int dims = coords->getNumVectors();  //2D or 3D?
      //get access to the coord data
      Teuchos::ArrayRCP<const double> xCoords = coords->getData(0);
      Teuchos::ArrayRCP<const double> yCoords = coords->getData(1);
      Teuchos::ArrayRCP<const double> zCoords;
      const Teuchos::RCP<Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> > procWinnersVec = aggregates->GetProcWinner();
      Teuchos::ArrayRCP<const LocalOrdinal> procWinners = procWinnersVec->getData(0);
      if(dims == 3)
        zCoords = coords->getData(2);
      //Get the sizes of the aggregates to speed up grabbing node IDs
      Teuchos::ArrayRCP<LocalOrdinal> aggSizes = aggregates->ComputeAggregateSizes();
      string aggStyle = "Point Cloud";
      try
      {
        aggStyle = pL.get<string>("aggregation: output file: agg style"); //Let "Point Cloud" be the default style
      }
      catch(exception& e) {}
      fout << "<!--Aggregates Visualization-->" << endl;
      string indent = "";
      int numAggsInt = int(numAggs);
      bool drawEdges = pL.get<bool>("aggregation: output file: graph edges");
      Teuchos::RCP<Matrix> Amat = drawEdges ? Get<RCP<Matrix> >(fineLevel, "A") : Teuchos::null;
      std::vector<int> graphEdges;
      if(drawEdges)
      {
        //Add all (global) edges to graphEdges
        graphEdges.reserve(2 * Amat->getGlobalNumEntries()); //two nodes for each edge, num edges = nonzeros in A
        //todo: fill graphEdges with global data
      }
      if(aggStyle == "Point Cloud")
      {
        doPointCloud_(fout, xCoords, yCoords, zCoords, numNodes, numAggsInt, aggSizes, dims, vertex2AggId, procWinner, aggregates, drawEdges, graphEdges);
      }
      else if(aggStyle == "Jacks")
      {
        doJacks_(fout, xCoords, yCoords, zCoords, numNodes, numAggsInt, aggSizes, dims, vertex2AggId, procWinner, aggregates, drawEdges, graphEdges);
      }
      else if(aggStyle == "Jacks++")
      {
        doJacksPlus_(fout, xCoords, yCoords, zCoords, numNodes, numAggsInt, aggSizes, dims, vertex2AggId, procWinner, aggregates, drawEdges, graphEdges);
      }
      else if(aggStyle == "Convex Hulls")
      {
        doConvexHulls_(fout, xCoords, yCoords, zCoords, numNodes, numAggsInt, aggSizes, dims, vertex2AggId, procWinner, aggregates, drawEdges, graphEdges);
      }
      else if(aggStyle == "Alpha Hulls")
      {
        doAlphaHulls_(fout, xCoords, yCoords, zCoords, numNodes, numAggsInt, aggSizes, dims, vertex2AggId, procWinner, aggregates, drawEdges, graphEdges);
      }
      else
      {
        std::cout << "Warning: Unrecognized agg style.\nPossible values are Point Cloud, Jacks, Jacks++, Convex Hulls and Alpha Hulls.\nDefaulting to Point Cloud." << std::endl;
        doPointCloud_(fout, xCoords, yCoords, zCoords, numNodes, numAggsInt, aggSizes, dims, vertex2AggId, procWinner, aggregates, drawEdges, graphEdges);
      }
    }
    fout.close();
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doPointCloud_(std::ofstream& fout, Teuchos::ArrayRCP<const double>& xCoords, Teuchos::ArrayRCP<const double>& yCoords, Teuchos::ArrayRCP<const double>& zCoords, int numNodes, int numAggs, Teuchos::ArrayRCP<LocalOrdinal>& aggSizes, int dims, Teuchos::ArrayRCP<LocalOrdinal>& vertex2AggId, Teuchos::ArrayRCP<LocalOrdinal>& procWinners, Teuchos::RCP<Aggregates>& aggregates, bool doGraphEdges, std::vector<int>& graphConnectionsA)
  {
    using namespace std;
    fout << "<!--Point Cloud-->" << endl;
    fout << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << endl;
    fout << "  <UnstructuredGrid>" << endl;
    //Number of points in each "piece" will be the number of nodes in each aggregate
    size_t numCells = (size_t) numNodes;
    fout << "    <Piece NumberOfPoints=\"" << numNodes << "\" NumberOfCells=\"" << numCells << "\">" << endl;
    fout << "      <PointData Scalars=\"Node Aggregate Processor\">" << endl;
    string indent = "          ";
    fout << "        <DataArray type=\"Int32\" Name=\"Node\" format=\"ascii\">" << endl;
    fout << indent;
    for(int node = 0; node < numNodes; node++)
    {
      fout << node << " ";
      if(node % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Aggregate\" format=\"ascii\">" << endl;
    fout << indent;
    for(int node = 0; node < numNodes; node++)
    {
      fout << vertex2AggId[node] << " ";
      if(node % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << endl;
    fout << indent;
    for(int node = 0; node < numNodes; node++)
    {
      fout << procWinners[node] << " ";
      if(node % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </PointData>" << endl;
    //Write the point coordinates
    fout << "      <Points>" << endl;
    indent = "          ";
    fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    fout << indent;
    for(int node = 0; node < numNodes; node++)
    {
      fout << xCoords[node] << " " << yCoords[node] << " ";
      if(dims == 3)
        fout << zCoords[node] << " ";
      else
        fout << "0 ";
      if(node % 4 == 3) //Put 4 sets of coordinates on each line
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Points>" << endl;
    fout << "      <Cells>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    indent = "          ";
    fout << indent;
    for(int node = 0; node < numNodes; node++)
    {
      fout << node << " ";
      if(node % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    fout << indent;
    for(int node = 0; node < numNodes; node++)
    {
      fout << node + 1 << " "; //offset is where the cell ends, not begins, so put a 1 first
      if(node % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
    fout << indent;
    for(int node = 0; node < numNodes; node++)
    {
      fout << "1 ";
      if(node % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Cells>" << endl;
    fout << "    </Piece>" << endl << endl;
    fout << "  </UnstructuredGrid>" << std::endl;
    fout << "</VTKFile>" << std::endl;
  }
      
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doJacks_(std::ofstream& fout, Teuchos::ArrayRCP<const double>& xCoords, Teuchos::ArrayRCP<const double>& yCoords, Teuchos::ArrayRCP<const double>& zCoords, int numNodes, int numAggs, Teuchos::ArrayRCP<LocalOrdinal>& aggSizes, int dims, Teuchos::ArrayRCP<LocalOrdinal>& vertex2AggId, Teuchos::ArrayRCP<LocalOrdinal>& procWinners, Teuchos::RCP<Aggregates>& aggregates, bool doGraphEdges, std::vector<int>& graphConnections)
  {
    using namespace std;
    fout << "<!--Jacks-->" << endl;
    fout << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << endl;
    int totalEdges = numNodes - numAggs; //assuming all nodes are in exactly one aggregate
    fout << "  <UnstructuredGrid>" << endl;
    fout << "    <Piece NumberOfPoints=\"" << numNodes << "\" NumberOfCells=\"" << totalEdges << "\">" << endl;
    fout << "      <PointData Scalars=\"Node Aggregate Processor\">" << endl;
    string indent = "          ";
    fout << "        <DataArray type=\"Int32\" Name=\"Node\" format=\"ascii\">" << endl;
    fout << indent;
    for(int node = 0; node < numNodes; node++)
    {
      fout << node << " ";
      if(node % 8 == 7)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Aggregate\" format=\"ascii\">" << endl;
    fout << indent;
    for(int node = 0; node < numNodes; node++)
    {
      fout << vertex2AggId[node] << " ";
      if(node % 8 == 7)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << endl;
    fout << indent;
    for(int node = 0; node < numNodes; node++)
    {
      fout << procWinners[node] << " ";
      if(node % 8 == 7)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </PointData>" << endl;
    fout << "      <Points>" << endl;
    fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    fout << indent;
    for(int node = 0; node < numNodes; node++)
    {
      fout << xCoords[node] << " " << yCoords[node] << " ";
      if(dims == 2)
        fout << "0 ";
      else
        fout << zCoords[node] << " ";
      if(node % 4 == 3)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Points>" << endl;
    //Form list of node pairs
    vector<int> connections;
    //For each aggregate, find the root node then connect it with the other nodes in the aggregate
    //Doesn't matter the order, as long as all edges are found.
    int root = 0;
    for(int i = 0; i < numAggs; i++)
    {
      while(!aggregates->IsRoot(root))
        root++;
      int numInAggFound = 0;
      for(int j = 0; j < numNodes; j++)
      {
        if(j == root) //don't make a connection from the root to itself
        {
          numInAggFound++;
          continue;
        }
        if(vertex2AggId[root] == vertex2AggId[j])
        {
          connections.push_back(root);
          connections.push_back(j);
          numInAggFound++;
          if(numInAggFound == aggSizes[vertex2AggId[root]]) //don't spend more time looking if done with that root
            break;
        }
      }
      root++; //get set up to look for the next root
    }
    fout << "      <Cells>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    indent = "          ";
    fout << indent;
    for(int i = 0; i < int(connections.size()); i++)
    {
      fout << connections[i] << " ";
      if(i % 8 == 7)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(connections.size() / 2); i++)
    {
      //i goes from [0 to n), where n is the number of lines/cells
      fout << i * 2 + 2 << " ";
      if(i % 16 == 14)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
    indent = "          ";
    fout << indent;
    for(int i = 0; i < totalEdges; i++)
    {
      fout << "3 "; //VTK geometry type code for lines
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Cells>" <<  endl;
    fout << "    </Piece>" << endl;
    fout << "  </UnstructuredGrid>" << endl;
    fout << "</VTKFile>" << endl;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doJacksPlus_(std::ofstream& fout, Teuchos::ArrayRCP<const double>& xCoords, Teuchos::ArrayRCP<const double>& yCoords, Teuchos::ArrayRCP<const double>& zCoords, int numNodes, int numAggs, Teuchos::ArrayRCP<LocalOrdinal>& aggSizes, int dims, Teuchos::ArrayRCP<LocalOrdinal>& vertex2AggId, Teuchos::ArrayRCP<LocalOrdinal>& procWinners, Teuchos::RCP<Aggregates>& aggregates, bool doGraphEdges, std::vector<int>& graphConnections)
  {
    //TODO
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doConvexHulls_(std::ofstream& fout, Teuchos::ArrayRCP<const double>& xCoords, Teuchos::ArrayRCP<const double>& yCoords, Teuchos::ArrayRCP<const double>& zCoords, int numNodes, int numAggs, Teuchos::ArrayRCP<LocalOrdinal>& aggSizes, int dims, Teuchos::ArrayRCP<LocalOrdinal>& vertex2AggId, Teuchos::ArrayRCP<LocalOrdinal>& procWinners, Teuchos::RCP<Aggregates>& aggregates, bool doGraphEdges, std::vector<int>& graphConnections)
  {
    if(dims == 2)
      doConvexHulls2D_(fout, xCoords, yCoords, numNodes, numAggs, aggSizes, vertex2AggId, procWinners, aggregates, doGraphEdges, graphConnections);
    else
      doConvexHulls3D_(fout, xCoords, yCoords, zCoords, numNodes, numAggs, aggSizes, vertex2AggId, procWinners, aggregates, doGraphEdges, graphConnections);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doConvexHulls2D_(std::ofstream& fout, Teuchos::ArrayRCP<const double>& xCoords, Teuchos::ArrayRCP<const double>& yCoords, int numNodes, int numAggs, Teuchos::ArrayRCP<LocalOrdinal>& aggSizes, Teuchos::ArrayRCP<LocalOrdinal>& vertex2AggId, Teuchos::ArrayRCP<LocalOrdinal>& procWinners, Teuchos::RCP<Aggregates>& aggregates, bool doGraphEdges, std::vector<int>& graphConnections)
  {
    using namespace std;
    vector<vector<int> > hulls; //outer vector contains aggregates, inner vector contains unique node IDs
    for(int agg = 0; agg < numAggs; agg++)
    {
      list<int> aggNodes;
      for(int i = 0; i < numNodes; i++)
      {
        if(vertex2AggId[i] == agg)
          aggNodes.push_back(i);
      }
      //have a list of nodes in the aggregate
      vector<int> thisHull;
      if(aggNodes.size() == 0)
      {
        hulls.push_back(thisHull);
        continue;
      }
      if(aggNodes.size() == 1)
      {
        thisHull.push_back(aggNodes.front());
        hulls.push_back(thisHull);
        continue;
      }
      if(aggNodes.size() == 2)
      {
        thisHull.push_back(aggNodes.front());
        thisHull.push_back(aggNodes.back());
        hulls.push_back(thisHull);
        continue;
      }
      //check if all points are collinear, need to explicitly draw a line in that case.
      bool collinear = true; //assume true at first, if a segment not parallel to others then clear
      {
        list<int>::iterator it = aggNodes.begin();
        vec3_ firstPoint(xCoords[*it], yCoords[*it], 0);
        it++;
        vec3_ secondPoint(xCoords[*it], yCoords[*it], 0);
        it++;  //it now points to third node in the aggregate
        vec3_ norm1(-(secondPoint.y - firstPoint.y), secondPoint.x - firstPoint.x, 0);
        do
        {
          vec3_ thisNorm(yCoords[*it] - firstPoint.y, firstPoint.x - xCoords[*it], 0);
          //rotate one of the vectors by 90 degrees so that dot product is 0 if the two are parallel
          double temp = thisNorm.x;
          thisNorm.x = thisNorm.y;
          thisNorm.y = temp;
          double comp = dotProduct_(norm1, thisNorm);
          if(-1e-8 > comp || comp > 1e-8)
          {
            collinear = false;
            break;
          }
          it++;
        }
        while(it != aggNodes.end());
      }
      if(collinear)
      {
        //find the most distant two points in the plane and use as endpoints of line representing agg
        list<int>::iterator min = aggNodes.begin();    //min X then min Y where x is a tie 
        list<int>::iterator max = aggNodes.begin(); //max X then max Y where x is a tie
        for(list<int>::iterator it = ++aggNodes.begin(); it != aggNodes.end(); it++)
        {
          if(xCoords[*it] < xCoords[*min])
            min = it;
          else if(xCoords[*it] == xCoords[*min])
          {
            if(yCoords[*it] < yCoords[*min])
              min = it;
          }
          if(xCoords[*it] > xCoords[*max])
            max = it;
          else if(xCoords[*it] == xCoords[*max])
          {
            if(yCoords[*it] > yCoords[*max])
              max = it;
          }
        }
        //Just set up a line between nodes *min and *max
        thisHull.push_back(*min);
        thisHull.push_back(*max);
        hulls.push_back(thisHull);
        continue; //jump to next aggregate in loop
      }
      list<int>::iterator min = aggNodes.begin();
      for(list<int>::iterator it = ++aggNodes.begin(); it != aggNodes.end(); it++)
      {
        if(xCoords[*it] < xCoords[*min])
          min = it;
        else if(xCoords[*it] == xCoords[*min])
        {
          if(yCoords[*it] < yCoords[*min])
            min = it;
        }
      }
      //this is the most common case: at least 3 nodes making up a polygon with positive area
      //do Jarvis march on these points to compute convex hull
      //start with "min" point
      thisHull.push_back(*min);
      //repeatedly sweep through aggNodes (with a pivot at thisHull.back()) counterclockwise. Connect up with last node found this way ("left-most"). If there is a tie for the left-most node, take the most distant one from thisHull.back().
      bool includeMin = false;
      while(1)
      {
        list<int>::iterator leftMost = aggNodes.begin();
        if(!includeMin && leftMost == min)
        {
          leftMost++;
        }
        list<int>::iterator it = leftMost;
        it++;
        while(it != aggNodes.end())
        {
          if(it == min && !includeMin) //don't compare to min on very first sweep
          {
            it++;
            continue;
          }
          //see if it is in front of line containing nodes thisHull.back() and leftMost
          //first get the left normal of leftMost - thisHull.back() (<dy, -dx>)
          vec3_ testNorm(yCoords[*leftMost] - yCoords[thisHull.back()], -(xCoords[*leftMost] - xCoords[thisHull.back()]), 0);
          //now dot testNorm with *it - leftMost. If dot is positive, leftMost becomes it. If dot is zero, take one further from thisHull.back().
          vec3_ itVec(xCoords[*it] - xCoords[*leftMost], yCoords[*it] - yCoords[*leftMost], 0);
          double dotProd = dotProduct_(testNorm, itVec);
          if(-1e-10 < dotProd && dotProd < 1e-10)
          {
            //thisHull.back(), it and leftMost are collinear.
            //Just sum the differences in x and differences in y for each and compare, don't need pythag
            double itDist = fabs(xCoords[*it] - xCoords[thisHull.back()]) + fabs(yCoords[*it] - yCoords[thisHull.back()]);
            double leftMostDist = fabs(xCoords[*leftMost] - xCoords[thisHull.back()]) + fabs(yCoords[*leftMost] - yCoords[thisHull.back()]);
            if(itDist > leftMostDist)
              leftMost = it;
          }
          else if(dotProd > 0)
            leftMost = it;
          it++;
        }
        //if leftMost is min, then the loop is complete.
        if(leftMost == min)
        {
          //for(int i = 0; i < thisHull.size(); i++)
            //cout << thisHull[i] << " ";
          //cout << endl << endl;
          hulls.push_back(thisHull);
          break; //this goes to the next aggregate
        }
        //add leftMost to thisHull, and remove it from aggNodes so later searches are faster
        thisHull.push_back(*leftMost);
        aggNodes.erase(leftMost);
        includeMin = true; //have found second point (the one after min) so now include min in the searches
      }
    }
    //Collapse set of points only to those that are used to define convex hull geometry
    vector<int> uniquePoints;
    for(int i = 0; i < int(hulls.size()); i++)
    {
      for(int j = 0; j < int(hulls[i].size()); j++)
      {
        uniquePoints.push_back(hulls[i][j]);
      }
    }
    sort(uniquePoints.begin(), uniquePoints.end());
    unique(uniquePoints.begin(), uniquePoints.end());
    //now go through each hull's set of nodes and replace global node values with indices in uniquePoints
    for(int i = 0; i < int(hulls.size()); i++)
    {
      for(int j = 0; j < int(hulls[i].size()); j++)
      {
        hulls[i][j] = distance(uniquePoints.begin(), find(uniquePoints.begin(), uniquePoints.end(), hulls[i][j]));
      }
    }
    //now have all the nodes that are hull vertices and the sizes of the hulls
    //write the PolyData and piece headers b/c now know how many elements will be in the piece
    fout << "<!--2D Convex Hulls-->" << endl;
    fout << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << endl;
    fout << "  <UnstructuredGrid>" << endl;
    fout << "    <Piece NumberOfPoints=\"" << uniquePoints.size() << "\" NumberOfCells=\"" << numAggs << "\">" << endl;
    string indent = "          ";
    fout << "      <PointData Scalars=\"Node Aggregate Processor\">" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Node\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(uniquePoints.size()); i++)
    {
      fout << uniquePoints[i] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Aggregate\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(uniquePoints.size()); i++)
    {
      fout << vertex2AggId[uniquePoints[i]] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(uniquePoints.size()); i++)
    {
      fout << procWinners[uniquePoints[i]] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </PointData>" << endl;
    //write out point coordinates (for each aggregate, draw each node)
    //hullPoints has those nodes in that order
    indent = "          ";
    fout << "      <Points>" << endl;
    fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(uniquePoints.size()); i++)
    {
      fout << xCoords[uniquePoints[i]] << " " << yCoords[uniquePoints[i]] << " 0 ";
      if(i % 4 == 3)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Points>" << endl;
    fout << "      <Cells>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    fout << indent;
    int totalItems = 0;
    for(int i = 0; i < int(hulls.size()); i++)
    {
      for(int j = 0; j < int(hulls[i].size()); j++)
      {
        fout << hulls[i][j] << " ";
        totalItems++;
        if(totalItems % 10 == 9)
          fout << endl << indent;
      }
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    //hullSizes has the size (in points) of each hull, so I need to add them together to get cumulative offset
    int offset = 0;
    fout << indent;
    for(int i = 0; i < int(hulls.size()); i++)
    {
      offset += hulls[i].size();
      fout << offset << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(hulls.size()); i++)
    {
      //Draw points and lines if polygon would be degenerate and invisible
      int thisHullSize = int(hulls[i].size());
      if(thisHullSize == 0)
        continue;
      else if(thisHullSize == 1)
        fout << "1 ";
      else if(thisHullSize == 2)
        fout << "3 ";
      else
        fout << "7 ";
      if(i % 30 == 29)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Cells>" << endl;
    fout << "    </Piece>" << endl;
    fout << "  </UnstructuredGrid>" << endl;
    fout << "</VTKFile>" << endl;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doConvexHulls3D_(std::ofstream& fout, Teuchos::ArrayRCP<const double>& xCoords, Teuchos::ArrayRCP<const double>& yCoords, Teuchos::ArrayRCP<const double>& zCoords, int numNodes, int numAggs, Teuchos::ArrayRCP<LocalOrdinal>& aggSizes, Teuchos::ArrayRCP<LocalOrdinal>& vertex2AggId, Teuchos::ArrayRCP<LocalOrdinal>& procWinners, Teuchos::RCP<Aggregates>& aggregates, bool doGraphEdges, std::vector<int>& graphConnections)
  {
    using namespace std;
    //Use 3D quickhull algo.
    //Vector of node indices representing triangle vertices
    //Note: Calculate the hulls first since will only include point data for points in the hulls
    vector<int> vertIndices;   //node ids (in triplets) that form triangles in hulls
    vector<int> geomSizes;    //list of sizes of geometry elements: 1 = point, 2 = line, 3 = triangle, >3 = polygon
    //Effectively the size() of vertIndices after each hull is added to it
    typedef list<int>::iterator Iter;
    for(int agg = 0; agg < numAggs; agg++)
    {
      list<int> aggNodes; //At first, list of all nodes in the aggregate. As nodes are enclosed or included by/in hull, remove them
      for(int i = 0; i < numNodes; i++)
      {
        if(vertex2AggId[i] == agg)
          aggNodes.push_back(i);
      }
      //First, check anomalous cases
      if(aggNodes.size() == 0)
        throw runtime_error("An aggregate has zero nodes in it!");
      else if(aggNodes.size() == 1)
      {
        vertIndices.push_back(aggNodes.front());
        geomSizes.push_back(1);
        continue;
      }
      else if(aggNodes.size() == 2)
      {
        vertIndices.push_back(aggNodes.front());
        vertIndices.push_back(aggNodes.back());
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
        vec3_ collinearVec(colCheckDX, colCheckDY, colCheckDZ);
        //normalize collinearVec
        double collinearVecMag = magnitude_(collinearVec);
        collinearVec.x /= collinearVecMag;
        collinearVec.y /= collinearVecMag;
        collinearVec.z /= collinearVecMag;
        vec3_ firstPoint(xCoords[aggNodes.front()], yCoords[aggNodes.front()], zCoords[aggNodes.front()]);
        for(Iter it = colCheck; it != aggNodes.end(); it++)
        {
          vec3_ thisVec(xCoords[*it], yCoords[*it], zCoords[*it]);
          vec3_ vecDiff = vecSubtract_(thisVec, firstPoint);
          //normalize vecDiff so that it can be directly compared to collinearVec
          double vecDiffMag = magnitude_(vecDiff);
          vecDiff.x /= vecDiffMag;
          vecDiff.y /= vecDiffMag;
          vecDiff.z /= vecDiffMag;
          //compare x, y and z separately and give some slack for rounding error
          vec3_ compare = vecSubtract_(vecDiff, collinearVec);
          if(compare.x < -1e-10 || compare.x > 1e-10 || compare.y < -1e-10 || compare.y > 1e-10 || compare.z < -1e-10 || compare.z > 1e-10)
          {
            areCollinear = false;
            break;
          }
        }
      }
      if(areCollinear)
      {
        //find the endpoints of segment describing all the points
        //compare x, if tie compare y, if tie compare z
        Iter min = aggNodes.begin();
        Iter max = aggNodes.begin();
        Iter it = ++aggNodes.begin();
        for(; it != aggNodes.end(); it++)
        {
          if(xCoords[*it] < xCoords[*min])
            min = it;
          else if(xCoords[*it] == xCoords[*min])
          {
            if(yCoords[*it] < yCoords[*min])
              min = it;
            else if(yCoords[*it] == yCoords[*min])
            {
              if(zCoords[*it] < zCoords[*min])
                min = it;
            }
          }
          if(xCoords[*it] > xCoords[*max])
            max = it;
          else if(xCoords[*it] == xCoords[*max])
          {
            if(yCoords[*it] > yCoords[*max])
              max = it;
            else if(yCoords[*it] == yCoords[*max])
            {
              if(zCoords[*it] > zCoords[*max])
                max = it;
            }
          }
        }
        vertIndices.push_back(*min);
        vertIndices.push_back(*max);
        geomSizes.push_back(2);
        continue;
      }
      Iter exIt = aggNodes.begin(); //iterator to be used for searching for min/max x/y/z
      int extremeSix[] = {*exIt, *exIt, *exIt, *exIt, *exIt, *exIt}; //nodes with minimumX, maxX, minY ...
      exIt++;
      for(; exIt != aggNodes.end(); exIt++)
      {
        Iter it = exIt;
        if(xCoords[*it] < xCoords[extremeSix[0]] || (xCoords[*it] == xCoords[extremeSix[0]] && yCoords[*it] < yCoords[extremeSix[0]]) || (xCoords[*it] == xCoords[extremeSix[0]] && yCoords[*it] == yCoords[extremeSix[0]] && zCoords[*it] < zCoords[extremeSix[0]]))
          extremeSix[0] = *it;
        if(xCoords[*it] > xCoords[extremeSix[1]] || (xCoords[*it] == xCoords[extremeSix[1]] && yCoords[*it] > yCoords[extremeSix[1]]) || (xCoords[*it] == xCoords[extremeSix[1]] && yCoords[*it] == yCoords[extremeSix[1]] && zCoords[*it] > zCoords[extremeSix[1]]))
          extremeSix[1] = *it;
        if(yCoords[*it] < yCoords[extremeSix[2]] || (yCoords[*it] == yCoords[extremeSix[2]] && zCoords[*it] < zCoords[extremeSix[2]]) || (yCoords[*it] == yCoords[extremeSix[2]] && zCoords[*it] == zCoords[extremeSix[2]] && xCoords[*it] < xCoords[extremeSix[2]]))
          extremeSix[2] = *it;
        if(yCoords[*it] > yCoords[extremeSix[3]] || (yCoords[*it] == yCoords[extremeSix[3]] && zCoords[*it] > zCoords[extremeSix[3]]) || (yCoords[*it] == yCoords[extremeSix[3]] && zCoords[*it] == zCoords[extremeSix[3]] && xCoords[*it] > xCoords[extremeSix[3]]))
          extremeSix[3] = *it;
        if(zCoords[*it] < zCoords[extremeSix[4]] || (zCoords[*it] == zCoords[extremeSix[4]] && xCoords[*it] < xCoords[extremeSix[4]]) || (zCoords[*it] == zCoords[extremeSix[4]] && xCoords[*it] == xCoords[extremeSix[4]] && yCoords[*it] < yCoords[extremeSix[4]]))
          extremeSix[4] = *it;
        if(zCoords[*it] > zCoords[extremeSix[5]] || (zCoords[*it] == zCoords[extremeSix[5]] && xCoords[*it] > xCoords[extremeSix[5]]) || (zCoords[*it] == zCoords[extremeSix[5]] && xCoords[*it] == xCoords[extremeSix[5]] && yCoords[*it] > zCoords[extremeSix[5]]))
          extremeSix[5] = *it;
      }
      vec3_ extremeVectors[6];
      for(int i = 0; i < 6; i++)
      {
        vec3_ thisExtremeVec(xCoords[extremeSix[i]], yCoords[extremeSix[i]], zCoords[extremeSix[i]]);
        extremeVectors[i] = thisExtremeVec;
      }
      double maxDist = 0;
      int base1 = 0; //ints from 0-5: which pair out of the 6 extreme points are the most distant? (indices in extremeSix and extremeVectors)
      int base2 = 0;
      for(int i = 0; i < 5; i++)
      {
        for(int j = i + 1; j < 6; j++)
        {
          double thisDist = distance_(extremeVectors[i], extremeVectors[j]);
          if(thisDist > maxDist)
          {
            maxDist = thisDist;
            base1 = i;
            base2 = j;
          }
        }
      }
      list<Triangle_> hullBuilding;    //each Triangle is a triplet of nodes (int IDs) that form a triangle
      //remove base1 and base2 iters from aggNodes, they are known to be in the hull
      aggNodes.remove(extremeSix[base1]);
      aggNodes.remove(extremeSix[base2]);
      //extremeSix[base1] and [base2] still have the vec3_ representation
      Triangle_ tri;
      tri.v1 = extremeSix[base1];
      tri.v2 = extremeSix[base2];
      //Now find the point that is furthest away from the line between base1 and base2
      maxDist = 0;
      //need the vectors to do "quadruple product" formula
      vec3_ b1 = extremeVectors[base1];
      vec3_ b2 = extremeVectors[base2];
      Iter thirdNode;
      for(Iter node = aggNodes.begin(); node != aggNodes.end(); node++)
      {
        vec3_ nodePos(xCoords[*node], yCoords[*node], zCoords[*node]);
        double dist = magnitude_(crossProduct_(vecSubtract_(nodePos, b1), vecSubtract_(nodePos, b2))) / magnitude_(vecSubtract_(b2, b1));
        if(dist > maxDist)
        {
          maxDist = dist;
          thirdNode = node;
        }
      }
      //Now know the last node in the first triangle
      tri.v3 = *thirdNode;
      hullBuilding.push_back(tri);
      vec3_ b3(xCoords[*thirdNode], yCoords[*thirdNode], zCoords[*thirdNode]);
      aggNodes.erase(thirdNode);
      //Find the fourth node (most distant from triangle) to form tetrahedron
      maxDist = 0;
      int fourthVertex = -1;
      for(Iter node = aggNodes.begin(); node != aggNodes.end(); node++)
      {
        vec3_ thisNode(xCoords[*node], yCoords[*node], zCoords[*node]);
        double nodeDist = pointDistFromTri_(thisNode, b1, b2, b3);
        if(nodeDist > maxDist)
        {
          maxDist = nodeDist;
          fourthVertex = *node;
        }
      }
      aggNodes.remove(fourthVertex);
      vec3_ b4(xCoords[fourthVertex], yCoords[fourthVertex], zCoords[fourthVertex]);
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
      vec3_ barycenter((b1.x + b2.x + b3.x + b4.x) / 4.0, (b1.y + b2.y + b3.y + b4.y) / 4.0, (b1.z + b2.z + b3.z + b4.z) / 4.0);
      for(list<Triangle_>::iterator tetTri = hullBuilding.begin(); tetTri != hullBuilding.end(); tetTri++)
      {
        vec3_ triVert1(xCoords[tetTri->v1], yCoords[tetTri->v1], zCoords[tetTri->v1]);
        vec3_ triVert2(xCoords[tetTri->v2], yCoords[tetTri->v2], zCoords[tetTri->v2]);
        vec3_ triVert3(xCoords[tetTri->v3], yCoords[tetTri->v3], zCoords[tetTri->v3]);
        vec3_ trinorm = getNorm_(triVert1, triVert2, triVert3);
        vec3_ ptInPlane = tetTri == hullBuilding.begin() ? b1 : b4; //first triangle definitely has b1 in it, other three definitely have b4
        if(isInFront_(barycenter, ptInPlane, trinorm))
        {
          //don't want the faces of the tetrahedron to face inwards (towards barycenter) so reverse orientation
          //by swapping two vertices
          int temp = tetTri->v1;
          tetTri->v1 = tetTri->v2;
          tetTri->v2 = temp;
        }
      }
//cout << "Done making the initial tetrahedron's faces face outward." << endl;
      //now, have starting polyhedron in hullBuilding (all faces are facing outwards according to getNorm_) and remaining nodes to process are in aggNodes
      //recursively, for each triangle, make a list of the points that are 'in front' of the triangle. Find the point with the maximum distance from the triangle.
      //Add three new triangles, each sharing one edge with the original triangle but now with the most distant point as a vertex. Remove the most distant point from
      //the list of all points that need to be processed. Also from that list remove all points that are in front of the original triangle but not in front of all three
      //new triangles, since they are now enclosed in the hull.
      //Construct point lists for each face of the tetrahedron individually.
      vec3_ trinorms[4]; //normals to the four tetrahedron faces, now oriented outwards
      int index = 0;
      for(list<Triangle_>::iterator tetTri = hullBuilding.begin(); tetTri != hullBuilding.end(); tetTri++)
      {
        vec3_ triVert1(xCoords[tetTri->v1], yCoords[tetTri->v1], zCoords[tetTri->v1]);
        vec3_ triVert2(xCoords[tetTri->v2], yCoords[tetTri->v2], zCoords[tetTri->v2]);
        vec3_ triVert3(xCoords[tetTri->v3], yCoords[tetTri->v3], zCoords[tetTri->v3]);
        trinorms[index] = getNorm_(triVert1, triVert2, triVert3);
        index++;
      }
      list<int> startPoints1;
      list<int> startPoints2;
      list<int> startPoints3;
      list<int> startPoints4;
      //scope this so that 'point' is not in function scope
      {
        Iter point = aggNodes.begin();
        while(!aggNodes.empty())  //this removes points one at a time as they are put in startPointsN or are already done
        {
          vec3_ pointVec(xCoords[*point], yCoords[*point], zCoords[*point]);
          //Note: Because of the way the tetrahedron faces are constructed above,
          //face 1 definitely contains b1 and faces 2-4 definitely contain b4.
          if(isInFront_(pointVec, b1, trinorms[0]))
          {
            startPoints1.push_back(*point);
            point = aggNodes.erase(point);
          }
          else if(isInFront_(pointVec, b4, trinorms[1]))
          {
            startPoints2.push_back(*point);
            point = aggNodes.erase(point);
          }
          else if(isInFront_(pointVec, b4, trinorms[2]))
          {
            startPoints3.push_back(*point);
            point = aggNodes.erase(point);
          }
          else if(isInFront_(pointVec, b4, trinorms[3]))
          {
            startPoints4.push_back(*point);
            point = aggNodes.erase(point);
          }
          else
          {
            point = aggNodes.erase(point); //points here are already inside tetrahedron.
          }
        }
        //Call processTriangle_ for each triangle in the initial tetrahedron, one at a time.
      }
      typedef list<Triangle_>::iterator TriIter;
      TriIter firstTri = hullBuilding.begin();
      Triangle_ start1 = *firstTri;
      firstTri++;
      Triangle_ start2 = *firstTri;
      firstTri++;
      Triangle_ start3 = *firstTri;
      firstTri++;
      Triangle_ start4 = *firstTri;
      //kick off depth-first recursive filling of hullBuilding list with all triangles in the convex hull
      if(!startPoints1.empty())
        processTriangle_(hullBuilding, start1, startPoints1, barycenter, xCoords, yCoords, zCoords);
      if(!startPoints2.empty())
        processTriangle_(hullBuilding, start2, startPoints2, barycenter, xCoords, yCoords, zCoords);
      if(!startPoints3.empty())
        processTriangle_(hullBuilding, start3, startPoints3, barycenter, xCoords, yCoords, zCoords);
      if(!startPoints4.empty())
        processTriangle_(hullBuilding, start4, startPoints4, barycenter, xCoords, yCoords, zCoords);
      //hullBuilding now has all triangles that make up this hull.
      //Dump hullBuilding info into the list of all triangles for the scene.
      vertIndices.reserve(vertIndices.size() + 3 * hullBuilding.size());
      for(TriIter hullTri = hullBuilding.begin(); hullTri != hullBuilding.end(); hullTri++)
      {
        vertIndices.push_back(hullTri->v1);
        vertIndices.push_back(hullTri->v2);
        vertIndices.push_back(hullTri->v3);
        geomSizes.push_back(3);
      }
    }
    //Make a set (sorted list w/no duplicates) of the node IDs that actually appear in convex hulls
    //Interior points won't be in VTK file
    vector<int> includePoints = vertIndices;
    sort(includePoints.begin(), includePoints.end());
    vector<int>::iterator uniqueEnd = unique(includePoints.begin(), includePoints.end());
    includePoints.erase(uniqueEnd, includePoints.end());
    int uniquePoints = int(includePoints.size());
    //Replace actual node values in vertIndices by the equivalent unique point index
    for(int i = 0; i < int(vertIndices.size()); i++)
    {
      //binary search for vertIndices[i] in includePoints, and get that index in includePoints
      int loInd = 0;
      int hiInd = includePoints.size() - 1;
      while(hiInd >= loInd)
      {
        int midInd = (loInd + hiInd) / 2;
        if(includePoints[midInd] == vertIndices[i])
        {
          vertIndices[i] = midInd;
          break;
        }
        else if(includePoints[midInd] < vertIndices[i]) //use upper half of previous range
          loInd = midInd + 1;
        else                                          //use lower half
          hiInd = midInd - 1;
      }
    }
    fout << "<!--3D Convex Hulls-->" << endl;
    fout << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << endl;
    fout << "  <UnstructuredGrid>" << endl;
    //Number of points in each "piece" will be the number of nodes in each aggregate
    fout << "    <Piece NumberOfPoints=\"" << uniquePoints << "\" NumberOfCells=\"" << geomSizes.size() << "\">" << endl;
    fout << "      <PointData Scalars=\"Node Aggregate Processor\">" << endl;
    string indent = "          ";
    fout << "        <DataArray type=\"Int32\" Name=\"Node\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < uniquePoints; i++)
    {
      fout << includePoints[i] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Aggregate\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < uniquePoints; i++)
    {
      fout << vertex2AggId[includePoints[i]] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < uniquePoints; i++)
    {
      fout << procWinners[includePoints[i]] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </PointData>" << endl;
    fout << "      <Points>" << endl;
    fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    indent = "          ";
    fout << indent;
    for(int i = 0; i < uniquePoints; i++)
    {
      fout << xCoords[includePoints[i]] << " " << yCoords[includePoints[i]] << " " << zCoords[includePoints[i]] << " ";
      if(i % 4 == 3)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Points>" << endl;
    fout << "      <Cells>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    fout << indent;
    for(int i = 0; i < int(vertIndices.size()); i++)
    {
      fout << vertIndices[i] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    fout << indent;
    int totalOffset = 0;
    for(int i = 0; i < int(geomSizes.size()); i++)
    {
      totalOffset += geomSizes[i];
      fout << totalOffset << " ";
      if(i % 10 == 9)
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
          fout << "1 ";
          break;
        case 2:
          fout << "3 ";
          break;
        case 3:
          fout << "5 ";
          break;
        default:
          fout << "7 ";
      }
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Cells>" << endl;
    fout << "    </Piece>" << endl;
    fout << "  </UnstructuredGrid>" << endl;
    fout << "</VTKFile>" << endl;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doAlphaHulls_(std::ofstream& fout, Teuchos::ArrayRCP<const double>& xCoords, Teuchos::ArrayRCP<const double>& yCoords, Teuchos::ArrayRCP<const double>& zCoords, int numNodes, int numAggs, Teuchos::ArrayRCP<LocalOrdinal>& aggSizes, int dims, Teuchos::ArrayRCP<LocalOrdinal>& vertex2AggId, Teuchos::ArrayRCP<LocalOrdinal>& procWinner, Teuchos::RCP<Aggregates>& aggregates, bool doGraphEdges, std::vector<int>& graphConnections)
  {
    using namespace std;
    //TODO
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceAll(std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) const {
    while(1) {
      const int pos = result.find(replaceWhat);
      if (pos == -1)
        break;
      result.replace(pos, replaceWhat.size(), replaceWithWhat);
    }
    return result;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  vec3_ AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::crossProduct_(vec3_ v1, vec3_ v2)
  {
    return vec3_(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dotProduct_(vec3_ v1, vec3_ v2)
  {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isInFront_(vec3_ point, vec3_ inPlane, vec3_ n)
  {
    vec3_ rel(point.x - inPlane.x, point.y - inPlane.y, point.z - inPlane.z); //position of the point relative to the plane
    return dotProduct_(rel, n) > 1e-12 ? true : false;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::magnitude_(vec3_ vec)
  {
    return sqrt(dotProduct_(vec, vec));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::distance_(vec3_ p1, vec3_ p2)
  {
    return magnitude_(vecSubtract_(p1, p2));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  vec3_ AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::vecSubtract_(vec3_ v1, vec3_ v2)
  {
    return vec3_(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);    
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  vec3_ AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNorm_(vec3_ v1, vec3_ v2, vec3_ v3) //normal to face of triangle (will be outward rel. to polyhedron) (v1, v2, v3 are in CCW order when normal is toward viewpoint)
  {
    return crossProduct_(vecSubtract_(v2, v1), vecSubtract_(v3, v1));
  }

  //get minimum distance from 'point' to plane containing v1, v2, v3 (or the triangle with v1, v2, v3 as vertices)
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::pointDistFromTri_(vec3_ point, vec3_ v1, vec3_ v2, vec3_ v3)
  {
    using namespace std;
    vec3_ norm = getNorm_(v1, v2, v3);
    //must normalize the normal vector
    double normScl = magnitude_(norm);
    norm.x /= normScl;
    norm.y /= normScl;
    norm.z /= normScl;
    //cout << "Have a normal vector <" << norm.x << "," << norm.y << "," << norm.z << ">" << endl;
    //cout << "The distance from the point (" << point.x << "," << point.y << "," << point.z << ")" << endl;
    //cout << "to the triangle containing points (" << v1.x << "," << v1.y << "," << v1.z << "), (" << v2.x << "," << v2.y << "," << v2.z << "), and (" << v3.x << "," << v3.y << "," << v3.z << ") is ";
    double rv = fabs(dotProduct_(norm, vecSubtract_(point, v1)));
    //cout << rv << endl;
    return rv;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::vector<Triangle_> AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::processTriangle_(std::list<Triangle_>& tris, Triangle_ tri, std::list<int>& pointsInFront, vec3_& barycenter, ArrayRCP<const double>& xCoords, ArrayRCP<const double>& yCoords, ArrayRCP<const double>& zCoords)
  {
    //*tri is in the tris list, and is the triangle to process here. tris is a complete list of all triangles in the hull so far. pointsInFront is only a list of the nodes in front of tri. Need coords also.
    //precondition: each triangle is already oriented so that getNorm_(v1, v2, v3) points outward (away from interior of hull)
    //First find the point furthest from triangle.
    using namespace std;
    typedef std::list<int>::iterator Iter;
    typedef std::list<Triangle_>::iterator TriIter;
    typedef list<pair<int, int> >::iterator EdgeIter;
//    cout << "Processing triangle with nodes: " << tri->v1 << ", " << tri->v2 << ", " << tri->v3 << "." << endl;
//    cout << "Nodes in front of this triangle: ";
//    for(Iter it = pointsInFront.begin(); it != pointsInFront.end() ;it++)
 //     cout << *it << " ";
 //   cout << endl;
    double maxDist = 0;
    //Need vector representations of triangle's vertices
    vec3_ v1(xCoords[tri.v1], yCoords[tri.v1], zCoords[tri.v1]);
    vec3_ v2(xCoords[tri.v2], yCoords[tri.v2], zCoords[tri.v2]);
    vec3_ v3(xCoords[tri.v3], yCoords[tri.v3], zCoords[tri.v3]);
    vec3_ farPointVec; //useful to have both the point's coordinates and it's position in the list
    Iter farPoint = pointsInFront.begin();
    for(Iter point = pointsInFront.begin(); point != pointsInFront.end(); point++)
    {
      vec3_ pointVec(xCoords[*point], yCoords[*point], zCoords[*point]);
      double dist = pointDistFromTri_(pointVec, v1, v2, v3);
      if(dist > maxDist)
      {
        dist = maxDist;
        farPointVec = pointVec;
        farPoint = point;
      }
    }
    //Find all the triangles that the point is in front of (can be more than 1)
    //At the same time, remove them from tris, as every one will be replaced later
    vector<Triangle_> visible; //use a list of iterators so that the underlying object is still in tris
    for(TriIter it = tris.begin(); it != tris.end();)
    {
      vec3_ vec1(xCoords[it->v1], yCoords[it->v1], zCoords[it->v1]);
      vec3_ vec2(xCoords[it->v2], yCoords[it->v2], zCoords[it->v2]);
      vec3_ vec3(xCoords[it->v3], yCoords[it->v3], zCoords[it->v3]);
      vec3_ norm = getNorm_(vec1, vec2, vec3);
      if(isInFront_(farPointVec, vec1, norm))
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
    for(vector<Triangle_>::iterator it = visible.begin(); it != visible.end(); it++)
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
    list<Triangle_> newTris;
    for(EdgeIter it = horizon.begin(); it != horizon.end(); it++)
    {
      Triangle_ t(it->first, it->second, *farPoint);
      newTris.push_back(t);
    }
    //Ensure every new triangle is oriented outwards, using the barycenter of the initial tetrahedron
    vector<Triangle_> trisToProcess;
    vector<list<int> > newFrontPoints;
    for(TriIter it = newTris.begin(); it != newTris.end(); it++)
    {
      vec3_ t1(xCoords[it->v1], yCoords[it->v1], zCoords[it->v1]);
      vec3_ t2(xCoords[it->v2], yCoords[it->v2], zCoords[it->v2]);
      vec3_ t3(xCoords[it->v3], yCoords[it->v3], zCoords[it->v3]);
      if(isInFront_(barycenter, t1, getNorm_(t1, t2, t3)))
      {
        //need to swap two vertices to flip orientation of triangle
        int temp = it->v1;
        vec3_ tempVec = t1;
        it->v1 = it->v2;
        t1 = t2;
        it->v2 = temp;
        t2 = tempVec;
      }
      vec3_ outwardNorm = getNorm_(t1, t2, t3); //now definitely points outwards
      //Add the triangle to tris
      tris.push_back(*it);
      trisToProcess.push_back(tris.back());
      //Make a list of the points that are in front of nextToProcess, to be passed in for processing
      list<int> newInFront;
      for(Iter point = pointsInFront.begin(); point != pointsInFront.end();)
      {
        vec3_ pointVec(xCoords[*point], yCoords[*point], zCoords[*point]);
        if(isInFront_(pointVec, t1, outwardNorm))
        {
          newInFront.push_back(*point);
          point = pointsInFront.erase(point);
        }
        else
          point++;
      }
      newFrontPoints.push_back(newInFront);
    }
    vector<Triangle_> allRemoved; //list of all invalid iterators that were erased by calls to processTriangle_ below
    for(int i = 0; i < int(trisToProcess.size()); i++)
    {
      if(!newFrontPoints[i].empty())
      {
        //Comparing the 'triangle to process' to the one for this call prevents infinite recursion/stack overflow.
        //TODO: Why was it doing that? Rounding error? Make more robust fix. But this does work for the time being.
        if(find(allRemoved.begin(), allRemoved.end(), trisToProcess[i]) == allRemoved.end() && !(trisToProcess[i] == tri))
        {
          //cout << "Going to process the nodes in front of triangle with vertices " << trisToProcess[i].v1 << ", " << trisToProcess[i].v2 << ", and " << trisToProcess[i].v3 << endl;
          vector<Triangle_> removedList = processTriangle_(tris, trisToProcess[i], newFrontPoints[i], barycenter, xCoords, yCoords, zCoords);
          for(int j = 0; j < int(removedList.size()); j++)
            allRemoved.push_back(removedList[j]);
        }
      }
    }
    return visible;
  }
} // namespace MueLu

#endif /* MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_ */
