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
#include "MueLu_Utilities.hpp"
#include <vector>
#include <list>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <cstdio>
#include <cmath>
//For alpha hulls (is optional feature requiring a third-party library)
#ifdef HAVE_MUELU_CGAL //Include all headers needed for both 2D and 3D fixed-alpha alpha shapes
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Delaunay_triangulation_2.h"
#include "CGAL/Delaunay_triangulation_3.h"
#include "CGAL/Alpha_shape_2.h"
#include "CGAL/Fixed_alpha_shape_3.h"
#include "CGAL/algorithm.h"
#endif

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    std::string output_msg = "Output filename template (%TIMESTEP is replaced by \'Output file: time step\' variable,"
        "%ITER is replaced by \'Output file: iter\' variable, %LEVELID is replaced level id, %PROCID is replaced by processor id)";
    std::string output_def = "aggs_level%LEVELID_proc%PROCID.out";

    validParamList->set< RCP<const FactoryBase> >("A", Teuchos::null, "Factory for A.");
    validParamList->set< RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory for Coordinates.");
    validParamList->set< RCP<const FactoryBase> >("Graph", Teuchos::null, "Factory for Graph.");
    validParamList->set< RCP<const FactoryBase> >("Aggregates", Teuchos::null, "Factory for Aggregates.");
    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Factory for UnAmalgamationInfo.");
    validParamList->set< RCP<const FactoryBase> >("DofsPerNode", Teuchos::null, "Factory for DofsPerNode.");
    // CMS/BMK: Old style factory-only options.  Deprecate me.
    validParamList->set< std::string >           ("Output filename",           output_def, output_msg);
    validParamList->set< int >                   ("Output file: time step",             0, "time step variable for output file name");
    validParamList->set< int >                   ("Output file: iter",                  0, "nonlinear iteration variable for output file name");

    // New-style master list options (here are same defaults as in masterList.xml)
    validParamList->set< std::string >           ("aggregation: output filename",                    "",                    "filename for VTK-style visualization output");
    validParamList->set< int >                   ("aggregation: output file: time step",             0,                     "time step variable for output file name");// Remove me?
    validParamList->set< int >                   ("aggregation: output file: iter",                  0,                     "nonlinear iteration variable for output file name");//Remove me?
    validParamList->set<std::string>             ("aggregation: output file: agg style",             "Point Cloud",         "style of aggregate visualization for VTK output");
    validParamList->set<bool>                    ("aggregation: output file: fine graph edges",      false,                 "Whether to draw all fine node connections along with the aggregates.");
    validParamList->set<bool>                    ("aggregation: output file: coarse graph edges",    false,                 "Whether to draw all coarse node connections along with the aggregates.");
    validParamList->set<bool>                    ("aggregation: output file: build colormap",        false,                 "Whether to output a random colormap for ParaView in a separate XML file.");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    Input(fineLevel, "Aggregates");         //< factory which created aggregates
    Input(fineLevel, "DofsPerNode");        //< CoalesceAndDropFactory (needed for DofsPerNode variable)
    Input(fineLevel, "UnAmalgamationInfo"); //< AmalgamationFactory (needed for UnAmalgamationInfo variable)

    const ParameterList & pL = GetParameterList();
    //Only pull in coordinates if the user explicitly requests direct VTK output, so as not to break uses of old code
    if(pL.isParameter("aggregation: output filename") && pL.get<std::string>("aggregation: output filename").length())
    {
      Input(fineLevel, "Coordinates");
      Input(fineLevel, "A");
      Input(fineLevel, "Graph");
      if(pL.get<bool>("aggregation: output file: coarse graph edges"))
      {
        Input(coarseLevel, "Coordinates");
        Input(coarseLevel, "A");
        Input(coarseLevel, "Graph");
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
    using namespace std;
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
    if(numProcs == 0)
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
      string indent = "";
      nodeMap_ = Amat->getMap();
      for(LocalOrdinal i = 0; i < numNodes_; i++)
      {
        isRoot_.push_back(aggregates->IsRoot(i));
      }
      //If problem is serial and not outputting fine nor coarse graph edges, don't make pvtu file
      //Otherwise create it
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
      else if(aggStyle == "Alpha Hulls")
      {
        #ifdef HAVE_MUELU_CGAL
        doAlphaHulls_(vertices, geomSizes);
        #else
        std::cout << "   Trilinos was not configured with CGAL so Alpha Hulls not available.\n   Using Convex Hulls instead." << std::endl;
        doConvexHulls_(vertices, geomSizes);
        #endif
      }
      else
      {
        std::cout << "   Warning: Unrecognized agg style.\n   Possible values are Point Cloud, Jacks, Jacks++, Convex Hulls and Alpha Hulls.\n   Defaulting to Point Cloud." << std::endl;
        aggStyle = "Point Cloud";
        doPointCloud_(vertices, geomSizes);
      }
      writeFile_(fout, aggStyle, vertices, geomSizes);
      if(doCoarseGraphEdges_)
      {
        string fname = filenameToWrite;
        string cEdgeFile = fname.insert(fname.rfind(".vtu"), "-coarsegraph");
        std::ofstream edgeStream(cEdgeFile.c_str());
        doGraphEdges_(edgeStream, Ac, coarseGraph, false, DofsPerNode);
      }
      if(doFineGraphEdges_)
      {
        string fname = filenameToWrite;
        string fEdgeFile = fname.insert(fname.rfind(".vtu"), "-finegraph");
        std::ofstream edgeStream(fEdgeFile.c_str());
        doGraphEdges_(edgeStream, Amat, fineGraph, true, DofsPerNode);
      }
      if(myRank == 0 && pL.get<bool>("aggregation: output file: build colormap"))
      {
        buildColormap_();
      }
    }
    fout.close();
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doPointCloud_(std::vector<int>& vertices, std::vector<int>& geomSize) const
  {
    vertices.reserve(vertices.size() + numNodes_);
    geomSize.reserve(geomSize.size() + numNodes_);
    for(int i = 0; i < numNodes_; i++)
    {
      vertices.push_back(i);
      geomSize.push_back(1);
    }
  }
      
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doJacks_(std::vector<int>& vertices, std::vector<int>& geomSizes) const
  {
    //For each aggregate, find the root node then connect it with the other nodes in the aggregate
    //Doesn't matter the order, as long as all edges are found.
    vertices.reserve(vertices.size() + 3 * (numNodes_ - numAggs_));
    geomSizes.reserve(vertices.size() + 2 * (numNodes_ - numAggs_));
    int root = 0;
    for(int i = 0; i < numAggs_; i++) //TODO: Replace this O(n^2) with a better way?
    {
      while(!isRoot_[root])
        root++;
      int aggSize = aggSizes_[vertex2AggId_[root]];
      if(aggSize == 1)
      {
        vertices.push_back(root);
        geomSizes.push_back(1);
      }
      else
      {
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
      }
      root++; //get set up to look for the next root
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doJacksPlus_(std::vector<int>& vertices, std::vector<int>& geomSizes) const
  {
    //TODO
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doConvexHulls_(std::vector<int>& vertices, std::vector<int>& geomSizes) const
  {
    if(dims_ == 2)
      doConvexHulls2D_(vertices, geomSizes);
    else
      doConvexHulls3D_(vertices, geomSizes);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doConvexHulls2D_(std::vector<int>& vertices, std::vector<int>& geomSizes) const
  {
    using namespace std;
    for(int agg = 0; agg < numAggs_; agg++)
    {
      list<int> aggNodes;
      for(int i = 0; i < numNodes_; i++)
      {
        if(vertex2AggId_[i] == agg)
          aggNodes.push_back(i);
      }
      //have a list of nodes in the aggregate
      if(aggNodes.size() == 0)
      {
        throw runtime_error("Aggregate contains zero nodes!");
      }
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
        geomSizes.push_back(2);
        continue;
      }
      //check if all points are collinear, need to explicitly draw a line in that case.
      bool collinear = true; //assume true at first, if a segment not parallel to others then clear
      {
        list<int>::iterator it = aggNodes.begin();
        vec3_ firstPoint(xCoords_[*it], yCoords_[*it], 0);
        it++;
        vec3_ secondPoint(xCoords_[*it], yCoords_[*it], 0);
        it++;  //it now points to third node in the aggregate
        vec3_ norm1(-(secondPoint.y - firstPoint.y), secondPoint.x - firstPoint.x, 0);
        do
        {
          vec3_ thisNorm(yCoords_[*it] - firstPoint.y, firstPoint.x - xCoords_[*it], 0);
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
          if(xCoords_[*it] < xCoords_[*min])
            min = it;
          else if(xCoords_[*it] == xCoords_[*min])
          {
            if(yCoords_[*it] < yCoords_[*min])
              min = it;
          }
          if(xCoords_[*it] > xCoords_[*max])
            max = it;
          else if(xCoords_[*it] == xCoords_[*max])
          {
            if(yCoords_[*it] > yCoords_[*max])
              max = it;
          }
        }
        //Just set up a line between nodes *min and *max
        vertices.push_back(*min);
        vertices.push_back(*max);
        geomSizes.push_back(2);
        continue; //jump to next aggregate in loop
      }
      vector<vec2_> points;
      vector<int> nodes;
      for(list<int>::iterator it = aggNodes.begin(); it != aggNodes.end(); it++)
      {
        points.push_back(vec2_(xCoords_[*it], yCoords_[*it]));
        nodes.push_back(*it);
      }
      vector<int> hull = giftWrap_(points, nodes);
      vertices.reserve(vertices.size() + hull.size());
      for(size_t i = 0; i < hull.size(); i++)
      {
        vertices.push_back(hull[i]);
      }
      geomSizes.push_back(hull.size());
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doConvexHulls3D_(std::vector<int>& vertices, std::vector<int>& geomSizes) const
  {
    using namespace std;
    //Use 3D quickhull algo.
    //Vector of node indices representing triangle vertices
    //Note: Calculate the hulls first since will only include point data for points in the hulls
    //Effectively the size() of vertIndices after each hull is added to it
    typedef list<int>::iterator Iter;
    for(int agg = 0; agg < numAggs_; agg++)
    {
      list<int> aggNodes; //At first, list of all nodes in the aggregate. As nodes are enclosed or included by/in hull, remove them
      for(int i = 0; i < numNodes_; i++)
      {
        if(vertex2AggId_[i] == agg)
          aggNodes.push_back(i);
      }
      //First, check anomalous cases
      if(aggNodes.size() == 0)
      {
        throw runtime_error("An aggregate has zero nodes in it!");
      }
      else if(aggNodes.size() == 1)
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
        Iter it = aggNodes.begin();
        vec3_ firstVec(xCoords_[*it], yCoords_[*it], zCoords_[*it]);
        vec3_ comp;
        {
          it++;
          vec3_ secondVec(xCoords_[*it], yCoords_[*it], zCoords_[*it]); //cross this with other vectors to compare
          comp = vecSubtract_(secondVec, firstVec);
          it++;
        }
        for(; it != aggNodes.end(); it++)
        {
          vec3_ thisVec(xCoords_[*it], yCoords_[*it], zCoords_[*it]);
          vec3_ cross = crossProduct_(vecSubtract_(thisVec, firstVec), comp);
          if(magnitude_(cross) > 1e-10)
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
          if(xCoords_[*it] < xCoords_[*min])
            min = it;
          else if(xCoords_[*it] == xCoords_[*min])
          {
            if(yCoords_[*it] < yCoords_[*min])
              min = it;
            else if(yCoords_[*it] == yCoords_[*min])
            {
              if(zCoords_[*it] < zCoords_[*min])
                min = it;
            }
          }
          if(xCoords_[*it] > xCoords_[*max])
            max = it;
          else if(xCoords_[*it] == xCoords_[*max])
          {
            if(yCoords_[*it] > yCoords_[*max])
              max = it;
            else if(yCoords_[*it] == yCoords_[*max])
            {
              if(zCoords_[*it] > zCoords_[*max])
                max = it;
            }
          }
        }
        vertices.push_back(*min);
        vertices.push_back(*max);
        geomSizes.push_back(2);
        continue;
      }
      bool areCoplanar = true;
      {
        //number of points is known to be >= 3
        Iter vert = aggNodes.begin();
        vec3_ v1(xCoords_[*vert], yCoords_[*vert], zCoords_[*vert]);
        vert++;
        vec3_ v2(xCoords_[*vert], yCoords_[*vert], zCoords_[*vert]);
        vert++;
        vec3_ v3(xCoords_[*vert], yCoords_[*vert], zCoords_[*vert]);
        vert++;
        //Make sure the first three points aren't also collinear (need a non-degenerate triangle to get a normal)
        while(magnitude_(crossProduct_(vecSubtract_(v1, v2), vecSubtract_(v1, v3))) < 1e-10)
        {
          //Replace the third point with the next point
          v3 = vec3_(xCoords_[*vert], yCoords_[*vert], zCoords_[*vert]);
          vert++;
        }
        for(; vert != aggNodes.end(); vert++)
        {
          vec3_ pt(xCoords_[*vert], yCoords_[*vert], zCoords_[*vert]);
          if(fabs(pointDistFromTri_(pt, v1, v2, v3)) > 1e-12)
          {
            areCoplanar = false;
            break;
          }
        }
        if(areCoplanar)
        {
          //do 2D convex hull
          vec3_ planeNorm = getNorm_(v1, v2, v3);
          planeNorm.x = fabs(planeNorm.x);
          planeNorm.y = fabs(planeNorm.y);
          planeNorm.z = fabs(planeNorm.z);
          vector<vec2_> points;
          vector<int> nodes;
          if(planeNorm.x >= planeNorm.y && planeNorm.x >= planeNorm.z)
          {
            //project points to yz plane to make hull
            for(Iter it = aggNodes.begin(); it != aggNodes.end(); it++)
            {
              nodes.push_back(*it);
              points.push_back(vec2_(yCoords_[*it], zCoords_[*it]));
            }
          }
          if(planeNorm.y >= planeNorm.x && planeNorm.y >= planeNorm.z)
          {
            //use xz
            for(Iter it = aggNodes.begin(); it != aggNodes.end(); it++)
            {
              nodes.push_back(*it);
              points.push_back(vec2_(xCoords_[*it], zCoords_[*it]));
            }           
          }
          if(planeNorm.z >= planeNorm.x && planeNorm.z >= planeNorm.y)
          {
            for(Iter it = aggNodes.begin(); it != aggNodes.end(); it++)
            {
              nodes.push_back(*it);
              points.push_back(vec2_(xCoords_[*it], yCoords_[*it]));
            }           
          }
          vector<int> convhull2d = giftWrap_(points, nodes);
          geomSizes.push_back(convhull2d.size());
          vertices.reserve(vertices.size() + convhull2d.size());
          for(size_t i = 0; i < convhull2d.size(); i++)
            vertices.push_back(convhull2d[i]);
          continue;
        }
      }
      Iter exIt = aggNodes.begin(); //iterator to be used for searching for min/max x/y/z
      int extremeSix[] = {*exIt, *exIt, *exIt, *exIt, *exIt, *exIt}; //nodes with minimumX, maxX, minY ...
      exIt++;
      for(; exIt != aggNodes.end(); exIt++)
      {
        Iter it = exIt;
        if(xCoords_[*it] < xCoords_[extremeSix[0]] ||
          (xCoords_[*it] == xCoords_[extremeSix[0]] && yCoords_[*it] < yCoords_[extremeSix[0]]) ||
          (xCoords_[*it] == xCoords_[extremeSix[0]] && yCoords_[*it] == yCoords_[extremeSix[0]] && zCoords_[*it] < zCoords_[extremeSix[0]]))
            extremeSix[0] = *it;
        if(xCoords_[*it] > xCoords_[extremeSix[1]] ||
          (xCoords_[*it] == xCoords_[extremeSix[1]] && yCoords_[*it] > yCoords_[extremeSix[1]]) ||
          (xCoords_[*it] == xCoords_[extremeSix[1]] && yCoords_[*it] == yCoords_[extremeSix[1]] && zCoords_[*it] > zCoords_[extremeSix[1]]))
            extremeSix[1] = *it;
        if(yCoords_[*it] < yCoords_[extremeSix[2]] ||
          (yCoords_[*it] == yCoords_[extremeSix[2]] && zCoords_[*it] < zCoords_[extremeSix[2]]) ||
          (yCoords_[*it] == yCoords_[extremeSix[2]] && zCoords_[*it] == zCoords_[extremeSix[2]] && xCoords_[*it] < xCoords_[extremeSix[2]]))
            extremeSix[2] = *it;
        if(yCoords_[*it] > yCoords_[extremeSix[3]] ||
          (yCoords_[*it] == yCoords_[extremeSix[3]] && zCoords_[*it] > zCoords_[extremeSix[3]]) ||
          (yCoords_[*it] == yCoords_[extremeSix[3]] && zCoords_[*it] == zCoords_[extremeSix[3]] && xCoords_[*it] > xCoords_[extremeSix[3]]))
            extremeSix[3] = *it;
        if(zCoords_[*it] < zCoords_[extremeSix[4]] ||
          (zCoords_[*it] == zCoords_[extremeSix[4]] && xCoords_[*it] < xCoords_[extremeSix[4]]) ||
          (zCoords_[*it] == zCoords_[extremeSix[4]] && xCoords_[*it] == xCoords_[extremeSix[4]] && yCoords_[*it] < yCoords_[extremeSix[4]]))
            extremeSix[4] = *it;
        if(zCoords_[*it] > zCoords_[extremeSix[5]] ||
          (zCoords_[*it] == zCoords_[extremeSix[5]] && xCoords_[*it] > xCoords_[extremeSix[5]]) ||
          (zCoords_[*it] == zCoords_[extremeSix[5]] && xCoords_[*it] == xCoords_[extremeSix[5]] && yCoords_[*it] > zCoords_[extremeSix[5]]))
            extremeSix[5] = *it;
      }
      vec3_ extremeVectors[6];
      for(int i = 0; i < 6; i++)
      {
        vec3_ thisExtremeVec(xCoords_[extremeSix[i]], yCoords_[extremeSix[i]], zCoords_[extremeSix[i]]);
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
        vec3_ nodePos(xCoords_[*node], yCoords_[*node], zCoords_[*node]);
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
      vec3_ b3(xCoords_[*thirdNode], yCoords_[*thirdNode], zCoords_[*thirdNode]);
      aggNodes.erase(thirdNode);
      //Find the fourth node (most distant from triangle) to form tetrahedron
      maxDist = 0;
      int fourthVertex = -1;
      for(Iter node = aggNodes.begin(); node != aggNodes.end(); node++)
      {
        vec3_ thisNode(xCoords_[*node], yCoords_[*node], zCoords_[*node]);
        double nodeDist = pointDistFromTri_(thisNode, b1, b2, b3);
        if(nodeDist > maxDist)
        {
          maxDist = nodeDist;
          fourthVertex = *node;
        }
      }
      aggNodes.remove(fourthVertex);
      vec3_ b4(xCoords_[fourthVertex], yCoords_[fourthVertex], zCoords_[fourthVertex]);
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
        vec3_ triVert1(xCoords_[tetTri->v1], yCoords_[tetTri->v1], zCoords_[tetTri->v1]);
        vec3_ triVert2(xCoords_[tetTri->v2], yCoords_[tetTri->v2], zCoords_[tetTri->v2]);
        vec3_ triVert3(xCoords_[tetTri->v3], yCoords_[tetTri->v3], zCoords_[tetTri->v3]);
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
        vec3_ triVert1(xCoords_[tetTri->v1], yCoords_[tetTri->v1], zCoords_[tetTri->v1]);
        vec3_ triVert2(xCoords_[tetTri->v2], yCoords_[tetTri->v2], zCoords_[tetTri->v2]);
        vec3_ triVert3(xCoords_[tetTri->v3], yCoords_[tetTri->v3], zCoords_[tetTri->v3]);
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
          vec3_ pointVec(xCoords_[*point], yCoords_[*point], zCoords_[*point]);
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
        processTriangle_(hullBuilding, start1, startPoints1, barycenter);
      if(!startPoints2.empty())
        processTriangle_(hullBuilding, start2, startPoints2, barycenter);
      if(!startPoints3.empty())
        processTriangle_(hullBuilding, start3, startPoints3, barycenter);
      if(!startPoints4.empty())
        processTriangle_(hullBuilding, start4, startPoints4, barycenter);
      //hullBuilding now has all triangles that make up this hull.
      //Dump hullBuilding info into the list of all triangles for the scene.
      vertices.reserve(vertices.size() + 3 * hullBuilding.size());
      for(TriIter hullTri = hullBuilding.begin(); hullTri != hullBuilding.end(); hullTri++)
      {
        vertices.push_back(hullTri->v1);
        vertices.push_back(hullTri->v2);
        vertices.push_back(hullTri->v3);
        geomSizes.push_back(3);
      }
    }
  }

#ifdef HAVE_MUELU_CGAL
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doAlphaHulls_(std::vector<int>& vertices, std::vector<int>& geomSizes) const
  {
    using namespace std;
    if(dims_ == 2)
      doAlphaHulls2D_(vertices, geomSizes);
    else if(dims_ == 3)
      doAlphaHulls3D_(vertices, geomSizes);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doAlphaHulls2D_(std::vector<int>& vertices, std::vector<int>& geomSizes) const
  {
    const double ALPHA_VAL = 2; //Make configurable?
    using namespace std;
    //CGAL setup
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::FT FT;
    typedef K::Point_2 Point;
    typedef K::Segment_2 Segment;
    typedef CGAL::Alpha_shape_vertex_base_2<K> Vb;
    typedef CGAL::Alpha_shape_face_base_2<K> Fb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
    typedef CGAL::Delaunay_triangulation_2<K,Tds> Triangulation_2;
    typedef CGAL::Alpha_shape_2<Triangulation_2> Alpha_shape_2;
    typedef Alpha_shape_2::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;
    for(int i = 0; i < numAggs_; i++)
    {
      //Populate a list of Point_2 for this aggregate
      list<Point> aggPoints;
      vector<int> aggNodes;
      for(int j = 0; j < numNodes_; j++)
      {
        if(vertex2AggId_[j] == i)
        {
          Point p(xCoords_[j], yCoords_[j]);
          aggPoints.push_back(p);
          aggNodes.push_back(j);
        }
      }
      Alpha_shape_2 hull(aggPoints.begin(), aggPoints.end(), FT(ALPHA_VAL), Alpha_shape_2::GENERAL);
      vector<Segment> segments;
      alpha_edges(hull, back_inserter(segments));
      vertices.reserve(vertices.size() + 2 * segments.size());
      geomSizes.reserve(geomSizes.size() + segments.size());
      for(size_t j = 0; j < segments.size(); j++)
      {
        for(size_t k = 0; k < aggNodes.size(); k++)
        {
          if(fabs(segments[j][0].x == xCoords_[aggNodes[k]]) < 1e-12 && fabs(segments[j][0].y == yCoords_[aggNodes[k]]) < 1e-12)
          {
            vertices.push_back(aggNodes[k]);
            break;
          }
        }
        for(size_t k = 0; k < aggNodes.size(); k++)
        {
          if(fabs(segments[j][1].x  == xCoords_[aggNodes[k]]) < 1e-12 && fabs(segments[j][1].y == yCoords_[aggNodes[k]]) < 1e-12)
          {
            vertices.push_back(aggNodes[k]);
            break;
          }
        }
        geomSizes.push_back(2); //all cells are line segments
      }
    }
  }
  
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doAlphaHulls3D_(std::vector<int>& vertices, std::vector<int>& geomSizes) const
  {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel Gt;
    typedef CGAL::Alpha_shape_cell_base_3<Gt> Fb;
    typedef CGAL::Triangulation_data_structure_3<Vb,Fb> Tds;
    typedef CGAL::Delaunay_triangulation_3<Gt,Tds> Triangulation_3;
    typedef Gt::Point_3 Point;
    typedef Alpha_shape_3::Alpha_iterator Alpha_iterator;
    typedef Alpha_shape_3::Cell_handle Cell_handle;
    typedef Alpha_shape_3::Vertex_handle Vertex_handle;
    typedef Alpha_shape_3::Facet Facet;
    typedef Alpha_shape_3::Edge Edge;
    typedef Gt::Weighted_point Weighted_point;
    typedef Gt::Bare_point Bare_point;
    const double ALPHA_VAL = 2; //Make configurable?
    using namespace std;
    for(int i = 0; i < numAggs_; i++)
    {
      list<Point> aggPoints;
      vector<int> aggNodes;
      for(int j = 0; j < numNodes_; j++)
      {
        if(vertex2AggId[j] == i)
        {
          Point p(xCoords_[j], yCoords_[j], zCoords_[j]);
          aggPoints.push_back(p);
          aggNodes.push_back(j);
        }
      }
      Fixed_alpha_shape_3 hull(aggPoints.begin(), aggPoints.end(), FT(ALPHA_VAL));
      list<Cell_handle> cells;
      list<Facet> facets;
      list<Edge> edges;
      hull.get_alpha_shape_cells(back_inserter(cells));
      hull.get_alpha_shape_facets(back_inserter(facets));
      hull.get_alpha_shape_edges(back_inserter(edges));
      for(size_t j = 0; j < cells.size(); j++)
      {
        Point tetPoints[4];
        tetPoints[0] = cells[j]->vertex(0);
        tetPoints[1] = cells[j]->vertex(1);
        tetPoints[2] = cells[j]->vertex(2);
        tetPoints[3] = cells[j]->vertex(3);
        for(int k = 0; k < 4; k++)
        {
          for(size_t l = 0; l < aggNodes.size(); l++)
          {
            if(fabs(tetPoints[k].x - xCoords_[aggNodes[l]]) < 1e-12 &&
               fabs(tetPoints[k].y - yCoords_[aggNodes[l]]) < 1e-12 &&
               fabs(tetPoints[k].z - zCoords_[aggNodes[l]]) < 1e-12)
            {
              vertices.push_back(aggNodes[l]);
              break;
            }
          }
        }
        geomSizes.push_back(-10); //tetrahedron
      }
      for(size_t j = 0; j < facets.size(); j++)
      {
        int indices[3];
        indices[0] = (facets[i].second + 1) % 4;
        indices[1] = (facets[i].second + 2) % 4;
        indices[2] = (facets[i].second + 3) % 4;
        if(facets[i].second % 2 == 0)
          swap(indices[0], indices[1]);
        Point facetPts[3];
        facetPts[0] = facets[i].first->vertex(indices[0])->point();
        facetPts[1] = facets[i].first->vertex(indices[1])->point();
        facetPts[2] = facets[i].first->vertex(indices[2])->point();
        //add triangles in terms of node indices
        for(size_t l = 0; l < aggNodes.size(); l++)
        {
          if(fabs(facetPts[k].x - xCoords_[aggNodes[l]]) < 1e-12 &&
             fabs(facetPts[k].y - yCoords_[aggNodes[l]]) < 1e-12 &&
             fabs(facetPts[k].z - zCoords_[aggNodes[l]]) < 1e-12)
          {
            vertices.push_back(aggNodes[l]);
            break;
          }
        }
        geomSizes.push_back(3);
      }
      for(size_t j = 0; j < edges.size(); j++)
      {
        
      }
    }
  }
#endif
  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceAll(std::string result, const std::string& replaceWhat, const std::string& replaceWithWhat) const
  {
    while(1)
    {
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
  double AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dotProduct_(vec2_ v1, vec2_ v2)
  {
    return v1.x * v2.x + v1.y * v2.y;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::isInFront_(vec3_ point, vec3_ inPlane, vec3_ n)
  {
    vec3_ rel(point.x - inPlane.x, point.y - inPlane.y, point.z - inPlane.z); //position of the point relative to the plane
    return dotProduct_(rel, n) > 1e-12;
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
  vec2_ AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::vecSubtract_(vec2_ v1, vec2_ v2)
  {
    return vec2_(v1.x - v2.x, v1.y - v2.y);    
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

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  vec2_ AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNorm_(vec2_ v) //"normal" to a 2D vector - just rotate 90 degrees to left
  {
    return vec2_(v.y, -v.x);
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
    double rv = fabs(dotProduct_(norm, vecSubtract_(point, v1)));
    //cout << rv << endl;
    return rv;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::vector<Triangle_> AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::processTriangle_(std::list<Triangle_>& tris, Triangle_ tri, std::list<int>& pointsInFront, vec3_& barycenter) const
  {
    //*tri is in the tris list, and is the triangle to process here. tris is a complete list of all triangles in the hull so far. pointsInFront is only a list of the nodes in front of tri. Need coords also.
    //precondition: each triangle is already oriented so that getNorm_(v1, v2, v3) points outward (away from interior of hull)
    //First find the point furthest from triangle.
    using namespace std;
    typedef std::list<int>::iterator Iter;
    typedef std::list<Triangle_>::iterator TriIter;
    typedef list<pair<int, int> >::iterator EdgeIter;
    double maxDist = 0;
    //Need vector representations of triangle's vertices
    vec3_ v1(xCoords_[tri.v1], yCoords_[tri.v1], zCoords_[tri.v1]);
    vec3_ v2(xCoords_[tri.v2], yCoords_[tri.v2], zCoords_[tri.v2]);
    vec3_ v3(xCoords_[tri.v3], yCoords_[tri.v3], zCoords_[tri.v3]);
    vec3_ farPointVec; //useful to have both the point's coordinates and it's position in the list
    Iter farPoint = pointsInFront.begin();
    for(Iter point = pointsInFront.begin(); point != pointsInFront.end(); point++)
    {
      vec3_ pointVec(xCoords_[*point], yCoords_[*point], zCoords_[*point]);
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
      vec3_ vec1(xCoords_[it->v1], yCoords_[it->v1], zCoords_[it->v1]);
      vec3_ vec2(xCoords_[it->v2], yCoords_[it->v2], zCoords_[it->v2]);
      vec3_ vec3(xCoords_[it->v3], yCoords_[it->v3], zCoords_[it->v3]);
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
      vec3_ t1(xCoords_[it->v1], yCoords_[it->v1], zCoords_[it->v1]);
      vec3_ t2(xCoords_[it->v2], yCoords_[it->v2], zCoords_[it->v2]);
      vec3_ t3(xCoords_[it->v3], yCoords_[it->v3], zCoords_[it->v3]);
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
        vec3_ pointVec(xCoords_[*point], yCoords_[*point], zCoords_[*point]);
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
    for(size_t i = 0; i < trisToProcess.size(); i++)
    {
      if(!newFrontPoints[i].empty())
      {
        //Comparing the 'triangle to process' to the one for this call prevents infinite recursion/stack overflow.
        //TODO: Why was it doing that? Rounding error? Make more robust fix. But this does work for the time being.
        if(find(allRemoved.begin(), allRemoved.end(), trisToProcess[i]) == allRemoved.end() && !(trisToProcess[i] == tri))
        {
          vector<Triangle_> removedList = processTriangle_(tris, trisToProcess[i], newFrontPoints[i], barycenter);
          for(size_t j = 0; j < removedList.size(); j++)
            allRemoved.push_back(removedList[j]);
        }
      }
    }
    return visible;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::vector<int> AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::giftWrap_(std::vector<vec2_>& points, std::vector<int>& nodes) const
  {
    using namespace std;
    if(points.size() < 3)
      throw runtime_error("Gift wrap algorithm input has to have at least 3 points.");
    vector<int> hull;
    vector<int>::iterator minNode = nodes.begin();
    vec2_ min = points[0];
    for(vector<int>::iterator it = nodes.begin(); it != nodes.end(); it++)
    {
      int i = it - nodes.begin();
      if(points[i].x < min.x || (fabs(points[i].x - min.x) < 1e-12 && points[i].y < min.y))
      {
        min = points[i];
        minNode = it;
      }
    }
    hull.push_back(*minNode);
    bool includeMin = false;
    while(1)
    {
      vector<int>::iterator leftMost = nodes.begin();
      if(!includeMin && leftMost == minNode)
      {
        leftMost++;
      }
      vector<int>::iterator it = leftMost;
      it++;
      for(; it != nodes.end(); it++)
      {
        if(it == minNode && !includeMin) //don't compare to min on very first sweep
          continue;
        if(*it == hull.back())
          continue;
        //see if it is in front of line containing nodes thisHull.back() and leftMost
        //first get the left normal of leftMost - thisHull.back() (<dy, -dx>)
        vec2_ leftMostVec = points[leftMost - nodes.begin()];
        vec2_ lastVec(xCoords_[hull.back()], yCoords_[hull.back()]);
        vec2_ testNorm = getNorm_(vecSubtract_(leftMostVec, lastVec));
        //now dot testNorm with *it - leftMost. If dot is positive, leftMost becomes it. If dot is zero, take one further from thisHull.back().
        vec2_ itVec(xCoords_[*it], yCoords_[*it]);
        double dotProd = dotProduct_(testNorm, vecSubtract_(itVec, lastVec));
        if(-1e-10 < dotProd && dotProd < 1e-10)
        {
          //thisHull.back(), it and leftMost are collinear.
          //Just sum the differences in x and differences in y for each and compare to get further one, don't need distance formula
          vec2_ itDist = vecSubtract_(itVec, lastVec);
          vec2_ leftMostDist = vecSubtract_(leftMostVec, lastVec);
          if(fabs(itDist.x) + fabs(itDist.y) > fabs(leftMostDist.x) + fabs(leftMostDist.y))
            leftMost = it;
        }
        else if(dotProd > 0)
          leftMost = it;
      }
      //if leftMost is min, then the loop is complete.
      if(*leftMost == *minNode)
        break;
      hull.push_back(*leftMost);
      includeMin = true; //have found second point (the one after min) so now include min in the searches
    }
    return hull;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::doGraphEdges_(std::ofstream& fout, Teuchos::RCP<Matrix>& A, Teuchos::RCP<GraphBase>& G, bool fine, int dofs) const
  {
    using namespace std;
    ArrayView<const Scalar> values;
    ArrayView<const LocalOrdinal> neighbors;
    //Allow two different colors of connections (by setting "aggregates" scalar to CONTRAST_1 or CONTRAST_2)
    vector<pair<int, int> > vert1; //vertices (node indices)
    vector<pair<int, int> > vert2; //size of every cell is assumed to be 2 vertices, since all edges are drawn as lines
    if(A->isGloballyIndexed())
    {
      ArrayView<const GlobalOrdinal> indices;
      for(GlobalOrdinal globRow = 0; globRow < GlobalOrdinal(A->getGlobalNumRows()); globRow++)
      {
        if(dofs == 1)
          A->getGlobalRowView(globRow, indices, values);
        neighbors = G->getNeighborVertices((LocalOrdinal) globRow);
        int gEdge = 0;
        int aEdge = 0;
        while(gEdge != int(neighbors.size()))
        {
          if(dofs == 1)
          {
            if(neighbors[gEdge] == indices[aEdge])
            {
              //graph and matrix both have this edge, wasn't filtered, show as color 1
              vert1.push_back(pair<int, int>(int(globRow), neighbors[gEdge]));
              gEdge++;
              aEdge++;
            }
            else
            {
              //graph contains an edge at gEdge which was filtered from A, show as color 2
              vert2.push_back(pair<int, int>(int(globRow), neighbors[gEdge]));
              gEdge++;
            }
          }
          else //for multiple DOF problems, don't try to detect filtered edges and ignore A
          {
            vert1.push_back(pair<int, int>(int(globRow), neighbors[gEdge]));
            gEdge++;
          }
        }
      }
    }
    else
    {
      ArrayView<const LocalOrdinal> indices;
      for(LocalOrdinal locRow = 0; locRow < LocalOrdinal(A->getNodeNumRows()); locRow++)
      {
        if(dofs == 1)
          A->getLocalRowView(locRow, indices, values);
        neighbors = G->getNeighborVertices(locRow);
        //Add those local indices (columns) to the list of connections (which are pairs of the form (localM, localN))
        int gEdge = 0;
        int aEdge = 0;
        while(gEdge != int(neighbors.size()))
        {
          if(dofs == 1)
          {
            if(neighbors[gEdge] == indices[aEdge])
            {
              vert1.push_back(pair<int, int>(locRow, neighbors[gEdge]));
              gEdge++;
              aEdge++;
            }
            else
            {
              vert2.push_back(pair<int, int>(locRow, neighbors[gEdge]));
              gEdge++;
            }
          }
          else
          {
            vert1.push_back(pair<int, int>(locRow, neighbors[gEdge]));
            gEdge++;
          }
        }
      }
    }
    for(size_t i = 0; i < vert1.size(); i ++)
    {
      if(vert1[i].first > vert1[i].second)
      {
        int temp = vert1[i].first;
        vert1[i].first = vert1[i].second;
        vert1[i].second = temp;
      }
    }
    for(size_t i = 0; i < vert2.size(); i++)
    {
      if(vert2[i].first > vert2[i].second)
      {
        int temp = vert2[i].first;
        vert2[i].first = vert2[i].second;
        vert2[i].second = temp;
      }
    }
    sort(vert1.begin(), vert1.end());
    vector<pair<int, int> >::iterator newEnd = unique(vert1.begin(), vert1.end()); //remove duplicate edges
    vert1.erase(newEnd, vert1.end());
    sort(vert2.begin(), vert2.end());
    newEnd = unique(vert2.begin(), vert2.end());
    vert2.erase(newEnd, vert2.end());
    vector<int> points1;
    points1.reserve(2 * vert1.size());
    for(size_t i = 0; i < vert1.size(); i++)
    {
      points1.push_back(vert1[i].first);
      points1.push_back(vert1[i].second);
    }
    vector<int> points2;
    points2.reserve(2 * vert2.size());
    for(size_t i = 0; i < vert2.size(); i++)
    {
      points2.push_back(vert2[i].first);
      points2.push_back(vert2[i].second);
    }
    vector<int> unique1 = makeUnique_(points1);
    vector<int> unique2 = makeUnique_(points2);
    fout << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << endl;
    fout << "  <UnstructuredGrid>" << endl;
    fout << "    <Piece NumberOfPoints=\"" << unique1.size() + unique2.size() << "\" NumberOfCells=\"" << vert1.size() + vert2.size() << "\">" << endl;
    fout << "      <PointData Scalars=\"Node Aggregate Processor\">" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Node\" format=\"ascii\">" << endl; //node and aggregate will be set to CONTRAST_1|2, but processor will have its actual value
    string indent = "          ";
    fout << indent;
    for(size_t i = 0; i < unique1.size(); i++)
    {
      fout << CONTRAST_1_ << " ";
      if(i % 25 == 24)
        fout << endl << indent;
    }
    for(size_t i = 0; i < unique2.size(); i++)
    {
      fout << CONTRAST_2_ << " ";
      if((i + 2 * vert1.size()) % 25 == 24)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Aggregate\" format=\"ascii\">" << endl;
    fout << indent;
    for(size_t i = 0; i < unique1.size(); i++)
    {
      fout << CONTRAST_1_ << " ";
      if(i % 25 == 24)
        fout << endl << indent;
    }
    for(size_t i = 0; i < unique2.size(); i++)
    {
      fout << CONTRAST_2_ << " ";
      if((i + 2 * vert2.size()) % 25 == 24)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << endl;
    fout << indent;
    for(size_t i = 0; i < unique1.size() + unique2.size(); i++)
    {
      fout << myRank_ << " ";
      if(i % 25 == 24)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </PointData>" << endl;
    fout << "      <Points>" << endl;
    fout << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    fout << indent;
    for(size_t i = 0; i < unique1.size(); i++)
    {
      if(fine)
      {
        fout << xCoords_[unique1[i]] << " " << yCoords_[unique1[i]] << " ";
        if(dims_ == 3)
          fout << zCoords_[unique1[i]] << " ";
        else
          fout << "0 ";
        if(i % 2)
          fout << endl << indent;
      }
      else
      {
        fout << cx_[unique1[i]] << " " << cy_[unique1[i]] << " ";
        if(dims_ == 3)
          fout << cz_[unique1[i]] << " ";
        else
          fout << "0 ";
        if(i % 2)
          fout << endl << indent;
      }
    }
    for(size_t i = 0; i < unique2.size(); i++)
    {
      if(fine)
      {
        fout << xCoords_[unique2[i]] << " " << yCoords_[unique2[i]] << " ";
        if(dims_ == 3)
          fout << zCoords_[unique2[i]] << " ";
        else
          fout << "0 ";
        if(i % 2)
          fout << endl << indent;
      }
      else
      {
        fout << cx_[unique2[i]] << " " << cy_[unique2[i]] << " ";
        if(dims_ == 3)
          fout << cz_[unique2[i]] << " ";
        else
          fout << "0 ";
        if((i + unique1.size()) % 2)
          fout << endl << indent;
      }
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Points>" << endl;
    fout << "      <Cells>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    fout << indent;
    for(size_t i = 0; i < points1.size(); i++)
    {
      fout << points1[i] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    for(size_t i = 0; i < points2.size(); i++)
    {
      fout << points2[i] + unique1.size() << " ";
      if((i + points1.size()) % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>"  << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    fout << indent;
    int offset = 0;
    for(size_t i = 0; i < vert1.size() + vert2.size(); i++)
    {
      offset += 2;
      fout << offset << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>"  << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
    fout << indent;
    for(size_t i = 0; i < vert1.size() + vert2.size(); i++)
    {
      fout << "3 ";
      if(i % 25 == 24)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>"  << endl;
    fout << "      </Cells>" << endl;
    fout << "    </Piece>" << endl;
    fout << "  </UnstructuredGrid>" << endl;
    fout << "</VTKFile>" << endl;
  }
  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writeFile_(std::ofstream& fout, std::string styleName, std::vector<int>& vertices, std::vector<int>& geomSizes) const
  {
    using namespace std;
    vector<int> uniqueFine = makeUnique_(vertices);
    string indent = "      ";
    fout << "<!--" << styleName << " Aggregates Visualization-->" << endl;
    fout << "<VTKFile type=\"UnstructuredGrid\" byte_order=\"LittleEndian\">" << endl;
    fout << "  <UnstructuredGrid>" << endl;
    fout << "    <Piece NumberOfPoints=\"" << uniqueFine.size() << "\" NumberOfCells=\"" << geomSizes.size() << "\">" << endl;
    fout << "      <PointData Scalars=\"Node Aggregate Processor\">" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Node\" format=\"ascii\">" << endl;
    indent = "          ";
    fout << indent;
    bool localIsGlobal = GlobalOrdinal(nodeMap_->getGlobalNumElements()) == GlobalOrdinal(nodeMap_->getNodeNumElements());
    for(size_t i = 0; i < uniqueFine.size(); i++)
    {
      if(localIsGlobal)
      {
        fout << uniqueFine[i] << " "; //if all nodes are on this processor, do not map from local to global
      }
      else
        fout << nodeMap_->getGlobalElement(uniqueFine[i]) << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Aggregate\" format=\"ascii\">" << endl;
    fout << indent;
    for(size_t i = 0; i < uniqueFine.size(); i++)
    {
      fout << aggsOffset_ + vertex2AggId_[uniqueFine[i]] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"Processor\" format=\"ascii\">" << endl;
    fout << indent;
    for(size_t i = 0; i < uniqueFine.size(); i++)
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
    for(size_t i = 0; i < uniqueFine.size(); i++)
    {
      fout << xCoords_[uniqueFine[i]] << " " << yCoords_[uniqueFine[i]] << " ";
      if(dims_ == 2)
        fout << "0 ";
      else
        fout << zCoords_[uniqueFine[i]] << " ";
      if(i % 3 == 2)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Points>" << endl;
    fout << "      <Cells>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
    fout << indent;
    for(size_t i = 0; i < vertices.size(); i++)
    {
      fout << vertices[i] << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
    fout << indent;
    int accum = 0;
    for(size_t i = 0; i < geomSizes.size(); i++)
    {
      accum += geomSizes[i];
      fout << accum << " ";
      if(i % 10 == 9)
        fout << endl << indent;
    }
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "        <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << endl;
    fout << indent;
    for(size_t i = 0; i < geomSizes.size(); i++)
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
    fout << endl;
    fout << "        </DataArray>" << endl;
    fout << "      </Cells>" << endl;
    fout << "    </Piece>" << endl;
    fout << "  </UnstructuredGrid>" << endl;
    fout << "</VTKFile>" << endl;
  } 
  
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::vector<int> AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::makeUnique_(std::vector<int>& vertices) const
  {
    using namespace std;
    vector<int> uniqueNodes = vertices;
    sort(uniqueNodes.begin(), uniqueNodes.end());
    vector<int>::iterator newUniqueFineEnd = unique(uniqueNodes.begin(), uniqueNodes.end());
    uniqueNodes.erase(newUniqueFineEnd, uniqueNodes.end());
    //uniqueNodes is now a sorted list of the nodes whose info actually goes in file
    //Now replace values in vertices with locations of the old values in uniqueFine
    for(size_t i = 0; i < vertices.size(); i++)
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
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::buildColormap_() const
  {
    using namespace std;
    try
    {
      ofstream color("random-colormap.xml");
      color << "<ColorMap name=\"MueLu-Random\" space=\"RGB\">" << endl;
      //Give -1, -2, -3 distinctive colors (so that the style functions can have constrasted geometry types)
      //Do red, orange, yellow to constrast with cool color spectrum for other types
      color << "  <Point x=\"" << CONTRAST_1_ << "\" o=\"1\" r=\"1\" g=\"0\" b=\"0\"/>" << endl;
      color << "  <Point x=\"" << CONTRAST_2_ << "\" o=\"1\" r=\"1\" g=\"0.6\" b=\"0\"/>" << endl;
      color << "  <Point x=\"" << CONTRAST_3_ << "\" o=\"1\" r=\"1\" g=\"1\" b=\"0\"/>" << endl;
      srand(time(NULL));
      for(int i = 0; i < 5000; i += 4)
      {
        color << "  <Point x=\"" << i << "\" o=\"1\" r=\"" << (rand() % 50) / 256.0 << "\" g=\"" << (rand() % 256) / 256.0 << "\" b=\"" << (rand() % 256) / 256.0 << "\"/>" << endl;
      }
      color << "</ColorMap>" << endl;
      color.close();
    }
    catch(exception& e)
    {
      cout << "   Error while building colormap file: " << e.what() << endl;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void AggregationExportFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::writePVTU_(std::ofstream& pvtu, std::string baseFname, int numProcs) const
  {
    using namespace std;
    //If using vtk, filenameToWrite now contains final, correct ***.vtu filename (for the current proc)
    //So the root proc will need to use its own filenameToWrite to make a list of the filenames of all other procs to put in
    //pvtu file.
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
    for(int i = 0; i < numProcs; i++)
    {
      //specify the piece for each proc (the replaceAll expression matches what the filenames will be on other procs)
      pvtu << "    <Piece Source=\"" << replaceAll(baseFname, "%PROCID", toString(i)) << "\"/>" << endl;
    }
    //reference files with graph pieces, if applicable
    if(doFineGraphEdges_)
    {
      for(int i = 0; i < numProcs; i++)
      {
        string fn = replaceAll(baseFname, "%PROCID", toString(i));
        pvtu << "    <Piece Source=\"" << fn.insert(fn.rfind(".vtu"), "-finegraph") << "\"/>" << endl;
      }
    }
    if(doCoarseGraphEdges_)
    {
      for(int i = 0; i < numProcs; i++)
      {
        string fn = replaceAll(baseFname, "%PROCID", toString(i));
        pvtu << "    <Piece Source=\"" << fn.insert(fn.rfind(".vtu"), "-coarsegraph") << "\"/>" << endl;
      }
    }
    pvtu << "  </PUnstructuredGrid>" << endl;
    pvtu << "</VTKFile>" << endl;
    pvtu.close();
  }
} // namespace MueLu

#endif /* MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_ */
