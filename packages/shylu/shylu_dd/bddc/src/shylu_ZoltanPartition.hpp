
//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef BDDC_ZOLTANPARTITION_H
#define BDDC_ZOLTANPARTITION_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <cmath>
#include <mpi.h>
#include "Tpetra_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"  
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "shylu_errorBDDC.hpp"
#include <zoltan_cpp.h>

using Teuchos::RCP;
using Teuchos::rcp;

enum PartitionOption{
  GRAPH = 0,
  GEOM = 1,
  PHG = 2,
  RCB = 3,
  RIB = 4
};

namespace bddc {

template <class LO,
	  class GO> class ZoltanPartition
{
public:
  //
  // Convenience typedefs
  //
  typedef Tpetra::Map<LO,GO>                                 Map;
  typedef Tpetra::CrsGraph<LO,GO>                            CrsGraph;
  typedef Tpetra::Export<LO,GO>                              Export;
  typedef Tpetra::Import<LO,GO>                              Import;
  
  ZoltanPartition(RCP<const CrsGraph> Graph,
		  RCP< Teuchos::ParameterList > Parameters) :
  m_Graph(Graph),
    m_Parameters(Parameters),
    m_Comm(GetMPIComm()),
    m_TComm(Graph->getComm()),
    m_IGO(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid())
  {
    ConstructGraphEntitiesForZoltan();
    Initialize();
    SetParameters();
  }

  ZoltanPartition(LO numRows,
		  const LO* rowBegin,
		  const GO* columnGIDs,
		  const GO* columnPIDs,
		  const double* coords,
		  const GO* globalIDs,
		  MPI_Comm Comm,
		  RCP< Teuchos::ParameterList > Parameters) :
  m_Graph(Teuchos::null),
    m_Parameters(Parameters),
    m_Comm(Comm),
    m_TComm(rcp( new Teuchos::MpiComm<int>(Comm) )),
    m_coords(coords)
  {
    ConstructGraphEntitiesForZoltan(numRows, rowBegin, columnGIDs, columnPIDs,
				    globalIDs);
    Initialize();
    SetParameters();
  }

  ~ZoltanPartition()
  {
  }

  const int* getParts() const
  {
    return m_parts.data();
  }

  void setDefaultWeights()
  {
    LO numEntries(0);
    if (m_LBMethod == GRAPH) {
      numEntries = m_ConnBegin[m_numMyObjects];
      if (m_edgeWeights.size() == 0) {
	m_edgeWeights.resize(numEntries, 1); 
      }
      m_vertexWeights.resize(m_numMyObjects, 1);
    }
    else {
      if (m_vertexWeights.size() == 0) {
	m_vertexWeights.resize(m_numMyObjects, 1);
      }
    }
  }

  void setVertexWeights(std::vector<float> & vertexWeights)
  {
    int numEntries = vertexWeights.size();
    BDDC_TEST_FOR_EXCEPTION(numEntries != m_numMyObjects, std::runtime_error, 
			    "numEntries inconsistent with m_numMyObjects");
    m_vertexWeights = vertexWeights;
  }

  void setEdgeWeights(std::vector<float> & edgeWeights)
  {
    if (m_LBMethod == GRAPH) {
      BDDC_TEST_FOR_EXCEPTION(edgeWeights.size() != m_ConnBegin[m_numMyObjects],
			      std::runtime_error, "edgeWeights size error");
    }
    m_edgeWeights = edgeWeights;
  }

  void doPartition()
  {
    int changes;
    int numGidEntries;
    int numLidEntries;
    int numImport;
    ZOLTAN_ID_PTR importGlobalIds;
    ZOLTAN_ID_PTR importLocalIds;
    int *importProcs;
    int *importToPart;
    int numExport;
    ZOLTAN_ID_PTR exportGlobalIds;
    ZOLTAN_ID_PTR exportLocalIds;
    int *exportProcs;
    int *exportToPart;
    
    int rc = m_zz->LB_Partition
      (changes, numGidEntries, numLidEntries,
       numImport, importGlobalIds, importLocalIds, importProcs, importToPart,
       numExport, exportGlobalIds, exportLocalIds, exportProcs, exportToPart);
    
    if (rc != ZOLTAN_OK) {
      printf("Partitioning failed on process %d\n", m_MyPID);
      delete m_zz;
      return;
    }

    m_parts.resize(m_numMyObjects, m_MyPID);
    for (int i=0; i<numExport; i++) {
      //      m_parts[i] = exportToPart[exportLocalIds[i]]; yada-like
      //      m_parts[exportLocalIds[i]] = exportToPart[i]; for RETURN_LISTS = PARTS
      const int index = exportGlobalIds[i] - m_objectGIDs[0];
      BDDC_TEST_FOR_EXCEPTION(index < 0, std::runtime_error, 
			      "index less than zero");
      BDDC_TEST_FOR_EXCEPTION(index >= m_numMyObjects, std::runtime_error, 
			      "index too large");
      m_parts[index] = exportToPart[i];
    }

    delete m_zz;
    Zoltan::LB_Free_Part(&importGlobalIds, &importLocalIds, &importProcs,
			 &importToPart);
    Zoltan::LB_Free_Part(&exportGlobalIds, &exportLocalIds, &exportProcs,
			 &exportToPart);
  }

  LO getNumMyObjects() 
  {
    return m_numMyObjects;
  }

  const double* getCoords()
  {
    return m_coords;
  }

  GO getObjectGIDs(int i) {
    return m_objectGIDs[i];
  }

  GO* getObjectGIDs() {
    &m_objectGIDs[0];
  }

  void getConnectivityBegin(int* & ConnBegin) {
    ConnBegin = &m_ConnBegin[0];
  }

  void getConnectivityGlobal(GO* & ConnGlobal) {
    ConnGlobal = &m_ConnGlobal[0];
  }

  void getConnectivityOwner(int* & ConnOwners) {
    ConnOwners = &m_ConnOwners[0];
  }

  void getEdgeWeights(float* & edgeWeights) {
    edgeWeights = &m_edgeWeights[0];
  }

  void getVertexWeights(float* & vertexWeights) {
    vertexWeights = &m_vertexWeights[0];
  }
  
  static int getNumberOfObjects(void *data, 
				int *ierr)
  {
    ZoltanPartition *objs = (ZoltanPartition *)data;
    *ierr = ZOLTAN_OK;
    return objs->getNumMyObjects();
  }

  static void getObjectList(void *data, 
			    int sizeGID, 
			    int sizeLID,
			    ZOLTAN_ID_PTR globalID, 
			    ZOLTAN_ID_PTR localID,
			    int wgt_dim, 
			    float *obj_wgts, 
			    int *ierr)
  {
    if ((sizeGID != 1) || (sizeLID != 1) || wgt_dim != 1) {
      *ierr = ZOLTAN_FATAL;
      return;
    }
    ZoltanPartition *objs = (ZoltanPartition *)data;
    int n = objs->getNumMyObjects();
    float* vertexWeights;
    objs->getVertexWeights(vertexWeights);
    for (int i=0; i<n; i++) {
      globalID[i] = objs->getObjectGIDs(i);
      localID[i] = i;
      obj_wgts[i] = vertexWeights[i];
    }
    *ierr = ZOLTAN_OK;
  }

  static int getSpatialDimension(void *data, 
				 int *ierr)
  {
    *ierr = ZOLTAN_OK;
    return 3;
  }

  static void getCoordinates(void *data, 
			     int sizeGID, 
			     int sizeLID,
			     int num_obj,
			     ZOLTAN_ID_PTR globalID, 
			     ZOLTAN_ID_PTR localID,
			     int num_dim,
			     double* coordinates,
			     int *ierr)
  {
    if ((sizeGID != 1) || (sizeLID != 1) || (num_dim != 3)) {
      *ierr = ZOLTAN_FATAL;
      return;
    }
    ZoltanPartition *objs = (ZoltanPartition *)data;
    const double* xyz = objs->getCoords();
    int numMyObjects = objs->getNumMyObjects();
    for (int i=0; i<num_obj; i++) {
      for (int j=0; j<num_dim; j++) {
	coordinates[num_dim*i+j] = xyz[j*numMyObjects+localID[i]];
      }
    }
    *ierr = ZOLTAN_OK;
  }

  static void getNumEdges(void *data, 
			  int sizeGID, 
			  int sizeLID,
			  int num_obj,
			  ZOLTAN_ID_PTR globalID, 
			  ZOLTAN_ID_PTR localID,
			  int *numEdges, 
			  int *ierr)
  {
    ZoltanPartition *objs = (ZoltanPartition *)data;
    int n = objs->getNumMyObjects();
    if ((sizeGID != 1) || (sizeLID != 1) || (num_obj != n)) {
      *ierr = ZOLTAN_FATAL;
      return;
    }
    int *Conn2(0);
    objs->getConnectivityBegin(Conn2);
    for (int i=0; i<n; i++){
      int idx = localID[i];
      numEdges[i] = Conn2[idx+1] - Conn2[idx];
    }
    *ierr = ZOLTAN_OK;
  }

  static void getEdgeList(void *data, 
			  int sizeGID, 
			  int sizeLID,
			  int num_obj, 
			  ZOLTAN_ID_PTR globalID, 
			  ZOLTAN_ID_PTR localID,
			  int *num_edges,
			  ZOLTAN_ID_PTR nborGID, 
			  int *nborProc,
			  int wgt_dim, 
			  float *ewgts, 
			  int *ierr)
  {
    ZoltanPartition *objs = (ZoltanPartition *)data;
    int n = objs->getNumMyObjects();
    if ((sizeGID != 1) || (sizeLID != 1) || (num_obj != n) || (wgt_dim != 1)) {
      *ierr = ZOLTAN_FATAL;
      return;
    }
    *ierr = ZOLTAN_OK;
    GO *ConnectivityGlobal(0);
    LO *ConnectivityOwners(0), *Conn2(0);
    objs->getConnectivityGlobal(ConnectivityGlobal);
    objs->getConnectivityOwner(ConnectivityOwners);
    objs->getConnectivityBegin(Conn2);
    float* edgeWeights;
    objs->getEdgeWeights(edgeWeights);
    for (int i=0; i<num_obj; i++) {
      int idx = localID[i];
      for (int j=Conn2[idx]; j<Conn2[idx+1]; j++) {
	nborGID[j] = ConnectivityGlobal[j];
	nborProc[j] = ConnectivityOwners[j];
	ewgts[j] = edgeWeights[j];
      }
    }
  }

private:
  void ConstructGraphEntitiesForZoltan()
  {
    double* nullDouble(0);
    m_coords = m_Parameters->get("Coordinates", nullDouble);
    int numMyRows = m_Graph->getNodeNumRows();
    int numEntries = m_Graph->getNodeNumEntries();
    RCP<const Map> RowMap = m_Graph->getRowMap();
    RCP<const Map> ColMap = m_Graph->getColMap();
    int numCols = ColMap->getNodeNumElements();
    std::vector<int> colLIDs(numCols), colPIDs(numCols);
    Teuchos::ArrayView<LO> ColLIDs(colLIDs);
    Teuchos::ArrayView<LO> ColPIDs(colPIDs);
    Teuchos::ArrayView<const GO> RowGIDs = RowMap->getNodeElementList();
    Teuchos::ArrayView<const GO> ColGIDs = ColMap->getNodeElementList();
    RowMap->getRemoteIndexList(ColGIDs, ColPIDs, ColLIDs);
    Teuchos::ArrayView<const LO> Indices;
    m_objectGIDs.resize(numMyRows);
    m_ConnOwners.resize(numEntries);
    m_ConnGlobal.resize(numEntries);
    m_ConnBegin.resize(numMyRows+1);
    numEntries = 0;
    std::vector<LO> count(numMyRows);
    int nodeNotOnAnyProcessors(-1);
    for (int i=0; i<numMyRows; i++) {
      m_objectGIDs[i] = RowMap->getNodeElementList()[i];
      m_Graph->getLocalRowView(i, Indices);
      for (int j=0; j<Indices.size(); j++) {
	BDDC_TEST_FOR_EXCEPTION(ColPIDs[Indices[j]] == nodeNotOnAnyProcessors, 
				std::runtime_error, "invalid ColPID");
	m_ConnOwners[numEntries] = ColPIDs[Indices[j]];
	m_ConnGlobal[numEntries] = ColGIDs[Indices[j]];
	numEntries++;
      }
      m_ConnBegin[i+1] = numEntries;
      count[i] = Indices.size();
    }
    m_numMyObjects = numMyRows;
  }

  void ConstructGraphEntitiesForZoltan(const LO numRows,
				       const LO* rowBegin,
				       const GO* columnGIDs,
				       const GO* columnPIDs,
				       const GO* globalIDs)
  {
    LO numMyRows = numRows;
    m_objectGIDs.resize(numMyRows);
    for (LO i=0; i<numMyRows; i++) {
      m_objectGIDs[i] = globalIDs[i];
    }
    m_numMyObjects = numMyRows;
    if (m_Parameters->get("LB Method", "Graph") == "Graph") {
      LO numEntries = rowBegin[numRows];
      m_ConnOwners.resize(numEntries);
      m_ConnGlobal.resize(numEntries);
      m_ConnBegin.resize(numMyRows+1);
      numEntries = 0;
      for (LO i=0; i<numMyRows; i++) {
	for (LO j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	  m_ConnOwners[numEntries] = columnPIDs[j];
	  m_ConnGlobal[numEntries] = columnGIDs[j];
	  numEntries++;
	}
	m_ConnBegin[i+1] = numEntries;
      }
    }
  }

  void Initialize()
  {
    float version;
    int argc(0);
    char **argv(0);
    int err = Zoltan_Initialize(argc, argv, &version);
    BDDC_TEST_FOR_EXCEPTION(err != ZOLTAN_OK, std::runtime_error, 
			    "Zoltan_Initialize error");
    m_zz = new Zoltan(m_Comm);
    m_MyPID = m_TComm->getRank();
    m_numProc = m_TComm->getSize();
    /*
    if (m_MyPID == 0) {
      std::cout << "Zoltan version = " << version << std::endl;
    } 
    */
    m_LBMethod = GRAPH;
    std::string aString = m_Parameters->get("LB Method", "Graph");
    if (aString == "RCB") m_LBMethod = RCB;
    if (aString == "RIB") m_LBMethod = RIB;
  }

  void SetParameters()
  {
    //
    // set Zoltan parameters
    //
    //    m_zz->Set_Param("LB_APPROACH", "PARTITION");// partition from scratch
    m_zz->Set_Param("LB_APPROACH", "REPARTITION");
    m_numParts = m_Parameters->get("Number of Parts", 1);
    std::ostringstream convert;
    convert << m_numParts;
    std::string numPartsString = convert.str();
    m_zz->Set_Param("GRAPH_PACKAGE", "PHG");
#ifdef HAVE_SHYLU_DDBDDC_PARMETIS
    const std::string graphPackage = m_Parameters->get("Graph Package", "PHG");
    if (graphPackage == "ParMETIS") {
      m_zz->Set_Param("GRAPH_PACKAGE", "ParMETIS");
    }
#endif
    m_zz->Set_Param("DEBUG_LEVEL", "0");
    m_zz->Set_Param("NUM_GLOBAL_PARTS", numPartsString);
    m_zz->Set_Param("RETURN_LISTS", "ALL");/* get Import & Export info */
    //    m_zz->Set_Param("RETURN_LISTS", "PARTS");   /* get Parts info */
    m_zz->Set_Param("NUM_GID_ENTRIES", "1");  /* global ID is 1 integer */
    m_zz->Set_Param("NUM_LID_ENTRIES", "1");  /* local ID is 1 integer */
    m_zz->Set_Param("OBJ_WEIGHT_DIM", "1");   /* weight is 1 float */
    m_zz->Set_Param("EDGE_WEIGHT_DIM", "1");  /* edge weights */
    //    m_zz->Set_Param("IMBALANCE_TOL", "1.05"); /* default is 1.1 */
    //  m_zz->Set_Param("IMBALANCE_TOL", "1.1"); /* default is 1.1 */
    switch (m_LBMethod) {
    case GRAPH:
      m_zz->Set_Param("LB_METHOD", "GRAPH"); // standard graph partitioning
      break;
    case RCB:
      m_zz->Set_Param("LB_METHOD", "RCB"); // recursive coordinate bisection
      break;
    case RIB:
      m_zz->Set_Param("LB_METHOD", "RIB"); // recursive inertial bisection
      break;
    default:
      std::cout << "Error: LB_APPROACH not supported\n";
      break;
    }
    //
    // register Zoltan query functions 
    //
    m_zz->Set_Num_Obj_Fn(ZoltanPartition::getNumberOfObjects, this);
    m_zz->Set_Obj_List_Fn(ZoltanPartition::getObjectList, this);
    m_zz->Set_Num_Geom_Fn(ZoltanPartition::getSpatialDimension, this);
    m_zz->Set_Geom_Multi_Fn(ZoltanPartition::getCoordinates, this);
    m_zz->Set_Num_Edges_Multi_Fn(ZoltanPartition::getNumEdges, this);
    m_zz->Set_Edge_List_Multi_Fn(ZoltanPartition::getEdgeList, this);
  }

  MPI_Comm GetMPIComm()
  {
    RCP<const Teuchos::Comm<int> > comm = m_Graph->getComm();
    const Teuchos::MpiComm<int>* mpiComm = 
      dynamic_cast<const Teuchos::MpiComm<int>* > (&(*comm));
    return *(mpiComm->getRawMpiComm());
  }

 private: // member variables
  RCP<const CrsGraph> m_Graph;
  RCP<Teuchos::ParameterList> m_Parameters;
  MPI_Comm m_Comm;
  RCP<const Teuchos::Comm<int> > m_TComm;
  RCP<CrsGraph> m_EdgeNodeGraph, m_NodeEdgeGraph, m_ElementNodeGraph;
  Tpetra::global_size_t m_IGO;
  enum PartitionOption m_inputGraphType, m_LBMethod, m_graphPackage,
    m_outputOption;
  int m_numParts{0}, m_numMyObjects{0}, m_numMyHyperEdges{0}, m_MyPID{-1}, 
     m_numProc{-1};
  std::vector<LO> m_ConnOwners, m_ConnBegin;
  std::vector<GO> m_objectGIDs, m_hyperEdgeGIDs, m_ConnGlobal,
    m_subdomainNodeGIDs;
  std::vector<float> m_edgeWeights, m_vertexWeights;
  std::vector<int> m_parts;
  bool m_debug;
  Zoltan *m_zz;
  const double* m_coords;
};

} // namespace bddc

#endif // BDDC_ZOLTANPARTITION_H
