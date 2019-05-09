
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef BDDC_PROBLEMMAKER_H
#define BDDC_PROBLEMMAKER_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <time.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <complex>

#ifdef HAVE_SHYLU_DDBDDC_METIS
#include "metis.h"
#endif

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

#include "shylu_enumsBDDC.hpp"
#include "shylu_errorBDDC.hpp"

#ifdef _OPENMP
#include "omp.h"
#endif

//#define TENSORPRODUCTBDDC
#ifdef TENSORPRODUCTBDDC
#include "TensorProductBDDC.hpp"
#endif

using Teuchos::RCP;
using Teuchos::rcp;

// Author: Clark R. Dohrmann
namespace bddc {

template <class LO, class GO, class SX, class SM> class ProblemMaker
{
 public:
  ProblemMaker(RCP<Teuchos::ParameterList> Parameters,
	       MPI_Comm Comm) :
  m_Parameters(Parameters),
    m_Comm(Comm),
    m_problemType(Parameters->get("Problem Type", SCALARPDE)),
    m_analysisType(Parameters->get("Analysis Type", STANDARD)),
    m_spatialDim(Parameters->get("Spatial Dimension", 1)),
    m_readModel(false)
  {
    MPI_Comm_size(Comm, &m_numProc);
    MPI_Comm_rank(Comm, &m_myPID);
    GenerateModel();
    ApplyEssentialBCs();
    setQuadraturePoints();
  }

  ProblemMaker(RCP<Teuchos::ParameterList> Parameters,
	       MPI_Comm Comm,
	       const std::vector<SM> & xFile,
	       const std::vector<SM> & yFile,
	       const std::vector<SM> & zFile,
	       const std::vector< std::vector<LO> > & elemConn,
	       const std::vector<GO> & nodeGIDs,
	       const std::vector< std::vector<LO> > & subElems,
	       const std::vector<LO> & constrainedNodes) :
  m_Parameters(Parameters),
    m_Comm(Comm),
    m_problemType(Parameters->get("Problem Type", SCALARPDE)),
    m_analysisType(Parameters->get("Analysis Type", STANDARD)),
    m_spatialDim(Parameters->get("Spatial Dimension", 1)),
    m_numNode(nodeGIDs.size()),
    m_numElem(elemConn.size()),
    m_x(xFile),
    m_y(yFile),
    m_z(zFile),
    m_nodeGIDs(nodeGIDs),
    m_elemConn(elemConn),
    m_readModel(true)
  {
    MPI_Comm_size(Comm, &m_numProc);
    MPI_Comm_rank(Comm, &m_myPID);
    determineNodeBegLocalDof();
    m_subElems = subElems;
    ApplyEssentialBCs(constrainedNodes);
    setQuadraturePoints();
  }

  ~ProblemMaker()
  {
#ifdef TENSORPRODUCTBDDC
    delete m_TP;
#endif
  }

  LO getNumNode() const {
    return m_numNode;
  }
  
  LO getNumElem() const {
    return m_numElem;
  }
  
  LO getNumDof() const {
    return m_nodeBegin[m_numNode];
  }
  
  int getNumDofPerNode() const
  {
    switch (m_problemType) {
    case SCALARPDE:
      return 1;
      break;
    case ELASTICITY:
      BDDC_TEST_FOR_EXCEPTION(m_spatialDim == 1, std::runtime_error, 
			      "spatialDim = 1 invalid for elasticity");
      return m_spatialDim;
      break;
    default:
      return 1;
      break;
    }
  }

  inline int power2(const int n) const
  {
    int retval = 1 << n;
    return retval;
  }

  const LO* getNodeBegin() const {
    return m_nodeBegin.data();
  }

  const LO* getLocalDofs() const {
    return m_localDofs.data();
  }

  const SM* getXcoord() const {
    return m_x.data();
  }

  const SM* getYcoord() const {
    return m_y.data();
  }

  const SM* getZcoord() const {
    return m_z.data();
  }

  const GO* getNodeGlobalIDs() const {
    return m_nodeGIDs.data();
  }

  const SX* getLoadVector() const {
    return m_loadVector.data();
  }

  void printElementMatrices(const int elem=0) const
  {
    std::vector<SM> elemStiffMatrix, elemMassMatrix;
    const double E(1.0), nu(0.3);
    getElementMatrices(m_elemConn[elem], E, nu, elemStiffMatrix, elemMassMatrix);
    printElementMatrix(elemStiffMatrix.data(), "stiffMatrix.dat");
    printElementMatrix(elemMassMatrix.data(), "massMatrix.dat");
  }

  void getElementMatrices(std::vector<SM> & elemStiffMatrix,
			  std::vector<SM> & elemMassMatrix,
			  const int elem=0)
  {
    const double E(1.0), nu(0.3);
    getElementMatrices(m_elemConn[elem], E, nu, elemStiffMatrix,
		       elemMassMatrix);
  }

  const std::vector< std::vector<LO> > getElementConnectivity() const
  {
    return m_elemConn;
  }

  SM norm2(const SX* values,
	   const LO numTerm) const
  {
    SM mag(0);
    for (LO i=0; i<numTerm; i++) {
      mag += squared(values[i]);
    }
    return sqrt(mag);
  }

  void getSubDomainElements(LO numSubDir1PerProc, 
			    LO numSubDir2PerProc,
			    LO numSubDir3PerProc, 
			    std::vector< std::vector<LO> > & subElems,
			    std::vector<LO> & subIndices) const
  {
    if (m_x.size() == 0) return;
    SM minX(m_x[0]), maxX(m_x[0]);
    SM minY(m_y[0]), maxY(m_y[0]);
    SM minZ(m_z[0]), maxZ(m_z[0]);
    for (LO i=1; i<m_numNode; i++) {
      if (m_x[i] < minX) minX = m_x[i];
      if (m_x[i] > maxX) maxX = m_x[i];
      if (m_y[i] < minY) minY = m_y[i];
      if (m_y[i] > maxY) maxY = m_y[i];
      if (m_z[i] < minZ) minZ = m_z[i];
      if (m_z[i] > maxZ) maxZ = m_z[i];
    }
    SM Hx = (maxX - minX)/numSubDir1PerProc;
    SM Hy = (maxY - minY)/numSubDir2PerProc;
    SM Hz = (maxZ - minZ)/numSubDir3PerProc;
    LO numElems = m_elemConn.size();
    std::vector<SM> xElem(numElems), yElem(numElems), zElem(numElems);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (LO i=0; i<numElems; i++) {
      SM sumX(0), sumY(0), sumZ(0);
      LO numNodePerElem = m_elemConn[i].size();
      for (LO j=0; j<numNodePerElem; j++) {
	LO node = m_elemConn[i][j];
	sumX += m_x[node];
	sumY += m_y[node];
	sumZ += m_z[node];
      }
      xElem[i] = sumX/numNodePerElem;
      yElem[i] = sumY/numNodePerElem;
      zElem[i] = sumZ/numNodePerElem;
    }
    LO numSub = numSubDir1PerProc;
    if (m_spatialDim > 1) numSub *= numSubDir2PerProc;
    if (m_spatialDim == 3) numSub *= numSubDir3PerProc;
    std::vector<LO> count(numSub, 0);
    subIndices.resize(3*numSub, -1);
    for (LO i=0; i<numElems; i++) {
      LO ix = (xElem[i] - minX)/Hx;
      LO iy = (yElem[i] - minY)/Hy;
      LO iz = (zElem[i] - minZ)/Hz;
      if (m_spatialDim == 1) {
	LO sub = ix;
	subIndices[3*sub] = ix;
	count[sub]++;
      }
      else if (m_spatialDim == 2) {
	LO sub = ix + numSubDir1PerProc*iy;
	subIndices[3*sub+0] = ix;
	subIndices[3*sub+1] = iy;
	count[sub]++;
      }
      else if (m_spatialDim == 3) {
	LO sub = ix + numSubDir1PerProc*iy +
	  numSubDir1PerProc*numSubDir2PerProc*iz;
	subIndices[3*sub+0] = ix;
	subIndices[3*sub+1] = iy;
	subIndices[3*sub+2] = iz;
	count[sub]++;
      }
    }
    subElems.resize(numSub);
    for (LO i=0; i<numSub; i++) {
      subElems[i].resize(count[i]);
      count[i] = 0;
    }
    for (LO i=0; i<numElems; i++) {
      LO ix = (xElem[i] - minX)/Hx;
      LO iy = (yElem[i] - minY)/Hy;
      LO iz = (zElem[i] - minZ)/Hz;
      if (m_spatialDim == 1) {
	LO sub = ix;
	subElems[sub][count[sub]++] = i;
      }
      else if (m_spatialDim == 2) {
	LO sub = ix + numSubDir1PerProc*iy;
	subElems[sub][count[sub]++] = i;
      }
      else if (m_spatialDim == 3) {
	LO sub = ix + numSubDir1PerProc*iy +
	  numSubDir1PerProc*numSubDir2PerProc*iz;
	subElems[sub][count[sub]++] = i;
      }
    }
  }

  void getSubDomainElements(LO numPartsPerSub,
			    std::vector< std::vector<LO> > & subElems) const
  {
#ifdef HAVE_SHYLU_DDBDDC_METIS
    // determine elements for each node
    std::vector< std::vector<LO> > nodeElems(m_numNode);
    LO numElems = m_elemConn.size();
    for (LO i=0; i<numElems; i++) {
      for (size_t j=0; j<m_elemConn[i].size(); j++) {
	LO node = m_elemConn[i][j];
	nodeElems[node].push_back(i);
      }
    }
    // determine adjacent elements
    std::vector<bool> elemFlag(numElems, false), nodeFlag(m_numNode, false);
    std::vector< std::vector<int> > adjElems(numElems);
    for (LO i=0; i<numElems; i++) {
      elemFlag[i] = true;
      for (size_t j=0; j<m_elemConn[i].size(); j++) {
	nodeFlag[m_elemConn[i][j]] = true;
      }
      LO numNodesPerFace = 2;
      if (m_spatialDim == 3) numNodesPerFace = 4;
      for (size_t j=0; j<m_elemConn[i].size(); j++) {
	LO node = m_elemConn[i][j];
	for (size_t k=0; k<nodeElems[node].size(); k++) {
	  LO elem = nodeElems[node][k];
	  if (elemFlag[elem] == false) {
	    LO numCommonNodes = 0;
	    for (size_t m=0; m<m_elemConn[elem].size(); m++) {
	      LO node2 = m_elemConn[elem][m];
	      if (nodeFlag[node2] == true) numCommonNodes++;
	    }
	    if (numCommonNodes == numNodesPerFace) {
	      adjElems[i].push_back(elem);
	      elemFlag[elem] = true;
	    }
	  }
	}
      }
      for (size_t j=0; j<m_elemConn[i].size(); j++) {
	nodeFlag[m_elemConn[i][j]] = false;
      }
      for (size_t j=0; j<adjElems[i].size(); j++) {
	elemFlag[adjElems[i][j]] = false;
      }
      elemFlag[i] = false;
    }
    // partition elements using Metis
    LO numTerms = 0;
    for (LO i=0; i<numElems; i++) numTerms += adjElems[i].size();
    std::vector<idx_t> C1(numTerms), C2(numElems+1, 0);
    numTerms = 0;
    for (LO i=0; i<numElems; i++) {
      for (size_t j=0; j<adjElems[i].size(); j++) {
	C1[numTerms++] = adjElems[i][j];
      }
      C2[i+1] = numTerms;
    }
    std::vector<idx_t> parts(numElems);
    idx_t* xadj   = C2.data();
    idx_t* adjncy = C1.data();
    int nparts = numPartsPerSub;
    idx_t* part = parts.data();
    int edgecut, ncon(1);
    std::vector<idx_t> options(METIS_NOPTIONS, 0);
    METIS_SetDefaultOptions(options.data());
    options[METIS_OPTION_NUMBERING] = 0;
    idx_t *vwgt(0), *vsize(0), *adjwgt(0);
    std::vector<real_t> ubvec(ncon, 1.10);
    std::vector<real_t> tpwgts(ncon*nparts, 1.0/nparts);
    int max_sub = 8;
    if (nparts < max_sub) {
      METIS_PartGraphRecursive(&numElems, &ncon, xadj, adjncy, vwgt, vsize, 
			       adjwgt, &nparts, &tpwgts[0], &ubvec[0], 
			       &options[0], &edgecut, part);
    }
    else {
      METIS_PartGraphKway(&numElems, &ncon, xadj, adjncy, vwgt, vsize, 
			  adjwgt, &nparts, &tpwgts[0], &ubvec[0], 
			  &options[0], &edgecut, part);
    }
    // assign elements to subdomains
    subElems.resize(numPartsPerSub);
    for (LO i=0; i<numElems; i++) {
      subElems[part[i]].push_back(i);
    }
#endif
  }

  void getSubDomainNodeData(std::vector< std::vector<LO> > & subNodes, 
			    std::vector< std::vector<LO> > & subNodeBegin,
			    std::vector< std::vector<LO> > & subLocalDofs) const
  {
    LO numSub = m_subElems.size();
    subNodes.resize(numSub);
    subNodeBegin.resize(numSub);
    subLocalDofs.resize(numSub);
    std::vector<LO> countNode(numSub, 0), countDof(numSub, 0);
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
    std::vector<LO> nodeFlag(m_numNode, -1), activeNodes(m_numNode);
#ifdef _OPENMP
#pragma omp for
#endif
    for (LO i=0; i<numSub; i++) {
      LO numElem = m_subElems[i].size();
      LO numNode(0), numDof(0);
      for (LO j=0; j<numElem; j++) {
	LO elem = m_subElems[i][j];
	for (size_t k=0; k<m_elemConn[elem].size(); k++) {
	  LO node = m_elemConn[elem][k];
	  if (nodeFlag[node] == -1) {
	    nodeFlag[node] = numNode;
	    activeNodes[numNode] = node;
	    numDof += m_nodeBegin[node+1] - m_nodeBegin[node];
	    numNode++;
	  }
	}
      }
      for (LO j=0; j<numNode; j++) nodeFlag[activeNodes[j]] = -1;
      countNode[i] = numNode;
      countDof[i] = numDof;
    }
    }
    for (LO i=0; i<numSub; i++) {
      subNodes[i].resize(countNode[i]);
      subNodeBegin[i].resize(countNode[i]+1, 0);
      subLocalDofs[i].resize(countDof[i]);
    }
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
    std::vector<LO> nodeFlag(m_numNode, -1), activeNodes(m_numNode);
#ifdef _OPENMP
#pragma omp for
#endif
    for (LO i=0; i<numSub; i++) {
      LO numElem = m_subElems[i].size();
      LO numNode(0), numDof(0);
      for (LO j=0; j<numElem; j++) {
	LO elem = m_subElems[i][j];
	for (size_t k=0; k<m_elemConn[elem].size(); k++) {
	  LO node = m_elemConn[elem][k];
	  if (nodeFlag[node] == -1) {
	    subNodes[i][numNode] = node;
	    nodeFlag[node] = numNode;
	    activeNodes[numNode] = node;
	    for (LO m=m_nodeBegin[node]; m<m_nodeBegin[node+1]; m++) {
	      subLocalDofs[i][numDof] = m_localDofs[m];
	      numDof++;
	    }
	    subNodeBegin[i][numNode+1] = numDof;
	    numNode++;
	  }
	}
      }
      for (LO j=0; j<numNode; j++) nodeFlag[activeNodes[j]] = -1;
    }
    }
    const bool sortNodes = m_Parameters->get("Sort Subdomain Nodes", false);
    if (sortNodes || (numSub == 1)) {
      // reorder nodes in ascending order
      for (LO i=0; i<numSub; i++) {
	std::vector<LO> nodes = subNodes[i];
	std::sort(nodes.begin(), nodes.end());
	const LO numNode = nodes.size();
	LO numDof(0);
	for (LO j=0; j<numNode; j++) {
	  const LO node = nodes[j];
	  subNodes[i][j] = node;
	  for (LO k=m_nodeBegin[node]; k<m_nodeBegin[node+1]; k++) {
	    subLocalDofs[i][numDof++] = m_localDofs[k];
	  }
	  subNodeBegin[i][j+1] = numDof;
	}
      }
    }
  }

  void adjustNodalCoords(const std::vector< std::vector<LO> > & subNodes)
  {
#ifdef TENSORPRODUCTBDDC
    const bool adjustNodalCoords = 
      m_Parameters->get("Adjust Nodal Coordinates", false);
    if (adjustNodalCoords == false) return;

    LO numElemPerSubDir1 = 
      m_Parameters->get("Num Elems Per Sub Dir 1", 1);
    LO numElemPerSubDir2 = 
      m_Parameters->get("Num Elems Per Sub Dir 2", 1);
    LO numElemPerSubDir3 = 
      m_Parameters->get("Num Elems Per Sub Dir 3", 1);
    if (m_spatialDim >= 2) {
      BDDC_TEST_FOR_EXCEPTION
	(numElemPerSubDir2 != numElemPerSubDir1, std::runtime_error, 
	 "numElems in each direction must be same for tensor product");
    }
    if (m_spatialDim == 3) {
      BDDC_TEST_FOR_EXCEPTION
	(numElemPerSubDir3 != numElemPerSubDir1, std::runtime_error, 
	 "numElems in each direction must be same for tensor product");
    }
    const int numSub = subNodes.size();
    std::vector<double> xNew(m_numNode), yNew(m_numNode), zNew(m_numNode);
    const int p = numElemPerSubDir1;
    m_TP = new TensorProduct(p);
    const std::vector<double> & xGL = m_TP->getGaussLobattoPoints();
    for (int i=0; i<numSub; i++) {
      double xMin(0), xMax(0), yMin(0), yMax(0), zMin(0), zMax(0);
      getMinMaxCoords(subNodes[i], xMin, xMax, yMin, yMax, zMin, zMax);
      const int numNode = subNodes[i].size();
      std::vector<int> xIndex(numNode), yIndex(numNode), zIndex(numNode);
      const double hX = (xMax - xMin)/numElemPerSubDir1;
      const double hY = (yMax - yMin)/numElemPerSubDir2;
      const double hZ = (zMax - zMin)/numElemPerSubDir3;
      for (int j=0; j<numNode; j++) {
	const int node = subNodes[i][j];
	double x = m_x[node] - xMin + 0.0001*hX;
	double y = m_y[node] - yMin + 0.0001*hY;
	double z = m_z[node] - zMin + 0.0001*hZ;
	xIndex[j] = int(x/hX);
	if (hY > 0) yIndex[j] = int(y/hY);
	if (hZ > 0) zIndex[j] = int(z/hZ);
	xNew[node] = xMin + (xGL[xIndex[j]] + 1)/2 * (xMax - xMin);
	yNew[node] = yMin + (xGL[yIndex[j]] + 1)/2 * (yMax - yMin);
	zNew[node] = zMin + (xGL[zIndex[j]] + 1)/2 * (zMax - zMin);
      }
    }
    m_x = xNew; m_y = yNew; m_z = zNew;
#endif
  }
  
#ifdef TENSORPRODUCTBDDC
  TensorProduct* getTP() {
    return m_TP;
  }
#endif

  void getMinMaxCoords(const std::vector<LO> & nodes, 
		       double & xMin, 
		       double & xMax, 
		       double & yMin, 
		       double & yMax, 
		       double & zMin, 
		       double & zMax)
  {
    const int numNode = nodes.size();
    if (numNode == 0) return;
    xMin = m_x[nodes[0]]; xMax = m_x[nodes[0]];
    yMin = m_y[nodes[0]]; yMax = m_y[nodes[0]];
    zMin = m_z[nodes[0]]; zMax = m_z[nodes[0]];
    for (int i=1; i<numNode; i++) {
      const int node = nodes[i];
      if (m_x[node] < xMin) xMin = m_x[node];
      if (m_x[node] > xMax) xMax = m_x[node];
      if (m_y[node] < yMin) yMin = m_y[node];
      if (m_y[node] > yMax) yMax = m_y[node];
      if (m_z[node] < zMin) zMin = m_z[node];
      if (m_z[node] > zMax) zMax = m_z[node];
    }
  }
  
  void getSubdomainMatrices(std::vector< std::vector<LO> > & subNodes,
			    std::vector< std::vector<LO> > & subNodeBegin,
			    std::vector< std::vector<LO> > & subRowBegin, 
			    std::vector< std::vector<LO> > & subColumns,
			    std::vector< std::vector<SX> > & subValues)
  {
    LO numSub = m_subElems.size();
    std::vector< std::vector<LO> > nodalConn;
    determineNodalConnectivity(m_elemConn, nodalConn);
    // first pass to determine memory requirements
    std::vector<LO> numNonzeros(numSub, 0);
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
    std::vector<LO> nodeMap(m_numNode, -1);
#ifdef _OPENMP
#pragma omp for
#endif
    for (LO sub=0; sub<numSub; sub++) {
      LO numNode = subNodes[sub].size();
      for (LO j=0; j<numNode; j++) nodeMap[subNodes[sub][j]] = j;
      const std::vector<LO> & nodeBegin = subNodeBegin[sub];
      LO nnz(0);
      for (LO i=0; i<numNode; i++) {
	LO node = subNodes[sub][i];
	for (LO j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
	  for (size_t k=0; k<nodalConn[node].size(); k++) {
	    LO node2 = nodeMap[nodalConn[node][k]];
	    if (node2 != -1) {
	      for (LO m=nodeBegin[node2]; m<nodeBegin[node2+1]; m++) {
		nnz++;
	      }
	    }
	  }
	}
      }
      numNonzeros[sub] = nnz;
      for (LO j=0; j<numNode; j++) nodeMap[subNodes[sub][j]] = -1;
    }
    }
    // allocate memory of std::vectors outside of parallel for region
    subRowBegin.resize(numSub);
    subColumns.resize(numSub);
    subValues.resize(numSub);
    LO maxNumDof(0);
    for (LO sub=0; sub<numSub; sub++) {
      LO numNode = subNodes[sub].size();
      const std::vector<LO> & nodeBegin = subNodeBegin[sub];
      LO numDof = nodeBegin[numNode];
      if (numDof > maxNumDof) maxNumDof = numDof;
      subRowBegin[sub].resize(numDof+1, 0);
      subColumns[sub].resize(numNonzeros[sub]);
      subValues[sub].resize(numNonzeros[sub]);
    }
    // determine matrix graph
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
    std::vector<LO> nodeMap(m_numNode, -1);
#ifdef _OPENMP
#pragma omp for
#endif
    for (LO sub=0; sub<numSub; sub++) {
      LO numNode = subNodes[sub].size();
      for (LO j=0; j<numNode; j++) nodeMap[subNodes[sub][j]] = j;
      const std::vector<LO> & nodeBegin = subNodeBegin[sub];
      std::vector<LO> & rowBegin = subRowBegin[sub];
      std::vector<LO> & columns = subColumns[sub];
      LO nnz(0);
      for (LO i=0; i<numNode; i++) {
	LO node = subNodes[sub][i];
	for (LO j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
	  //	  const int start = nnz;
	  for (size_t k=0; k<nodalConn[node].size(); k++) {
	    LO node2 = nodeMap[nodalConn[node][k]];
	    if (node2 != -1) {
	      for (LO m=nodeBegin[node2]; m<nodeBegin[node2+1]; m++) {
		columns[nnz++] = m;
	      }
	    }
	  }
	  rowBegin[j+1] = nnz;
	  // sort columns of row
	  //	  std::sort(columns.begin()+start, columns.begin()+nnz);
	}
      }
      for (LO j=0; j<numNode; j++) nodeMap[subNodes[sub][j]] = -1;
    }
    }
    // assemble sparse matrix
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
    std::vector<LO> nodeMap(m_numNode, -1);
    std::vector<LO> imap(maxNumDof, -1);
    std::vector<SX> rowValues(m_numDofElem);
#ifdef _OPENMP
#pragma omp for
#endif
    for (LO sub=0; sub<numSub; sub++) {
      LO numNode = subNodes[sub].size();
      for (LO j=0; j<numNode; j++) nodeMap[subNodes[sub][j]] = j;
      const std::vector<LO> & nodeBegin = subNodeBegin[sub];
      LO numElem = m_subElems[sub].size();
      const std::vector<LO> & rowBegin = subRowBegin[sub];
      const std::vector<LO> & columns = subColumns[sub];
      std::vector<SX> & values = subValues[sub];
      std::vector<SM> elemStiffMatrix, elemMassMatrix, elementLoadVector;
      double E, nu;
      int propNum;
      getMaterialProperties(sub, E, nu, propNum);
      for (LO i=0; i<numElem; i++) {
	LO elem = m_subElems[sub][i];
	getElementMatrices(m_elemConn[elem], E, nu, elemStiffMatrix, 
			   elemMassMatrix);
	determineElementLoadVector(elemMassMatrix, elementLoadVector);
	LO elementRow(0);
	for (size_t j=0; j<m_elemConn[elem].size(); j++) {
	  LO node = nodeMap[m_elemConn[elem][j]];
	  BDDC_TEST_FOR_EXCEPTION(node == -1, std::runtime_error, 
				  "invalid node number");
	  for (LO k=nodeBegin[node]; k<nodeBegin[node+1]; k++) {
	    for (LO m=rowBegin[k]; m<rowBegin[k+1]; m++) {
	      imap[columns[m]] = m;
	    }
	    int elementCol(0);
	    getElementMatrixRow
	      (elementRow, elemStiffMatrix, rowValues);
	    for (size_t n=0; n<m_elemConn[elem].size(); n++) {
	      LO node2 = nodeMap[m_elemConn[elem][n]];
	      BDDC_TEST_FOR_EXCEPTION(node2 == -1, std::runtime_error, 
				      "invalid node number");
	      for (LO p=nodeBegin[node2]; p<nodeBegin[node2+1]; p++) {
		LO col = p;
		BDDC_TEST_FOR_EXCEPTION(imap[col] == -1, std::runtime_error, 
				      "invalid imap[col] number");
		values[imap[col]] += rowValues[elementCol];
		elementCol++;
	      }
	    }
	    elementRow++;
	    for (LO m=rowBegin[k]; m<rowBegin[k+1]; m++) {
	      imap[columns[m]] = -1;
	    }
	  }
	}
      }
      for (LO j=0; j<numNode; j++) nodeMap[subNodes[sub][j]] = -1;
    }
    }
  }

  void addDiagonalStiffness(std::vector< std::vector<LO> > & subRowBegin, 
			    std::vector< std::vector<LO> > & subColumns,
			    std::vector< std::vector<SX> > & subValues,
			    const SM diagScaleFactor) const
  {
    if (diagScaleFactor == 1) return;
    LO numSub = subRowBegin.size();
    for (LO i=0; i<numSub; i++) {
      const std::vector<LO> & rowBegin = subRowBegin[i];
      const std::vector<LO> & columns = subColumns[i];
      std::vector<SX> & values = subValues[i];
      LO numRows = rowBegin.size()-1;
      for (LO j=0; j<numRows; j++) {
	for (LO k=rowBegin[j]; k<rowBegin[j+1]; k++) {
	  if (columns[k] == j) values[k] *= diagScaleFactor;
	}
      }
    }
  }

  void addAsymmetry(std::vector< std::vector<LO> > & subRowBegin, 
		    std::vector< std::vector<LO> > & subColumns,
		    std::vector< std::vector<SX> > & subValues) const
  {
    LO numSub = subRowBegin.size();
    for (LO i=0; i<numSub; i++) {
      const std::vector<LO> & rowBegin = subRowBegin[i];
      const std::vector<LO> & columns = subColumns[i];
      std::vector<SX> & values = subValues[i];
      LO numRows = rowBegin.size()-1;
      for (LO j=0; j<numRows; j++) {
	for (LO k=rowBegin[j]; k<rowBegin[j+1]; k++) {
	  if (columns[k] > j) values[k] *= 0.98;
	}
      }
    }
  }

 private:
  void outputModel()
  {
    const bool writeModel = m_Parameters->get("Output Model", false);
    if (writeModel == false) return;
    BDDC_TEST_FOR_EXCEPTION(m_numProc != 1, std::runtime_error, 
			    "numProc must be 1");
    std::ofstream fout;
    // coordinates
    fout.open("coords.dat");
    for (LO i=0; i<m_numNode; i++) {
      fout << std::setw(23) << std::setprecision(16);
      fout << m_x[i] << " " << m_y[i] << " " << m_z[i] << '\n';
    }
    fout.close();
    // element connectivity and blockIds
    const int numSub = m_subElems.size();
    fout.open("elemConnAndBlocks.dat");
    for (int i=0; i<numSub; i++) {
      double E, nu;
      int blockNum;
      getMaterialProperties(i, E, nu, blockNum);
      for (size_t j=0; j<m_subElems[i].size(); j++) {
	const int elem = m_subElems[i][j];
	for (size_t k=0; k<m_elemConn[elem].size(); k++) {
	  fout << m_elemConn[elem][k]+1 << " ";
	}
	fout << blockNum << "\n";
      }
    }
  }

  void getMaterialProperties(const LO sub, 
			     double & E, 
			     double & nu,
			     int & propNum)
  {
    const double E1 = m_Parameters->get("E1", 1.0);
    const double nu1 = m_Parameters->get("nu1", 0.3);
    const double E2 = m_Parameters->get("E2", 1.0);
    const double nu2 = m_Parameters->get("nu2", 0.3);
    const int matPropOption = m_Parameters->get("Material Property Option", 0);
    propNum = -1;
    if (m_readModel) {
      BDDC_TEST_FOR_EXCEPTION(matPropOption != 0, std::runtime_error, 
			      "matProOption must be 0");
      propNum = 1;
    }
    else {
      const LO ix = m_subIndices[3*sub+0];
      const LO iy = m_subIndices[3*sub+1];
      const LO iz = m_subIndices[3*sub+2];
      if (matPropOption == 0) {
	propNum = 1;
      }
      else if (matPropOption == 1) {
	if (ix%2 == 0) propNum = 1;
	else propNum = 2;
      }
      else if (matPropOption == 2) {
	BDDC_TEST_FOR_EXCEPTION(m_spatialDim < 2, std::runtime_error, 
			      "spatial dimension must be at least 2");
	if ((ix+iy)%2 == 0) propNum = 1;
	else propNum = 2;
      }
      else if (matPropOption == 3) {
	BDDC_TEST_FOR_EXCEPTION(m_spatialDim != 3, std::runtime_error, 
			      "spatial dimension must be 3");
	if ((ix+iy+iz)%2 == 0) propNum = 1;
	else propNum = 2;
      }
    }
    if (propNum == 1) {
      E = E1;
      nu = nu1;
    }
    else {
      E = E2;
      nu = nu2;
    }
    BDDC_TEST_FOR_EXCEPTION(propNum == -1, std::runtime_error, 
			    "invalid propNum");
  }

  void GenerateModel()
  {
    switch(m_spatialDim) {
    case 1:
      getModel1D(m_elemConn);
      break;
    case 2:
      getModel2D(m_elemConn);
      break;
    case 3:
      getModel3D(m_elemConn);
      break;
    default:
      break;
    }
    determineNodeBegLocalDof();
    LO numSubDir1PerProc = m_Parameters->get
      ("Number Subdomains Per Subregion Direction 1", 1);
    LO numSubDir2PerProc = m_Parameters->get
      ("Number Subdomains Per Subregion Direction 2", 1);
    LO numSubDir3PerProc = m_Parameters->get
      ("Number Subdomains Per Subregion Direction 3", 1);
    getSubDomainElements(numSubDir1PerProc, numSubDir2PerProc,
			 numSubDir3PerProc, m_subElems, m_subIndices);
    outputModel();
  }

  void ApplyEssentialBCs()
  {
    double xMin(0), xMax(0);
    if (m_x.size() > 0) {
      xMin = m_x[0];
      xMax = m_x[0];
    }
    for (LO i=1; i<m_numNode; i++) {
      if (m_x[i] < xMin) xMin = m_x[i];
      if (m_x[i] > xMax) xMax = m_x[i];
    }
    double xMinAll, xMaxAll;
    int err = MPI_Allreduce(&xMin, &xMinAll, 1, MPI_DOUBLE, MPI_MIN, m_Comm);
    err += MPI_Allreduce(&xMax, &xMaxAll, 1, MPI_DOUBLE, MPI_MAX, m_Comm);
    BDDC_TEST_FOR_EXCEPTION(err != 0, std::runtime_error, 
			    "MPI_Allreduce error");
    double length = xMaxAll - xMinAll;
    if (m_x.size() == 0) return;
    if (m_Parameters->get("Apply Left Side Essential BCs", false)) {
      ApplyEssentialBCs(xMinAll, length);
    }
    if (m_Parameters->get("Apply Right Side Essential BCs", false)) {
      ApplyEssentialBCs(xMaxAll, length);
    }
  }

  void ApplyEssentialBCs(const std::vector<LO> & constrainedNodes)
  {
    std::vector<bool> flagEBC(m_numNode, false);
    for (size_t i=0; i<constrainedNodes.size(); i++) {
      const LO node = constrainedNodes[i];
      flagEBC[node] = true;
    }
    LO numDofActive(0);
    std::vector<LO> nodeBeginNew(m_numNode+1, 0);
    for (LO i=0; i<m_numNode; i++) {
      for (LO j=m_nodeBegin[i]; j<m_nodeBegin[i+1]; j++) {
	if (flagEBC[i] == false) {
	  m_localDofs[numDofActive++] = m_localDofs[j];
	}
      }
      nodeBeginNew[i+1] = numDofActive;
    }
    m_nodeBegin = nodeBeginNew;
    m_localDofs.resize(numDofActive);
  }

  void ApplyEssentialBCs(const double xCoord,
			 const double length)
  {
    LO numDof = m_nodeBegin[m_numNode];
    std::vector<LO> dofMap(numDof);
    LO numDofActive(0);
    double tol(1e-10);
    std::vector<LO> nodeBeginNew(m_numNode+1, 0);
    for (LO i=0; i<m_numNode; i++) {
      bool flagEBC(false);
      if (fabs(m_x[i]-xCoord) < tol*length) flagEBC = true;
      for (LO j=m_nodeBegin[i]; j<m_nodeBegin[i+1]; j++) {
	if (flagEBC) dofMap[j] = -1;
	else {
	  dofMap[j] = numDofActive;
	  numDofActive++;
	}
      }
      nodeBeginNew[i+1] = numDofActive;
    }
    m_nodeBegin = nodeBeginNew;
  }

  void determineNodeBegLocalDof()
  {
    m_nodeBegin.resize(m_numNode+1, 0);
    int numDofPerNode = getNumDofPerNode();
    m_localDofs.resize(m_numNode*numDofPerNode);
    LO count(0);
    for (LO i=0; i<m_numNode; i++) {
      m_nodeBegin[i+1] = m_nodeBegin[i] + numDofPerNode;
      for (LO j=0; j<numDofPerNode; j++) {
	m_localDofs[count] = j;
	//	if (m_problemType == SCALARPDE) m_localDofs[count] = 6;
	if (m_problemType == SCALARPDE) m_localDofs[count] = 0;
	count++;
      }
    }
  }

  void shapeHex8(const double n1, 
		 const double n2, 
		 const double n3,
		 double* phi,
		 double* phi_1,
		 double* phi_2,
		 double* phi_3)
  {
    phi[0] = (1-n1)*(1-n2)*(1-n3)/8;
    phi[1] = (1+n1)*(1-n2)*(1-n3)/8;
    phi[2] = (1+n1)*(1+n2)*(1-n3)/8;
    phi[3] = (1-n1)*(1+n2)*(1-n3)/8;
    phi[4] = (1-n1)*(1-n2)*(1+n3)/8;
    phi[5] = (1+n1)*(1-n2)*(1+n3)/8;
    phi[6] = (1+n1)*(1+n2)*(1+n3)/8;
    phi[7] = (1-n1)*(1+n2)*(1+n3)/8;
    //  
    phi_1[0] = -(-1+n2)*(-1+n3)/8;
    phi_1[1] = (-1+n2)*(-1+n3)/8;
    phi_1[2] = -(1+n2)*(-1+n3)/8;
    phi_1[3] = (1+n2)*(-1+n3)/8;
    phi_1[4] = (-1+n2)*(1+n3)/8;
    phi_1[5] = -(-1+n2)*(1+n3)/8;
    phi_1[6] = (1+n2)*(1+n3)/8;
    phi_1[7] = -(1+n2)*(1+n3)/8;
    //  
    phi_2[0] = -(-1+n1)*(-1+n3)/8;
    phi_2[1] = (1+n1)*(-1+n3)/8;
    phi_2[2] = -(1+n1)*(-1+n3)/8;
    phi_2[3] = (-1+n1)*(-1+n3)/8;
    phi_2[4] = (-1+n1)*(1+n3)/8;
    phi_2[5] = -(1+n1)*(1+n3)/8;
    phi_2[6] = (1+n1)*(1+n3)/8;
    phi_2[7] = -(-1+n1)*(1+n3)/8;
    //  
    phi_3[0] = -(-1+n1)*(-1+n2)/8;
    phi_3[1] = (1+n1)*(-1+n2)/8;
    phi_3[2] = -(1+n1)*(1+n2)/8;
    phi_3[3] = (-1+n1)*(1+n2)/8;
    phi_3[4] = (-1+n1)*(-1+n2)/8;
    phi_3[5] = -(1+n1)*(-1+n2)/8;
    phi_3[6] = (1+n1)*(1+n2)/8;
    phi_3[7] = -(-1+n1)*(1+n2)/8;
  }

  void shapeQuad4(const double n1, 
		  const double n2, 
		  double* phi,
		  double* phi_1,
		  double* phi_2)
  {
    phi[0]=(1-n1)*(1-n2)/4;
    phi[1]=(1+n1)*(1-n2)/4;
    phi[2]=(1+n1)*(1+n2)/4;
    phi[3]=(1-n1)*(1+n2)/4;
    //
    phi_1[0]=-(1-n2)/4;
    phi_1[1]= (1-n2)/4;
    phi_1[2]= (1+n2)/4;
    phi_1[3]=-(1+n2)/4;
    //
    phi_2[0]=-(1-n1)/4;
    phi_2[1]=-(1+n1)/4;
    phi_2[2]= (1+n1)/4;
    phi_2[3]= (1-n1)/4;
  }

  void shapeBar2(const double n1, 
		 double* phi,
		 double* phi_1)
  {
    phi[0]=(1-n1)/2;
    phi[1]=(1+n1)/2;
    //
    phi_1[0]= -0.5;
    phi_1[1]=  0.5;
  }

  void calculateGradientsAndJacobian3D(const double* x, 
				       const double* y, 
				       const double* z, 
				       const double eta1,
				       const double eta2,
				       const double eta3,
				       double* phi,
				       double* phid, 
				       double & detJ)
  {
    double phi_1[8], phi_2[8], phi_3[8];
    shapeHex8(eta1, eta2, eta3, phi, phi_1, phi_2, phi_3);
    Teuchos::BLAS<int, double>  BLAS;
    int INCX(1), INCY(1);
    double x_1 = BLAS.DOT(8, phi_1, INCX, x, INCY);
    double y_1 = BLAS.DOT(8, phi_1, INCX, y, INCY);
    double z_1 = BLAS.DOT(8, phi_1, INCX, z, INCY);
    double x_2 = BLAS.DOT(8, phi_2, INCX, x, INCY);
    double y_2 = BLAS.DOT(8, phi_2, INCX, y, INCY);
    double z_2 = BLAS.DOT(8, phi_2, INCX, z, INCY);
    double x_3 = BLAS.DOT(8, phi_3, INCX, x, INCY);
    double y_3 = BLAS.DOT(8, phi_3, INCX, y, INCY);
    double z_3 = BLAS.DOT(8, phi_3, INCX, z, INCY);
    double J[] = {x_1, x_2, x_3, y_1, y_2, y_3, z_1, z_2, z_3};
    detJ = x_1*(y_2*z_3 - z_2*y_3) +
           y_1*(z_2*x_3 - x_2*z_3) +
           z_1*(x_2*y_3 - y_2*x_3);
    int IPIV[3], NRHS(8), INFO;
    for (int i=0; i<NRHS; i++) {
      phid[3*i+0] = phi_1[i];
      phid[3*i+1] = phi_2[i];
      phid[3*i+2] = phi_3[i];
    }
    Teuchos::LAPACK<int, double> LAPACK;
    LAPACK.GESV(3, NRHS, J, 3, IPIV, phid, 3, &INFO);
  }

  void calculateGradientsAndJacobian2D(const double* x, 
				       const double* y, 
				       const double eta1,
				       const double eta2,
				       double* phi,
				       double* phid, 
				       double & detJ)
  {
    double phi_1[4], phi_2[4];
    shapeQuad4(eta1, eta2, phi, phi_1, phi_2);
    Teuchos::BLAS<int, double>  BLAS;
    int INCX(1), INCY(1);
    double x_1 = BLAS.DOT(4, phi_1, INCX, x, INCY);
    double y_1 = BLAS.DOT(4, phi_1, INCX, y, INCY);
    double x_2 = BLAS.DOT(4, phi_2, INCX, x, INCY);
    double y_2 = BLAS.DOT(4, phi_2, INCX, y, INCY);
    double J[] = {x_1, x_2, y_1, y_2};
    detJ = x_1*y_2 - y_1*x_2;
    int IPIV[2], NRHS(4), INFO;
    for (int i=0; i<NRHS; i++) {
      phid[2*i+0] = phi_1[i];
      phid[2*i+1] = phi_2[i];
    }
    Teuchos::LAPACK<int, double> LAPACK;
    LAPACK.GESV(2, NRHS, J, 2, IPIV, phid, 2, &INFO);
  }

  void calculateGradientsAndJacobian1D(const double* x, 
				       const double eta1,
				       double* phi,
				       double* phid, 
				       double & detJ)
  {
    double phi_1[2];
    shapeBar2(eta1, phi, phi_1);
    Teuchos::BLAS<int, double>  BLAS;
    int INCX(1), INCY(1);
    double x_1 = BLAS.DOT(2, phi_1, INCX, x, INCY);
    double J = x_1;
    detJ = J;
    int IPIV[1], NRHS(2), INFO;
    for (int i=0; i<NRHS; i++) {
      phid[i+0] = phi_1[i];
    }
    Teuchos::LAPACK<int, double> LAPACK;
    LAPACK.GESV(1, NRHS, &J, 1, IPIV, phid, 1, &INFO);
  }

  void elemMatricesPoisson3D(const double* x, 
			     const double* y, 
			     const double* z, 
			     const double eta1, 
			     const double eta2, 
			     const double eta3, 
			     const double E,
			     const double rho,
			     double* dK, 
			     double* dM)
  {
    double phi[8], phid[24], detJ;
    calculateGradientsAndJacobian3D(x, y, z, eta1, eta2, eta3, phi, phid, detJ);
    double ALPHA = detJ*E; // diffusion coefficient is E
    double BETA = 0;
    Teuchos::BLAS<int, double>  BLAS;
    BLAS.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, 8, 8, 3, ALPHA, phid, 3, 
	      phid, 3, BETA, dK, 8);
    ALPHA = detJ*rho; // density is rho
    BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, 8, 8, 1, ALPHA, phi, 8, 
	      phi, 8, BETA, dM, 8);
  }

  void elemMatricesPoisson2D(const double* x, 
			     const double* y, 
			     const double eta1, 
			     const double eta2, 
			     const double E,
			     const double rho,
			     double* dK, 
			     double* dM)
  {
    double phi[4], phid[8], detJ;
    calculateGradientsAndJacobian2D(x, y, eta1, eta2, phi, phid, detJ);
    double ALPHA = detJ*E; // diffusion coefficient is E
    double BETA = 0;
    Teuchos::BLAS<int, double>  BLAS;
    BLAS.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, 4, 4, 2, ALPHA, phid, 2, 
	      phid, 2, BETA, dK, 4);
    ALPHA = detJ*rho; // density is rho
    BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, 4, 4, 1, ALPHA, phi, 4, 
	      phi, 4, BETA, dM, 4);
  }

  void elemMatricesPoisson1D(const double* x, 
			     const double eta1, 
			     const double E,
			     const double rho,
			     double* dK, 
			     double* dM)
  {
    double phi[2], phid[2], detJ;
    calculateGradientsAndJacobian1D(x, eta1, phi, phid, detJ);
    double ALPHA = detJ*E; // diffusion coefficient is E
    double BETA = 0;
    Teuchos::BLAS<int, double>  BLAS;
    BLAS.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, 2, 2, 1, ALPHA, phid, 1, 
	      phid, 1, BETA, dK, 2);
    ALPHA = detJ*rho; // density is rho
    BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, 2, 2, 1, ALPHA, phi, 2, 
	      phi, 2, BETA, dM, 2);
  }

  void calculateKmat(double E, 
		     double nu, 
		     double* Kmat)
  {
    memset(Kmat, 0, 36*sizeof(double));
    double G = E/2.0/(1+nu);
    double lambda = E*nu/(1+nu)/(1-2*nu);
    for (int i=0; i<3; i++) {
      Kmat[i+6*i] = 2*G;
      for (int j=0; j<3; j++) {
	Kmat[i+6*j] += lambda;
      }
    }
    for (int i=3; i<6; i++) {
      Kmat[i+6*i] = G;
    }
  }

  void calculateKmat2D(double E, 
		       double nu, 
		       double* Kmat)
  {
    const double sfac = E/(1-nu*nu);
    Kmat[0] = sfac;
    Kmat[1] = sfac*nu;
    Kmat[2] = 0;
    Kmat[3] = sfac*nu;
    Kmat[4] = sfac;
    Kmat[5] = 0;
    Kmat[6] = 0;
    Kmat[7] = 0;
    Kmat[8] = sfac*(1 - nu)/2;
  }

  void elemMatricesElasticity3D(const double* x, 
				const double* y, 
				const double* z, 
				const double eta1, 
				const double eta2, 
				const double eta3, 
				const double E,
				const double nu,
				const double rho,
				double* dK, 
				double* dM)
  {
    double phi[8], phid[24], detJ;
    calculateGradientsAndJacobian3D(x, y, z, eta1, eta2, eta3, phi, phid, detJ);
    int i1[8], i2[8], i3[8];
    for (int i=0; i<8; i++) {
      i1[i] = 3*i + 0;
      i2[i] = 3*i + 1;
      i3[i] = 3*i + 2;
    }
    double A[6*24];
    memset(A, 0, 6*24*sizeof(double));
    for (int i=0; i<8; i++) {
      double row1 = phid[0+3*i];
      double row2 = phid[1+3*i];
      double row3 = phid[2+3*i];
      A[0+6*i1[i]] = row1;
      A[1+6*i2[i]] = row2;
      A[2+6*i3[i]] = row3;
      A[3+6*i1[i]] = row2;
      A[3+6*i2[i]] = row1;
      A[4+6*i2[i]] = row3;
      A[4+6*i3[i]] = row2;
      A[5+6*i3[i]] = row1;
      A[5+6*i1[i]] = row3;
    }
    double Kmat[36];
    calculateKmat(E, nu, Kmat);
    Teuchos::BLAS<int, double>  BLAS;
    double ALPHA = detJ;
    double BETA = 0;
    double KmatA[6*24];
    BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 6, 24, 6, ALPHA, Kmat, 6, 
	      A, 6, BETA, KmatA, 6);
    ALPHA = 1.0;
    BLAS.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, 24, 24, 6, ALPHA, A, 6, 
	      KmatA, 6, BETA, dK, 24);
    for (int i=0; i<8; i++) {
      for (int j=0; j<8; j++) {
	double value = phi[i]*phi[j]*detJ;
	dM[i1[i]+24*i1[j]] = value;
	dM[i2[i]+24*i2[j]] = value;
	dM[i3[i]+24*i3[j]] = value;
      }
    }
  }

  void elemMatricesElasticity2D(const double* x, 
				const double* y, 
				const double eta1, 
				const double eta2, 
				const double E,
				const double nu,
				const double rho,
				double* dK, 
				double* dM)
  {
    double phi[4], phid[8], detJ;
    calculateGradientsAndJacobian2D(x, y, eta1, eta2, phi, phid, detJ);
    int i1[4], i2[4];
    for (int i=0; i<4; i++) {
      i1[i] = 2*i + 0;
      i2[i] = 2*i + 1;
    }
    double A[3*8];
    memset(A, 0, 3*8*sizeof(double));
    for (int i=0; i<4; i++) {
      double row1 = phid[0+2*i];
      double row2 = phid[1+2*i];
      A[0+3*i1[i]] = row1;
      A[1+3*i2[i]] = row2;
      A[2+3*i1[i]] = row2;
      A[2+3*i2[i]] = row1;
    }
    double Kmat[9];
    calculateKmat2D(E, nu, Kmat);
    Teuchos::BLAS<int, double>  BLAS;
    double ALPHA = detJ;
    double BETA = 0;
    double KmatA[3*8];
    BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 3, 8, 3, ALPHA, Kmat, 3, 
	      A, 3, BETA, KmatA, 3);
    ALPHA = 1.0;
    BLAS.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, 8, 8, 3, ALPHA, A, 3, 
	      KmatA, 3, BETA, dK, 8);
    for (int i=0; i<4; i++) {
      for (int j=0; j<4; j++) {
	double value = phi[i]*phi[j]*detJ;
	dM[i1[i]+8*i1[j]] = value;
	dM[i2[i]+8*i2[j]] = value;
      }
    }
  }

  void setQuadraturePoints()
  {
    if (m_Parameters->get("Quad Type 2 Point", 0) == 0) {
      // standard Gauss-Legendre integration points
      m_Aa_quad[0] = -0.577350269189626;
      m_Aa_quad[1] =  0.577350269189626;
    }
    else {
      // Gauss-Lobatto integration points
      m_Aa_quad[0] = -1.0;
      m_Aa_quad[1] =  1.0;
    }
  }

  void calculateElementMatrices3D(const std::vector<LO> & nodes,
				  const double E,
				  const double nu,
				  std::vector<double> & K, 
				  std::vector<double> & M)
  {
    int numNode = nodes.size();
    BDDC_TEST_FOR_EXCEPTION(numNode != 8, std::runtime_error, 
			    "numNode must be 8");
    double x[8], y[8], z[8];
    for (int i=0; i<numNode; i++) {
      x[i] = m_x[nodes[i]];
      y[i] = m_y[nodes[i]];
      z[i] = m_z[nodes[i]];
    }
    double Ha_quad[] = {1.0, 1.0};
    int numDof = 8;
    if (m_problemType == ELASTICITY) numDof = 24;
    std::vector<double> dK(numDof*numDof), dM(numDof*numDof);
    int numQuad = 2;
    const double rho(1.0);
    for (int i=0; i<numQuad; i++) {
      for (int j=0; j<numQuad; j++) {
	for (int k=0; k<numQuad; k++) {
	  double eta1 = m_Aa_quad[i];
	  double eta2 = m_Aa_quad[j];
	  double eta3 = m_Aa_quad[k];
	  if (m_problemType == SCALARPDE) {
	    elemMatricesPoisson3D(x, y, z, eta1, eta2, eta3, E, rho,
				  dK.data(), dM.data());
	  }
	  else {
	    elemMatricesElasticity3D(x, y, z, eta1, eta2, eta3, E, nu, rho,
				     dK.data(), dM.data());
	  }
	  for (int m=0; m<numDof*numDof; m++) {
	    K[m] += dK[m]*Ha_quad[i]*Ha_quad[j]*Ha_quad[k];
	    M[m] += dM[m]*Ha_quad[i]*Ha_quad[j]*Ha_quad[k];
	  }
	}
      }
    }
  }

  void calculateElementMatrices2D(const std::vector<LO> & nodes, 
				  const double E,
				  const double nu,
				  std::vector<double> & K, 
				  std::vector<double> & M)
  {
    int numNode = nodes.size();
    BDDC_TEST_FOR_EXCEPTION(numNode != 4, std::runtime_error, 
			    "numNode must be 4");
    double x[4], y[4];
    for (int i=0; i<numNode; i++) {
      x[i] = m_x[nodes[i]];
      y[i] = m_y[nodes[i]];
    }
    double Ha_quad[] = {1.0, 1.0};
    int numDof = 4;
    if (m_problemType == ELASTICITY) numDof = 8;
    std::vector<double> dK(numDof*numDof), dM(numDof*numDof);
    int numQuad = 2;
    const double rho(1.0);
    for (int i=0; i<numQuad; i++) {
      for (int j=0; j<numQuad; j++) {
	double eta1 = m_Aa_quad[i];
	double eta2 = m_Aa_quad[j];
	if (m_problemType == SCALARPDE) {
	  elemMatricesPoisson2D(x, y, eta1, eta2, E, rho,
				dK.data(), dM.data());
	}
	else {
	  elemMatricesElasticity2D(x, y, eta1, eta2, E, nu, rho,
				   dK.data(), dM.data());
	}
	for (int m=0; m<numDof*numDof; m++) {
	  K[m] += dK[m]*Ha_quad[i]*Ha_quad[j];
	  M[m] += dM[m]*Ha_quad[i]*Ha_quad[j];
	}
      }
    }
  }

  void calculateElementMatrices1D(const std::vector<LO> & nodes, 
				  const double E,
				  std::vector<double> & K, 
				  std::vector<double> & M)
  {
    int numNode = nodes.size();
    BDDC_TEST_FOR_EXCEPTION(numNode != 2, std::runtime_error, 
			    "numNode must be 2");
    double x[2];
    for (int i=0; i<numNode; i++) {
      x[i] = m_x[nodes[i]];
    }
    double Ha_quad[] = {1.0, 1.0};
    int numDof = 2;
    std::vector<double> dK(numDof*numDof), dM(numDof*numDof);
    int numQuad = 2;
    const double rho(1.0);
    for (int i=0; i<numQuad; i++) {
      double eta1 = m_Aa_quad[i];
      elemMatricesPoisson1D(x, eta1, E, rho, dK.data(), dM.data());
      for (int m=0; m<numDof*numDof; m++) {
	K[m] += dK[m]*Ha_quad[i];
	M[m] += dM[m]*Ha_quad[i];
      }
    }
  }

  void getElementMatrices(const std::vector<LO> & nodes,
			  const double E,
			  const double nu,
			  std::vector<SM> & elemStiffMatrix,
			  std::vector<SM> & elemMassMatrix)
  {
    int numDofPerNode = getNumDofPerNode();
    int numNodePerElem = power2(m_spatialDim);
    m_numDofElem = numNodePerElem*numDofPerNode;
    elemStiffMatrix.assign(m_numDofElem*m_numDofElem, 0);
    elemMassMatrix.assign(m_numDofElem*m_numDofElem, 0);
    if (m_spatialDim == 1) {
      calculateElementMatrices1D(nodes, E, elemStiffMatrix, elemMassMatrix);
    }
    else if (m_spatialDim == 2) {
      calculateElementMatrices2D(nodes, E, nu, elemStiffMatrix, elemMassMatrix);
    }
    else if (m_spatialDim == 3) {
      calculateElementMatrices3D(nodes, E, nu, elemStiffMatrix, elemMassMatrix);
    }
  }

  void addElasticity(int i, 
		     int j, 
		     std::vector<LO> & nodes,
		     int numDofElem,
		     double scaleFactorStiff,
		     std::vector<double> & elemMatrix) const
  {
    double dirCosine[3], delta[6], L;
    getDirectionCosines(nodes[i], nodes[j], dirCosine, L);
    int dofvec[6];
    for (int ldof=0; ldof<m_spatialDim; ldof++) {
      int ii = ldof+0*m_spatialDim;
      int jj = ldof+1*m_spatialDim;
      delta[ii] = -dirCosine[ldof];
      delta[jj] =  dirCosine[ldof];
      dofvec[ii] = ldof+i*m_spatialDim;
      dofvec[jj] = ldof+j*m_spatialDim;
    }
    for (int i=0; i<2*m_spatialDim; i++) {
      for (int j=0; j<2*m_spatialDim; j++) {
	int row = dofvec[i];
	int col = dofvec[j];
	elemMatrix[row+col*numDofElem] += scaleFactorStiff*delta[i]*delta[j];
      }
    }
  }

  void getDirectionCosines(int nodei, 
			   int nodej, 
			   double* dirCosine,
			   double & L) const
  {
    double dx = m_x[nodej] - m_x[nodei];
    double dy = m_y[nodej] - m_y[nodei];
    double dz = 0;
    if (m_spatialDim == 3) dz = m_z[nodej] - m_z[nodei];
    L = sqrt(dx*dx + dy*dy + dz*dz);
    dirCosine[0] = dx/L;
    dirCosine[1] = dy/L;
    dirCosine[2] = dz/L;
  }

  int estimateMaxEntriesPerRow() const
  {
    int numDofPerNode = getNumDofPerNode();
    int maxEntriesPerRow = numDofPerNode*myIpow(3, m_spatialDim);
    return maxEntriesPerRow;
  }

  void getColumns(std::vector<LO> & nodes, 
		  std::vector<LO> & columns) const
  {
    int count(0);
    for (size_t i=0; i<nodes.size(); i++) {
      LO node = nodes[i];
      for (LO j=m_nodeBegin[node]; j<m_nodeBegin[node+1]; j++) {
	columns[count++] = j;
      }
    }
  }

  void getElementMatrixRow(const int localRow, 
			   const std::vector<SM> & elemStiffMatrix,
			   std::vector<SX> & rowValues) const
  {
    SX a(0);
    int index;
    rowValues.resize(m_numDofElem);
    switch(m_analysisType) {
    case STANDARD:
      for (LO i=0; i<m_numDofElem; i++) {
	index = localRow+m_numDofElem*i;
	rowValues[i] = getNumber(elemStiffMatrix[index], 0, a);
      }
      break;
      /*
    case HELMHOLTZA:
      for (LO i=0; i<m_numDofElem; i++) {
	index = localRow+m_numDofElem*i;
	rowValues[i] =  
	  getNumber(elemStiffMatrix[index], 
		    betaDamping*elemStiffMatrix[index]*omega, a);
	rowValues[i] -= getNumber(elemMassMatrix[index]*omega*omega, 0, a);
      }
      break;
      */
    default:
      break;
    }
  }

  double getNumber(double realPart,
		   double imagPart,
		   double a) const
  {
    return realPart;
  }

  std::complex<double> getNumber(double realPart,
				 double imagPart,
				 std::complex<double> a) const
  {
    return std::complex<double>(realPart, imagPart);
  }

  void getModel1D(std::vector< std::vector<LO> > & elemConn)
  {
    LO numSubDir1 = m_Parameters->get("Number of Subdomains Direction 1", 1);
    double lengthDir1 = m_Parameters->get("Length Direction 1", double(1));
    LO numElemPerSubDir1 = 
      m_Parameters->get("Number of Elements Per Subregion Direction 1", 1);
    BDDC_TEST_FOR_EXCEPTION(m_numProc < numSubDir1, std::runtime_error, 
			    "not enough processors");
    LO numActiveProc = numSubDir1;
    if (m_myPID > numActiveProc-1) {
      m_numNode = 0;
      m_numElem = 0;
      return;
    }
    m_numNode = numElemPerSubDir1+1;
    m_numElem = numElemPerSubDir1;
    double lengthSubDir1 = lengthDir1/numSubDir1;
    const double hElem = lengthSubDir1/numElemPerSubDir1;
    m_x.resize(m_numNode, 0);
    m_y.resize(m_numNode, 0);
    m_z.resize(m_numNode, 0);
    m_nodeGIDs.resize(m_numNode);
    m_x[0] = lengthSubDir1*m_myPID;
    m_nodeGIDs[0] = numElemPerSubDir1*m_myPID;
    for (LO i=0; i<numElemPerSubDir1; i++) {
      m_x[i+1] = m_x[i] + hElem;
      m_nodeGIDs[i+1] = m_nodeGIDs[i] + 1;
    }
    elemConn.resize(m_numElem);
    for (LO i=0; i<numElemPerSubDir1; i++) {
      elemConn[i].resize(2);
      elemConn[i][0] = i;
      elemConn[i][1] = i+1;
    }
  }
  
  void getModel2D(std::vector< std::vector<LO> > & elemConn)
  {
    LO numSubDir1 = m_Parameters->get("Number of Subdomains Direction 1", 1);
    LO numSubDir2 = m_Parameters->get("Number of Subdomains Direction 2", 1);
    double lengthDir1 = m_Parameters->get("Length Direction 1", double(1));
    double lengthDir2 = m_Parameters->get("Length Direction 2", double(1));
    LO numElemPerSubDir1 = 
      m_Parameters->get("Number of Elements Per Subregion Direction 1", 1);
    LO numElemPerSubDir2 = 
      m_Parameters->get("Number of Elements Per Subregion Direction 2", 1);
    BDDC_TEST_FOR_EXCEPTION(m_numProc < numSubDir1*numSubDir2, 
			    std::runtime_error, "not enough processors");
    LO numActiveProc = numSubDir1*numSubDir2;
    if (m_myPID > numActiveProc-1) {
      m_numNode = 0;
      m_numElem = 0;
      return;
    }
    LO nnodeDir1 = numElemPerSubDir1+1;
    LO nnodeDir2 = numElemPerSubDir2+1;
    m_numNode = nnodeDir1*nnodeDir2;
    m_numElem = numElemPerSubDir1*numElemPerSubDir2;
    double lengthSubDir1 = lengthDir1/numSubDir1;
    double lengthSubDir2 = lengthDir2/numSubDir2;
    double hElem1 = lengthSubDir1/numElemPerSubDir1;
    double hElem2 = lengthSubDir2/numElemPerSubDir2;
    int i1 = m_myPID%numSubDir1;
    int i2 = m_myPID/numSubDir1;
    GO nnodeDir1All = numElemPerSubDir1*numSubDir1+1;
    LO nodeGIDll = i1*numElemPerSubDir1 + i2*nnodeDir1All*numElemPerSubDir2;
    SM xLL = i1*lengthSubDir1;
    SM yLL = i2*lengthSubDir2;
    m_x.resize(m_numNode, 0);
    m_y.resize(m_numNode, 0);
    m_z.resize(m_numNode, 0);
    m_nodeGIDs.resize(m_numNode);
    for (LO i=0; i<nnodeDir1; i++) {
      for (LO j=0; j<nnodeDir2; j++) {
	LO index = i+j*nnodeDir1;
	m_x[index] = xLL + hElem1*i;
	m_y[index] = yLL + hElem2*j;
	m_nodeGIDs[index] = nodeGIDll + i + nnodeDir1All*j;
      }
    }
    elemConn.resize(m_numElem);
    for (LO i=0; i<numElemPerSubDir1; i++) {
      for (LO j=0; j<numElemPerSubDir2; j++) {
	LO index = i+j*numElemPerSubDir1;
	elemConn[index].resize(4);
	elemConn[index][0] = i + j*nnodeDir1;
	elemConn[index][1] = elemConn[index][0] + 1;
	elemConn[index][3] = elemConn[index][0] + nnodeDir1;
	elemConn[index][2] = elemConn[index][3] + 1;
      }
    }
  }

  void getModel3D(std::vector< std::vector<LO> > & elemConn)
  {
    LO numSubDir1 = m_Parameters->get("Number of Subdomains Direction 1", 1);
    LO numSubDir2 = m_Parameters->get("Number of Subdomains Direction 2", 1);
    LO numSubDir3 = m_Parameters->get("Number of Subdomains Direction 3", 1);
    double lengthDir1 = m_Parameters->get("Length Direction 1", double(1));
    double lengthDir2 = m_Parameters->get("Length Direction 2", double(1));
    double lengthDir3 = m_Parameters->get("Length Direction 3", double(1));
    LO numElemPerSubDir1 = 
      m_Parameters->get("Number of Elements Per Subregion Direction 1", 1);
    LO numElemPerSubDir2 = 
      m_Parameters->get("Number of Elements Per Subregion Direction 2", 1);
    LO numElemPerSubDir3 = 
      m_Parameters->get("Number of Elements Per Subregion Direction 3", 1);
    BDDC_TEST_FOR_EXCEPTION(m_numProc < numSubDir1*numSubDir2*numSubDir3, 
			    std::runtime_error, "not enough processors");
    LO numActiveProc = numSubDir1*numSubDir2*numSubDir3;
    if (m_myPID > numActiveProc-1) {
      m_numNode = 0;
      m_numElem = 0;
      return;
    }
    LO nnodeDir1 = numElemPerSubDir1+1;
    LO nnodeDir2 = numElemPerSubDir2+1;
    LO nnodeDir3 = numElemPerSubDir3+1;
    m_numNode = nnodeDir1*nnodeDir2*nnodeDir3;
    m_numElem = numElemPerSubDir1*numElemPerSubDir2*numElemPerSubDir3;
    double lengthSubDir1 = lengthDir1/numSubDir1;
    double lengthSubDir2 = lengthDir2/numSubDir2;
    double lengthSubDir3 = lengthDir3/numSubDir3;
    double hElem1 = lengthSubDir1/numElemPerSubDir1;
    double hElem2 = lengthSubDir2/numElemPerSubDir2;
    double hElem3 = lengthSubDir3/numElemPerSubDir3;
    int i1 = m_myPID%numSubDir1;
    int i3 = m_myPID/(numSubDir1*numSubDir2);
    int i2 = (m_myPID-i3*numSubDir1*numSubDir2)/numSubDir1;
    GO nnodeDir1All = numElemPerSubDir1*numSubDir1+1;
    GO nnodeDir2All = numElemPerSubDir2*numSubDir2+1;
    LO nodeGIDll = i1*numElemPerSubDir1 + 
      i2*numElemPerSubDir2*nnodeDir1All +
      i3*numElemPerSubDir3*nnodeDir1All*nnodeDir2All;
    SM xLL = i1*lengthSubDir1;
    SM yLL = i2*lengthSubDir2;
    SM zLL = i3*lengthSubDir3;
    m_x.resize(m_numNode, 0);
    m_y.resize(m_numNode, 0);
    m_z.resize(m_numNode, 0);
    m_nodeGIDs.resize(m_numNode);
    GO nnode12 = nnodeDir1All*nnodeDir2All;
    for (LO i=0; i<nnodeDir1; i++) {
      for (LO j=0; j<nnodeDir2; j++) {
	for (LO k=0; k<nnodeDir3; k++) {
	  LO index = i+j*nnodeDir1+k*nnodeDir1*nnodeDir2;
	  m_x[index] = xLL + hElem1*i;
	  m_y[index] = yLL + hElem2*j;
	  m_z[index] = zLL + hElem3*k;
	  m_nodeGIDs[index] = nodeGIDll + i + nnodeDir1All*j + nnode12*k;
	}
      }
    }
    elemConn.resize(m_numElem);
    LO elem12 = numElemPerSubDir1*numElemPerSubDir2;
    for (LO i=0; i<numElemPerSubDir1; i++) {
      for (LO j=0; j<numElemPerSubDir2; j++) {
	for (LO k=0; k<numElemPerSubDir3; k++) {
	  LO index = i + j*numElemPerSubDir1 + k*elem12;
	  elemConn[index].resize(8);
	  elemConn[index][0] = i + j*nnodeDir1 + k*nnodeDir1*nnodeDir2 ;
	  elemConn[index][1] = elemConn[index][0] + 1;
	  elemConn[index][3] = elemConn[index][0] + nnodeDir1;
	  elemConn[index][2] = elemConn[index][3] + 1;
	  elemConn[index][4] = elemConn[index][0] + nnodeDir1*nnodeDir2;
	  elemConn[index][5] = elemConn[index][1] + nnodeDir1*nnodeDir2;
	  elemConn[index][6] = elemConn[index][2] + nnodeDir1*nnodeDir2;
	  elemConn[index][7] = elemConn[index][3] + nnodeDir1*nnodeDir2;
	}
      }
    }
  }

  void determineNodalConnectivity
    (const std::vector< std::vector<LO> > & elemConn,
     std::vector< std::vector<LO> > & nodalConn)
  {
    std::vector< std::vector<LO> > nodeElements(m_numNode);
    LO numElem = elemConn.size();
    for (LO i=0; i<numElem; i++) {
      for (size_t j=0; j<elemConn[i].size(); j++) {
	LO node = elemConn[i][j];
	nodeElements[node].push_back(i);
      }
    }
    nodalConn.resize(m_numNode);
    std::vector<LO> anode(m_numNode);
    std::vector<bool> nodeFlag(m_numNode, false);
    for (LO i=0; i<m_numNode; i++) {
      LO nanode = 0;
      for (size_t j=0; j<nodeElements[i].size(); j++) {
	LO elem = nodeElements[i][j];
	for (size_t k=0; k<elemConn[elem].size(); k++) {
	  LO node = elemConn[elem][k];
	  if (nodeFlag[node] == false) {
	    nodalConn[i].push_back(node);
	    anode[nanode] = node;
	    nodeFlag[node] = true;
	    nanode++;
	  }
	}
      }
      for (LO j=0; j<nanode; j++) {
	nodeFlag[anode[j]] = false;
      }
    }
  }
  
  void printElementMatrix(const double* matrix,
			  const char* fname) const
  {
    std::ofstream fout;
    fout.open(fname);
    for (int i=0; i<m_numDofElem; i++) {
      for (int j=0; j<m_numDofElem; j++) {
	fout << i+1 << " " << j+1 << std::setw(22) << std::setprecision(15)
	     << matrix[i+j*m_numDofElem] << std::endl;
      }
    }
    fout.close();
  }

  void getParts(const double value,
		double & realPart,
		double & imagPart) const
  {
    realPart = value;
    imagPart = 0;
  }

  void getParts(const std::complex<double> value,
		double & realPart,
		double & imagPart) const
  {
    realPart = std::real(value);
    imagPart = std::imag(value);
  }

  SM squared(const double value) const
  {
    return value*value;
  }

  SM squared(const std::complex<double> value) const
  {
    return value*std::conj(value);
  }

  void determineElementLoadVector(const std::vector<SM> & elemMassMatrix,
				  std::vector<SX> & elementLoadVector) const
  {
    int numDofPerNode = getNumDofPerNode();
    int numNodePerElem = power2(m_spatialDim);
    BDDC_TEST_FOR_EXCEPTION(m_numDofElem != numNodePerElem*numDofPerNode, 
			    std::runtime_error, "not enough processors");
    elementLoadVector.assign(m_numDofElem, 0);
    int loadDirection = m_Parameters->get("Load Direction", 1);
    for (int i=0; i<numNodePerElem; i++) {
      for (int j=0; j<numDofPerNode; j++) {
	if (j == loadDirection) {
	  LO index = i*numDofPerNode+j;
	  SX rowSum(0);
	  for (int k=0; k<m_numDofElem; k++) {
	    LO indexElem = index+k*m_numDofElem;
	    rowSum += elemMassMatrix[indexElem];
	  }
	  elementLoadVector[index] = rowSum;
	}
      }
    }
  }

  RCP<Teuchos::ParameterList> m_Parameters;
  MPI_Comm m_Comm;
  enum ProblemType m_problemType;
  enum AnalysisType m_analysisType;
  LO m_spatialDim;
  int m_myPID{-1}, m_numProc{-1}, m_numDofElem{0};
  LO m_numNode{0}, m_numElem{0};
  std::vector<double> m_x, m_y, m_z;
  std::vector<LO> m_nodeBegin, m_localDofs, m_subIndices;
  std::vector<GO> m_nodeGIDs;
  std::vector<SX> m_loadVector;
  std::vector< std::vector<LO> > m_elemConn;
  std::vector< std::vector<LO> > m_subElems;
  bool m_readModel{false};
#ifdef TENSORPRODUCTBDDC
  TensorProduct* m_TP{nullptr};
#endif
  double m_Aa_quad[2];

};

} // namespace bddc

#endif // BDDC_PROBLEMMAKER_H
