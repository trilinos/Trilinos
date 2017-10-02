
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

#ifndef PROBLEMMAKERBDDC_H
#define PROBLEMMAKERBDDC_H
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
#include <assert.h>

#include "metis.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

#include "shylu_enumsBDDC.h"

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
    m_spatialDim(Parameters->get("Spatial Dimension", 1))
  {
    MPI_Comm_size(Comm, &m_numProc);
    MPI_Comm_rank(Comm, &m_myPID);
    GenerateModel();
    ApplyEssentialBCs();
    GenerateConstraintEquations();
  }

  ~ProblemMaker()
  {
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
      assert (m_spatialDim != 1);
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
    return &m_nodeBegin[0];
  }

  const LO* getLocalDofs() const {
    return &m_localDofs[0];
  }

  const SM* getXcoord() const {
    return &m_x[0];
  }

  const SM* getYcoord() const {
    return &m_y[0];
  }

  const SM* getZcoord() const {
    return &m_z[0];
  }

  const GO* getNodeGlobalIDs() const {
    return &m_nodeGIDs[0];
  }

  const LO* getRowBegin() const {
    return &m_rowBegin[0];
  }

  const LO* getColumns() const {
    return &m_columns[0];
  }

  const SX* getValues() const {
    return &m_values[0];
  }

  const SX* getLoadVector() const {
    return &m_loadVector[0];
  }

  void printElementMatrices() const
  {
    printElementMatrix(&m_elemStiffMatrix[0], "stiffMatrix.dat");
    printElementMatrix(&m_elemMassMatrix[0], "massMatrix.dat");
  }

  void printSubdomainMatrix(const char* fnameBase,
			    int ID) const
  {
    char fname[101];
    sprintf(fname, "%s%d.dat", fnameBase, ID);
    LO i, j;
    std::ofstream fout;
    fout.open(fname);
    LO numDof = m_nodeBegin[m_numNode];
    for (LO i=0; i<numDof; i++) {
      for (LO j=m_rowBegin[i]; j<m_rowBegin[i+1]; j++) {
	fout << i+1 << " ";
	fout << m_columns[j]+1 << " ";
	SM realPart, imagPart;
	getParts(m_values[j], realPart, imagPart);
	fout << std::setw(22) << std::setprecision(15);
	fout << realPart;
	if (fabs(imagPart) != 0) {
	  fout << " " << imagPart;
	}
	fout << std::endl;
      }
    }
    fout.close();
  }

  void getConstrainEquations(LO & numRowsCT,
			     LO & numColsCT,
			     const LO* rowsCT,
			     const LO* colsCT,
			     const LO* rowBeginCT,
			     const LO* columnsCT,
			     const SM* matrixCT) const
  {
    numRowsCT = m_rowsCT.size();
    numColsCT = m_colsCT.size();
    rowsCT = &m_rowsCT[0];
    colsCT = &m_colsCT[0];
    rowBeginCT = &m_rowBeginCT[0];
    columnsCT = &m_columnsCT[0];
    matrixCT = &m_matrixCT[0];
  }

  void matrixVectorProduct(const SX* X, 
			   int numVectors, 
			   SX* AX) const
  {
    LO numDof = m_nodeBegin[m_numNode];
    for (int j=0; j<numVectors; j++) {
      const SX* x = &X[j*numDof];
      SX* Ax = &AX[j*numDof];
      for (LO i=0; i<numDof; i++) {
	Ax[i] = 0;
	for (LO j=m_rowBegin[i]; j<m_rowBegin[i+1]; j++) {
	  Ax[i] += m_values[j]*x[m_columns[j]];
	}
      }
    }
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
			    std::vector< std::vector<LO> > & subElems) const
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
    subElems.resize(numSub);
    for (LO i=0; i<numElems; i++) {
      LO ix = (xElem[i] - minX)/Hx;
      LO iy = (yElem[i] - minY)/Hy;
      LO iz = (zElem[i] - minZ)/Hz;
      if (m_spatialDim == 1) {
	LO sub = ix;
	subElems[sub].push_back(i);
      }
      else if (m_spatialDim == 2) {
	LO sub = ix + numSubDir1PerProc*iy;
	subElems[sub].push_back(i);
      }
      else if (m_spatialDim == 3) {
	LO sub = ix + numSubDir1PerProc*iy +
	  numSubDir1PerProc*numSubDir2PerProc*iz;
	subElems[sub].push_back(i);
      }
    }
  }

  void getSubDomainElements(LO numPartsPerSub,
			    std::vector< std::vector<LO> > & subElems) const
  {
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
    int* xadj   = C2.data();
    int* adjncy = C1.data();
    int nparts = numPartsPerSub;
    idx_t* part = parts.data();
    int edgecut, ncon(1);
    std::vector<idx_t> options(METIS_NOPTIONS, 0);
    METIS_SetDefaultOptions(&options[0]);
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
  }

  void getSubDomainNodeData(std::vector< std::vector<LO> > & subElems, 
			    std::vector< std::vector<LO> > & subNodes, 
			    std::vector< std::vector<LO> > & subNodeBegin,
			    std::vector< std::vector<LO> > & subLocalDofs) const
  {
    LO numSub = subElems.size();
    subNodes.resize(numSub);
    subNodeBegin.resize(numSub);
    subLocalDofs.resize(numSub);
    std::vector<LO> nodeFlag(m_numNode, -1), activeNodes(m_numNode);
    for (LO i=0; i<numSub; i++) {
      LO numElem = subElems[i].size();
      LO numNode(0), numDof(0);
      for (LO j=0; j<numElem; j++) {
	LO elem = subElems[i][j];
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
      subNodes[i].resize(numNode);
      subNodeBegin[i].resize(numNode+1, 0);
      subLocalDofs[i].resize(numDof);
      for (LO j=0; j<numNode; j++) nodeFlag[activeNodes[j]] = -1;
      numNode = numDof = 0;
      for (LO j=0; j<numElem; j++) {
	LO elem = subElems[i][j];
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

  void getSubdomainMatrices(std::vector< std::vector<LO> > & subElems, 
			    std::vector< std::vector<LO> > & subNodes,
			    std::vector< std::vector<LO> > & subNodeBegin,
			    std::vector< std::vector<LO> > & subRowBegin, 
			    std::vector< std::vector<LO> > & subColumns,
			    std::vector< std::vector<SX> > & subValues)
  {
    LO numSub = subElems.size();
    subRowBegin.resize(numSub);
    subColumns.resize(numSub);
    subValues.resize(numSub);
    std::vector< std::vector<LO> > nodalConn;
    determineNodalConnectivity(m_elemConn, nodalConn);
    double omega = m_Parameters->get("omega", 0.0);
    double betaDamping = m_Parameters->get("betaDamping", 0.0);
    std::vector<LO> nodeMap(m_numNode, -1);
    for (LO sub=0; sub<numSub; sub++) {
      LO numNode = subNodes[sub].size();
      for (LO j=0; j<numNode; j++) nodeMap[subNodes[sub][j]] = j;
      std::vector<LO> & nodeBegin = subNodeBegin[sub];
      std::vector<LO> & rowBegin = subRowBegin[sub];
      std::vector<LO> & columns = subColumns[sub];
      std::vector<SX> & values = subValues[sub];
      LO numDof = nodeBegin[numNode];
      rowBegin.resize(numDof+1, 0);
      columns.resize(0);
      LO nnz(0);
      for (LO i=0; i<numNode; i++) {
	LO node = subNodes[sub][i];
	for (LO j=nodeBegin[i]; j<nodeBegin[i+1]; j++) {
	  for (size_t k=0; k<nodalConn[node].size(); k++) {
	    LO node2 = nodeMap[nodalConn[node][k]];
	    if (node2 != -1) {
	      for (LO m=nodeBegin[node2]; m<nodeBegin[node2+1]; m++) {
		columns.push_back(m);
		nnz++;
	      }
	    }
	  }
	  rowBegin[j+1] = nnz;
	}
      }
      std::vector<LO> imap(numDof, -1);
      values.resize(nnz, 0);
      std::vector<SX> rowValues(m_numDofElem);
      LO numElem = subElems[sub].size();
      for (LO i=0; i<numElem; i++) {
	LO elementRow(0);
	LO elem = subElems[sub][i];
	for (size_t j=0; j<m_elemConn[elem].size(); j++) {
	  LO node = nodeMap[m_elemConn[elem][j]];
	  assert (node != -1);
	  for (LO k=nodeBegin[node]; k<nodeBegin[node+1]; k++) {
	    for (LO m=rowBegin[k]; m<rowBegin[k+1]; m++) {
	      imap[columns[m]] = m;
	    }
	    int elementCol(0);
	    getElementMatrixRow(elementRow, omega, betaDamping, rowValues);
	    for (size_t n=0; n<m_elemConn[elem].size(); n++) {
	      LO node2 = nodeMap[m_elemConn[elem][n]];
	      assert (node2 != -1);
	      for (LO p=nodeBegin[node2]; p<nodeBegin[node2+1]; p++) {
		LO col = p;
		assert (imap[col] != -1);
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

  void addDiagonalStiffness(std::vector< std::vector<LO> > & subRowBegin, 
			    std::vector< std::vector<LO> > & subColumns,
			    std::vector< std::vector<SX> > & subValues,
			    SM diagScaleFactor) const
  {
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
  void GenerateModel()
  {
    switch(m_spatialDim) {
    case 1:
      getModel1D(m_hElem, m_elemConn);
      break;
    case 2:
      getModel2D(m_hElem, m_elemConn);
      break;
    case 3:
      getModel3D(m_hElem, m_elemConn);
      break;
    default:
      break;
    }
    assembleSubdomainMatrix(m_hElem, m_elemConn);
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
    assert (err == 0);
    double length = xMaxAll - xMinAll;
    if (m_x.size() == 0) return;
    if (m_Parameters->get("Apply Left Side Essential BCs", false)) {
      ApplyEssentialBCs(xMinAll, length);
    }
    if (m_Parameters->get("Apply Right Side Essential BCs", false)) {
      ApplyEssentialBCs(xMaxAll, length);
    }
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
	  m_loadVector[numDofActive] = m_loadVector[j];
	  numDofActive++;
	}
      }
      nodeBeginNew[i+1] = numDofActive;
    }
    m_nodeBegin = nodeBeginNew;
    m_loadVector.resize(numDofActive);
    adjustProblemMatrix(numDof, numDofActive, dofMap);
  }

  void adjustProblemMatrix(const int numDof, 
			   int numDofActive, 
			   const std::vector<LO> & dofMap)
  {
    if (numDofActive == numDof) return;
    std::vector<LO> rowBeginNew(numDofActive+1, 0);
    numDofActive = 0;
    LO numTerm(0);
    for (LO i=0; i<numDof; i++) {
      if (dofMap[i] != -1) {
	for (LO j=m_rowBegin[i]; j<m_rowBegin[i+1]; j++) {
	  LO col = dofMap[m_columns[j]];
	  if (col != -1) {
	    m_columns[numTerm] = col;
	    m_values[numTerm] = m_values[j];
	    numTerm++;
	  }
	}
	numDofActive++;
	rowBeginNew[numDofActive] = numTerm;
      }
    }
    m_rowBegin = rowBeginNew;
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
	if (m_problemType == SCALARPDE) m_localDofs[count] = 6;
	count++;
      }
    }
  }

  void shapeHex8(double n1, 
		 double n2, 
		 double n3,
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

  void calculateGradientsAndJacobian(const double* x, 
				     const double* y, 
				     const double* z, 
				     double eta1,
				     double eta2,
				     double eta3,
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
    int IPIV[8], NRHS(8), INFO;
    for (int i=0; i<NRHS; i++) {
      phid[3*i+0] = phi_1[i];
      phid[3*i+1] = phi_2[i];
      phid[3*i+2] = phi_3[i];
    }
    Teuchos::LAPACK<int, double> LAPACK;
    LAPACK.GESV(3, NRHS, J, 3, IPIV, phid, 3, &INFO);
  }

  void elemMatricesPoisson(double* x, 
			   double* y, 
			   double* z, 
			   double eta1, 
			   double eta2, 
			   double eta3, 
			   double E,
			   double rho,
			   double* dK, 
			   double* dM)
  {
    double phi[8], phid[24], detJ;
    calculateGradientsAndJacobian(x, y, z, eta1, eta2, eta3, phi, phid, detJ);
    double ALPHA = detJ*E; // diffusion coefficient is E
    double BETA = 0;
    Teuchos::BLAS<int, double>  BLAS;
    BLAS.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, 8, 8, 3, ALPHA, phid, 3, 
	      phid, 3, BETA, dK, 8);
    ALPHA = detJ*rho; // density is rho
    BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::TRANS, 8, 8, 1, ALPHA, phi, 8, 
	      phi, 8, BETA, dM, 8);
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

  void elemMatricesElasticity(double* x, 
			      double* y, 
			      double* z, 
			      double eta1, 
			      double eta2, 
			      double eta3, 
			      double E,
			      double nu,
			      double rho,
			      double* dK, 
			      double* dM)
  {
    double phi[8], phid[24], detJ;
    calculateGradientsAndJacobian(x, y, z, eta1, eta2, eta3, phi, phid, detJ);
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

  void calculateElementMatrices(std::vector<LO> & nodes, 
				std::vector<double> & K, 
				std::vector<double> & M)
  {
    int numNode = nodes.size();
    assert (numNode == 8);
    double x[8], y[8], z[8];
    for (int i=0; i<numNode; i++) {
      x[i] = m_x[nodes[i]];
      y[i] = m_y[nodes[i]];
      z[i] = m_z[nodes[i]];
    }
    double Aa_quad[] = {-0.577350269189626, 0.577350269189626};
    double Ha_quad[] = {1.0, 1.0};
    int numDof = 8;
    if (m_problemType == ELASTICITY) numDof = 24;
    std::vector<double> dK(numDof*numDof), dM(numDof*numDof);
    int numQuad = 2;
    for (int i=0; i<numQuad; i++) {
      for (int j=0; j<numQuad; j++) {
	for (int k=0; k<numQuad; k++) {
	  double eta1 = Aa_quad[i];
	  double eta2 = Aa_quad[j];
	  double eta3 = Aa_quad[k];
	  if (m_problemType == SCALARPDE) {
	    double E(1.0), rho(1.0);
	    elemMatricesPoisson(x, y, z, eta1, eta2, eta3, E, rho,
				dK.data(), dM.data());
	  }
	  else {
	    double E(1.0), nu(0.3), rho(1.0);
	    elemMatricesElasticity(x, y, z, eta1, eta2, eta3, E, nu, rho,
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

  void getElementMatrices(double hElem,
			  std::vector<LO> & nodes,
			  std::vector<double> & elemStiffMatrix,
			  std::vector<double> & elemMassMatrix)
  {
    int numDofPerNode = getNumDofPerNode();
    int numNodePerElem = power2( m_spatialDim);
    m_numDofElem = numNodePerElem*numDofPerNode;
    elemStiffMatrix.resize(m_numDofElem*m_numDofElem, 0);
    elemMassMatrix.resize(m_numDofElem*m_numDofElem, 0);
    // calculate actual finite element matrices for 3D problems
    if (m_spatialDim == 3) {
      calculateElementMatrices(nodes, elemStiffMatrix, elemMassMatrix);
      return;
    }
    double scaleFactorMass = pow(hElem, m_spatialDim);
    double artificialStiffness = 
      m_Parameters->get("Artificial Foundation Stiffness", double(0));
    for (int i=0; i<m_numDofElem; i++) {
      int index = i+m_numDofElem*i;
      elemMassMatrix[index] = scaleFactorMass/numNodePerElem;
      elemStiffMatrix[index] = artificialStiffness*elemMassMatrix[index];
    }
    double scaleFactorStiff = pow(hElem, m_spatialDim-2);
    switch (m_problemType) {
    case SCALARPDE:
      for (int i=0; i<m_numDofElem; i++) {
	for (int j=i+1; j<m_numDofElem; j++) {
	  elemStiffMatrix[i + i*m_numDofElem] +=  scaleFactorStiff;
	  elemStiffMatrix[j + j*m_numDofElem] +=  scaleFactorStiff;
	  elemStiffMatrix[i + j*m_numDofElem] += -scaleFactorStiff;
	  elemStiffMatrix[j + i*m_numDofElem] += -scaleFactorStiff;
	}
      }
      break;
    case ELASTICITY:
      assert (m_spatialDim != 1);
      for (int i=0; i<numNodePerElem; i++) {
	for (int j=0; j<numNodePerElem; j++) {
	  if (i != j) {
	    addElasticity(i, j, nodes, m_numDofElem, 
			  scaleFactorStiff, elemStiffMatrix);
	  }
	}
      }
      break;
    default:
      break;
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

  void assembleSubdomainMatrix(double hElem,
			       std::vector< std::vector<LO> > &elemConn)
  {
    std::vector< std::vector<LO> > nodalConn;
    determineNodalConnectivity(elemConn, nodalConn);
    determineNodeBegLocalDof();
    if (elemConn.size() == 0) return;
    getElementMatrices(hElem, elemConn[0], m_elemStiffMatrix, m_elemMassMatrix);
    std::vector<SX> elementLoadVector;
    determineElementLoadVector(elementLoadVector);
    double omega = m_Parameters->get("omega", 0.0);
    double betaDamping = m_Parameters->get("betaDamping", 0.0);
    LO numDof = m_nodeBegin[m_numNode];
    m_loadVector.resize(numDof);
    // determine potentially non-zero columns in each row of matrix
    m_rowBegin.resize(numDof+1, 0);
    m_columns.resize(0);
    LO nnz(0);
    for (LO i=0; i<m_numNode; i++) {
      for (LO j=m_nodeBegin[i]; j<m_nodeBegin[i+1]; j++) {
	for (size_t k=0; k<nodalConn[i].size(); k++) {
	  LO node2 = nodalConn[i][k];
	  for (LO m=m_nodeBegin[node2]; m<m_nodeBegin[node2+1]; m++) {
	    m_columns.push_back(m);
	    nnz++;
	  }
	}
	m_rowBegin[j+1] = nnz;
      }
    }
    // element and force assembly
    std::vector<LO> imap(numDof, -1);
    m_values.resize(nnz, 0);
    std::vector<SX> rowValues(m_numDofElem);
    for (LO i=0; i<m_numElem; i++) {
      LO elementRow(0);
      for (size_t j=0; j<elemConn[i].size(); j++) {
	LO node = elemConn[i][j];
	for (LO k=m_nodeBegin[node]; k<m_nodeBegin[node+1]; k++) {
	  m_loadVector[k] += elementLoadVector[k-m_nodeBegin[node]];
	  for (LO m=m_rowBegin[k]; m<m_rowBegin[k+1]; m++) {
	    imap[m_columns[m]] = m;
	  }
	  int elementCol(0);
	  getElementMatrixRow(elementRow, omega, betaDamping, rowValues);
	  for (size_t n=0; n<elemConn[i].size(); n++) {
	    LO node2 = elemConn[i][n];
	    for (LO p=m_nodeBegin[node2]; p<m_nodeBegin[node2+1]; p++) {
	      LO col = p;
	      assert (imap[col] != -1);
	      m_values[imap[col]] += rowValues[elementCol];
	      elementCol++;
	    }
	  }
	  elementRow++;
	  for (LO m=m_rowBegin[k]; m<m_rowBegin[k+1]; m++) {
	    imap[m_columns[m]] = -1;
	  }
	}
      }
    }
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

  void getElementMatrixRow(int localRow, 
			   double omega,
			   double betaDamping,
			   std::vector<SX> & rowValues) const
  {
    SX a(0);
    int index;
    switch(m_analysisType) {
    case STANDARD:
      for (LO i=0; i<m_numDofElem; i++) {
	index = localRow+m_numDofElem*i;
	rowValues[i] = getNumber(m_elemStiffMatrix[index], 0, a);
      }
      break;
    case HELMHOLTZA:
      for (LO i=0; i<m_numDofElem; i++) {
	index = localRow+m_numDofElem*i;
	rowValues[i] =  
	  getNumber(m_elemStiffMatrix[index], 
		    betaDamping*m_elemStiffMatrix[index]*omega, a);
	rowValues[i] -= getNumber(m_elemMassMatrix[index]*omega*omega, 0, a);
      }
      break;
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

  void getModel1D(double & hElem, 
		  std::vector< std::vector<LO> > & elemConn)
  {
    LO numSubDir1 = m_Parameters->get("Number of Subdomains Direction 1", 1);
    double lengthDir1 = m_Parameters->get("Length Direction 1", double(1));
    LO numElemPerSubDir1 = 
      m_Parameters->get("Number of Elements Per Subdomain Direction 1", 1);
    assert (m_numProc >= numSubDir1);
    LO numActiveProc = numSubDir1;
    if (m_myPID > numActiveProc-1) {
      m_numNode = 0;
      m_numElem = 0;
      return;
    }
    m_numNode = numElemPerSubDir1+1;
    m_numElem = numElemPerSubDir1;
    double lengthSubDir1 = lengthDir1/numSubDir1;
    hElem = lengthSubDir1/numElemPerSubDir1;
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
  
  void getModel2D(double & hElem, 
		  std::vector< std::vector<LO> > & elemConn)
  {
    LO numSubDir1 = m_Parameters->get("Number of Subdomains Direction 1", 1);
    LO numSubDir2 = m_Parameters->get("Number of Subdomains Direction 2", 1);
    double lengthDir1 = m_Parameters->get("Length Direction 1", double(1));
    double lengthDir2 = m_Parameters->get("Length Direction 2", double(1));
    LO numElemPerSubDir1 = 
      m_Parameters->get("Number of Elements Per Subdomain Direction 1", 1);
    LO numElemPerSubDir2 = 
      m_Parameters->get("Number of Elements Per Subdomain Direction 2", 1);
    assert (m_numProc >= numSubDir1*numSubDir2);
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
    hElem = std::max(hElem1, hElem2);
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

  void getModel3D(double & hElem, 
		  std::vector< std::vector<LO> > & elemConn)
  {
    LO numSubDir1 = m_Parameters->get("Number of Subdomains Direction 1", 1);
    LO numSubDir2 = m_Parameters->get("Number of Subdomains Direction 2", 1);
    LO numSubDir3 = m_Parameters->get("Number of Subdomains Direction 3", 1);
    double lengthDir1 = m_Parameters->get("Length Direction 1", double(1));
    double lengthDir2 = m_Parameters->get("Length Direction 2", double(1));
    double lengthDir3 = m_Parameters->get("Length Direction 3", double(1));
    LO numElemPerSubDir1 = 
      m_Parameters->get("Number of Elements Per Subdomain Direction 1", 1);
    LO numElemPerSubDir2 = 
      m_Parameters->get("Number of Elements Per Subdomain Direction 2", 1);
    LO numElemPerSubDir3 = 
      m_Parameters->get("Number of Elements Per Subdomain Direction 3", 1);
    assert (m_numProc >= numSubDir1*numSubDir2*numSubDir3);
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
    hElem = std::max(hElem1, hElem2);
    hElem = std::max(hElem, hElem3);
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

  void determineNodalConnectivity(
		  std::vector< std::vector<LO> > & elemConn,
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
  
  void GenerateConstraintEquations()
  {
    //
    // Notes: All constraints are for subdomain 2. Let u_L denote the vector of 
    // dofs for nodes on the left side of the second subdomain. Similarly, Let 
    // u_L+1 denote the vector of dofs for nodes one layer to the right of the 
    // left side of the second subdomain. The constraints are u_L+1 = u_L.
    //
    if (m_myPID != 2) return;
    if (!m_Parameters->get("Generate Constraint Equations", false)) {
      return;
    }
    SM xMin = m_x[0];
    std::vector<LO> leftSideNodes;
    for (int i=0; i<m_numNode; i++) {
      if (m_x[i] == xMin) leftSideNodes.push_back(i);
    }
    LO numConstraints(0);
    m_rowBeginCT.resize(1, 0);
    for (size_t i=0; i<leftSideNodes.size(); i++) {
      LO node1 = leftSideNodes[i];
      size_t numDofNode1 = m_nodeBegin[node1+1]-m_nodeBegin[node1];
      LO node2 = node1+1;
      size_t numDofNode2 = m_nodeBegin[node2+1]-m_nodeBegin[node2];
      assert (numDofNode2 == numDofNode1);
      for (size_t j=0; j<numDofNode1; j++) {
	LO dof = m_nodeBegin[node1]+j;
	m_rowsCT.push_back(dof);
	m_matrixCT.push_back(1);
	m_rowBeginCT.push_back(m_rowBeginCT.size());
	LO dof2 = m_nodeBegin[node2]+j;
	m_rowsCT.push_back(dof2);
	m_matrixCT.push_back(-1);
	m_colsCT.push_back(numConstraints);
	numConstraints++;
      }
    }
  }

  void printElementMatrix(const SM* matrix,
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

  void determineElementLoadVector(std::vector<SX> & elementLoadVector) const
  {
    int numDofPerNode = getNumDofPerNode();
    int numNodePerElem = power2(m_spatialDim);
    assert (m_numDofElem == numNodePerElem*numDofPerNode);
    elementLoadVector.resize(m_numDofElem);
    int loadDirection = m_Parameters->get("Load Direction", 1);
    for (int i=0; i<numNodePerElem; i++) {
      for (int j=0; j<numDofPerNode; j++) {
	if (j == loadDirection) {
	  LO index = i*numDofPerNode+j;
	  SX rowSum(0);
	  for (int k=0; k<m_numDofElem; k++) {
	    LO indexElem = index+k*m_numDofElem;
	    rowSum += m_elemMassMatrix[indexElem];
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
  int m_myPID, m_numProc, m_numDofElem;
  LO m_numNode, m_numElem;
  std::vector<double> m_x, m_y, m_z, m_elemStiffMatrix, m_elemMassMatrix;
  std::vector<LO> m_nodeBegin, m_localDofs, m_rowBegin, m_columns;
  std::vector<GO> m_nodeGIDs;
  std::vector<SX> m_values, m_loadVector;
  std::vector<LO> m_rowsCT, m_colsCT, m_rowBeginCT, m_columnsCT;
  std::vector<SM> m_matrixCT;
  std::vector< std::vector<LO> > m_elemConn;
  double m_hElem;
};

} // namespace bddc

#endif // PROBLEMMAKERBDDC_H
