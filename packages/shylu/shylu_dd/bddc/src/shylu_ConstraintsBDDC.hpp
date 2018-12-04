
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

#ifndef BDDC_CONSTRAINTS_H
#define BDDC_CONSTRAINTS_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "shylu_enumsBDDC.hpp"
#include "shylu_PartitionOfUnityBDDC.hpp"
#include "shylu_SubdomainBDDC.hpp"
#include "shylu_UtilBDDC.hpp"
#include "shylu_errorBDDC.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace bddc {
  
template <class SX, class SM, class LO, class GO> 
  class ConstraintsBDDC
{
public:

  ConstraintsBDDC()
  {
  }
  ConstraintsBDDC
    (LO numNodes,
     const LO* nodeBegin,
     const LO* localDofs,
     const SM* xCoord,
     const SM* yCoord,
     const SM* zCoord,
     const std::vector< std::vector<LO> > & subNodes,
     std::vector< SubdomainBDDC<SX,SM,LO,GO>* > & Subdomain,
     RCP< PartitionOfUnity<SX,SM,LO,GO> > & Partition,
     const std::vector<SM> & diagBoundary,
     RCP<Teuchos::ParameterList> & Parameters) :
  m_numNodes(numNodes),
    m_nodeBegin(nodeBegin),
    m_localDofs(localDofs),
    m_xCoord(xCoord),
    m_yCoord(yCoord),
    m_zCoord(zCoord),
    m_subNodes(subNodes),
    m_Subdomain(Subdomain),
    m_Partition(Partition),
    m_diagBoundary(diagBoundary),
    m_Parameters(Parameters),
    m_numSub(Subdomain.size()),
    m_numDofs(nodeBegin[numNodes]),
    m_spatialDim(Parameters->get("Spatial Dimension", 3)),
    m_problemType(Parameters->get("Problem Type BDDC", SCALARPDE)),
    m_IGO(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid()),
    m_useCorners(Parameters->get("Use Corners", false)),
    m_useEdges(Parameters->get("Use Edges", false)),
    m_useFaces(Parameters->get("Use Faces", false))
  {
  }

  void determineBaseConstraints()
  {
    const std::vector< std::vector<LO> > & subdomainEquivClasses = 
      m_Partition->getSubdomainEquivClasses();
    const std::vector< std::vector<LO> > & equivClasses = 
      m_Partition->getEquivClasses();
    const std::vector<LO> & equivCard = m_Partition->getEquivCardinality();
    std::vector<LO> dofMap(m_numDofs, -1), dofMapEquiv(m_numDofs, -1);
    for (LO i=0; i<m_numSub; i++) {
      LO numEquiv = subdomainEquivClasses[i].size();
      setDofMap(m_subNodes[i], dofMap);
      std::vector< std::vector<LO> > equivClassDofs(numEquiv),
	equivConstraintsLocalDofs(numEquiv);
      std::vector< std::vector<SX> > equivConstraints(numEquiv);
      std::vector<LO> numEquivConstraints(numEquiv);
      LO numCorners(0), numEdges(0), numFaces(0);
      getNumCornersEdgesFaces(subdomainEquivClasses[i], equivClasses,
			      equivCard, numCorners, numEdges, numFaces);
      for (LO j=0; j<numEquiv; j++) {
	LO equiv = subdomainEquivClasses[i][j];
	const std::vector<LO> & equivNodes = equivClasses[equiv];
	getDofsEquiv(equivNodes, dofMap, equivClassDofs[j]);
	setDofMap(equivNodes, dofMapEquiv);
	determineEquivConstraints
	  (i, equivNodes, equiv, dofMapEquiv, numCorners, numEdges, numFaces,
	   equivConstraintsLocalDofs[j], equivConstraints[j]);
	resetDofMap(equivNodes, dofMapEquiv);
      }
      m_Subdomain[i]->setEquivalenceClasses(equivClassDofs);
      m_Subdomain[i]->setInitialConstraints(equivConstraintsLocalDofs,
					    equivConstraints);
      resetDofMap(m_subNodes[i], dofMap);
    }
  }

private:
  LO m_numNodes;
  const LO* m_nodeBegin;
  const LO* m_localDofs;
  const SM *m_xCoord, *m_yCoord, *m_zCoord;
  const std::vector< std::vector<LO> > & m_subNodes;
  std::vector< SubdomainBDDC<SX,SM,LO,GO>* > & m_Subdomain;
  RCP< PartitionOfUnity<SX,SM,LO,GO> > & m_Partition;
  const std::vector<SM> m_diagBoundary;
  RCP<Teuchos::ParameterList> & m_Parameters;
  LO m_numSub, m_numDofs, m_spatialDim;
  enum ProblemType m_problemType;
  Tpetra::global_size_t m_IGO;
  bool m_useCorners, m_useEdges, m_useFaces;

  void setDofMap(const std::vector<LO> & nodes, 
		 std::vector<LO> & dofMap)
  {
    LO numDof = 0;
    for (size_t i=0; i<nodes.size(); i++) {
      LO node = nodes[i];
      for (LO j=m_nodeBegin[node]; j<m_nodeBegin[node+1]; j++) {
	dofMap[j] = numDof++;
      }
    }
  }

  void resetDofMap(const std::vector<LO> & nodes, 
		   std::vector<LO> & dofMap)
  {
    for (size_t i=0; i<nodes.size(); i++) {
      LO node = nodes[i];
      for (LO j=m_nodeBegin[node]; j<m_nodeBegin[node+1]; j++) {
	dofMap[j] = -1;
      }
    }
  }

  void getNumCornersEdgesFaces
    (const std::vector<LO> & subEquivClasses, 
     const std::vector< std::vector<LO> > & equivClasses,
     const std::vector<LO> & equivCard, 
     LO & numCorners, 
     LO & numEdges, 
     LO & numFaces)
  {
    numCorners = numEdges = numFaces = 0;
    LO numEquiv = subEquivClasses.size();
    for (LO i=0; i<numEquiv; i++) {
      LO equiv = subEquivClasses[i];
      enum EquivType equivType = m_Partition->getEquivType(equiv);
      if (equivType == CORNER) numCorners++;
      if (equivType == EDGE) numEdges++;
      if (equivType == FACE) numFaces++;
    }
  }

  void getDofsEquiv(const std::vector<LO> & equivClass,
		    std::vector<LO> & dofMap,
		    std::vector<LO> & equivDofs)
  {
    LO numNodes = equivClass.size();
    for (LO i=0; i<numNodes; i++) {
      LO node = equivClass[i];
      for (LO j=m_nodeBegin[node]; j<m_nodeBegin[node+1]; j++) {
	LO dof = dofMap[j];
	BDDC_TEST_FOR_EXCEPTION(dof == -1, std::runtime_error, "invalid dof");
	equivDofs.push_back(dof);
      }
    }
  }

  void getLocalDofsAndCoordinates(const std::vector<LO> & equivNodes, 
				  std::vector<LO> & equivLocalDofs, 
				  std::vector<SM> & xCoords,
				  std::vector<SM> & yCoords,
				  std::vector<SM> & zCoords)
  {
    LO numRows(0);
    for (size_t i=0; i<equivNodes.size(); i++) {
      LO node = equivNodes[i];
      numRows += m_nodeBegin[node+1] - m_nodeBegin[node];
    }
    equivLocalDofs.resize(numRows);
    xCoords.resize(numRows);
    yCoords.resize(numRows); 
    zCoords.resize(numRows);
    numRows = 0;
    SM xSum(0), ySum(0), zSum(0);
    for (size_t i=0; i<equivNodes.size(); i++) {
      LO node = equivNodes[i];
      for (LO j=m_nodeBegin[node]; j<m_nodeBegin[node+1]; j++) {
	int localDof = m_localDofs[j];
	equivLocalDofs[numRows] = localDof;
	xCoords[numRows] = m_xCoord[node]; xSum += m_xCoord[node];
	yCoords[numRows] = m_yCoord[node]; ySum += m_yCoord[node];
	zCoords[numRows] = m_zCoord[node]; zSum += m_zCoord[node];
	numRows++;
      }
    }
    SM xCent = xSum/numRows;
    SM yCent = ySum/numRows;
    SM zCent = zSum/numRows;
    for (LO i=0; i<numRows; i++) {
      xCoords[i] -= xCent;
      yCoords[i] -= yCent;
      zCoords[i] -= zCent;
    }
  }

  int getNumNeededAncestors()
  {
    int numNeededAncestors(0);
    if (m_problemType == SCALARPDE) {
      numNeededAncestors = 1;
    }
    else if (m_problemType == ELASTICITY) {
      if (m_spatialDim == 2) numNeededAncestors = 2;
      if (m_spatialDim == 3) numNeededAncestors = 3;
    }
    return numNeededAncestors;
  }

  void getNullSpace(const std::vector<LO> & localDofs, 
		    const std::vector<SM> & xCoords, 
		    const std::vector<SM> & yCoords, 
		    const std::vector<SM> & zCoords, 
		    std::vector<SX> & nullSpace)
  {
    const LO numRows = localDofs.size();
    if (m_problemType == SCALARPDE) {
      nullSpace.resize(numRows, 1);
    }
    else if (m_problemType == ELASTICITY) {
      if (m_spatialDim == 2) {
	int numRBM = 3;
	LO numTerms = numRBM*numRows;
	nullSpace.resize(numTerms, 0);
	numTerms = 0;
	for (int j=0; j<numRBM; j++) {
	  for (LO i=0; i<numRows; i++) {
	    nullSpace[numTerms++] = 
	      getElasticityRBM2D(j, localDofs[i],
				 xCoords[i], yCoords[i], zCoords[i]);
	  }
	}
      }
      else if (m_spatialDim == 3) {
	int numRBM = 6;
	LO numTerms = numRBM*numRows;
	nullSpace.resize(numTerms, 0);
	numTerms = 0;
	for (int j=0; j<numRBM; j++) {
	  for (LO i=0; i<numRows; i++) {
	    nullSpace[numTerms++] = 
	      getElasticityRBM3D(j, localDofs[i],
				 xCoords[i], yCoords[i], zCoords[i]);
	  }
	}
      }
    }
  }

  SM getElasticityRBM2D(int rbm, 
			int localDof,
			const SM xCoord, 
			const SM yCoord, 
			const SM zCoord)
  {
    SM value = 0.0;
    if (localDof == rbm) value = 1.0;
    if (rbm == 2) {
      if (localDof == 0)      value = -yCoord;
      else if (localDof == 1) value =  xCoord;
    }
    return value;
  }

  SM getElasticityRBM3D(int rbm, 
			int localDof,
			const SM xCoord, 
			const SM yCoord, 
			const SM zCoord)
  {
    SM value = 0.0;
    if (localDof == rbm) value = 1.0;
    switch(rbm) {
    case 3:
      if (localDof == 1)      value = -zCoord;
      else if (localDof == 2) value =  yCoord;
      break;
    case 4:
      if (localDof == 2)      value = -xCoord;
      else if (localDof == 0) value =  zCoord;
      break;
    case 5:
      if (localDof == 0)      value = -yCoord;
      else if (localDof == 1) value =  xCoord;
      break;
    default:
      break;
    }
    return value;
  }

  void determineIndependentColumns(const LO numRows, 
				   const std::vector<SX> & nullSpace, 
				   const int numMissing,
				   const enum EquivType equivType,
				   std::vector<int> & independentCols)
  {
    if (numRows == 0) return;
    int numColsMax = nullSpace.size()/numRows;
    if (m_problemType == SCALARPDE) {
      numColsMax = 1;
    }
    else if (m_problemType == ELASTICITY) {
      if (numMissing <= 1) numColsMax = m_spatialDim;
      if ((m_spatialDim == 3) && (equivType == EDGE)) {
	numColsMax = std::min(numColsMax, m_spatialDim);
      }
    }
    std::vector<SX> A = nullSpace;
    SM tol(1e-10);
    Teuchos::BLAS<int, SX>  BLAS;
    int INCX(1);
    for (int i=0; i<numColsMax; i++) {
      SX* col = &A[numRows*i];
      SM norm = BLAS.NRM2(numRows, col, INCX);
      if (norm > tol) {
	independentCols.push_back(i);
	for (LO j=0; j<numRows; j++) col[j] /= norm;
	for (int j=i+1; j<numColsMax; j++) {
	  SX* col2 = &A[numRows*j];
	  SX dotProd = 0;
	  for (LO k=0; k<numRows; k++) dotProd += col[k]*col2[k];
	  for (LO k=0; k<numRows; k++) col2[k] -= dotProd*col[k];
	}
      }
    }
  }

  void determineEquivConstraints(const LO numRows, 
				 const std::vector<SX> & nullSpace, 
				 const std::vector<int> & independentCols, 
				 std::vector<SX> & equivConstraints)
  {
    const int numCols = independentCols.size();
    std::vector<SX> A(numRows*numCols);
    for (int j=0; j<numCols; j++) {
      const int col = independentCols[j];
      for (int i=0; i<numRows; i++) {
	A[j*numRows+i] = nullSpace[col*numRows+i];
      }
    }
    // The the pseudo-inverse of A (transpose), obtained using the singular 
    // value decomposition, is used for the equivalence class constraints
    int M = numRows;
    int N = numCols;
    std::vector<SM> S(std::max(M, N));
    std::vector<SM> RWORK(5*std::min(M, N));
    std::vector<SX> VT(N*N);
    char JOBU('O'), JOBVT('A');
    SX U[1], dumWORK[1];
    int LDU(1), INFO, LDVT(N), LWORK(-1);
    Teuchos::LAPACK<int, SX> LAPACK;
    LAPACK.GESVD(JOBU, JOBVT, M, N, A.data(), M, S.data(), U, LDU, 
		 VT.data(), LDVT, dumWORK, LWORK, RWORK.data(), &INFO);
    BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GESVD error");
    LWORK = dumWORK[0];
    std::vector<SX> WORK(LWORK);
    LAPACK.GESVD(JOBU, JOBVT, M, N, A.data(), M, S.data(), U, LDU, 
		 VT.data(), LDVT, WORK.data(), LWORK, RWORK.data(), &INFO);
    // Note: for A = U * S * V^T, the matrix V^T is stored one column at
    //       a time in VT (regular Fortran column ordering for transpose of V)
    BDDC_TEST_FOR_EXCEPTION(INFO != 0, std::runtime_error, "GESVD error");
    // scale columns of U by inverse of singular values
    for (int j=0; j<N; j++) {
      SX* colU = &A[M*j];
      for (LO i=0; i<numRows; i++) {
	colU[i] /= S[j];
      }
    }
    equivConstraints.resize(M*N);
    Teuchos::BLAS<int, SX>  BLAS;
    SX ALPHA(1), BETA(0);
    BLAS.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, M, N, N, ALPHA, A.data(), 
	      M, VT.data(), N, BETA, equivConstraints.data(), M);
  }

  void determineEquivConstraints(LO sub,
				 const std::vector<LO> & equivNodes,
				 LO equiv,
				 std::vector<LO> & dofMapEquiv,
				 LO numCorners, 
				 LO numEdges, 
				 LO numFaces,
				 std::vector<LO> & equivConstraintsLocalDofs, 
				 std::vector<SX> & equivConstraints)
  {
    enum EquivType equivType = m_Partition->getEquivType(equiv);
    int numAncestors = m_Partition->getNumActiveAncestors(equiv);
    int numNeeded = getNumNeededAncestors();
    int numMissing = numNeeded - numAncestors;
    bool isActive = false;

    if (equivType == CORNER) {
      numMissing = 1;
      if (m_useCorners == true) isActive = true;
    }
    else if (equivType == EDGE) {
      if (m_useEdges == true) isActive = true;
    }
    else if (equivType == FACE) {
      if (m_useFaces == true) isActive = true;
      if (numMissing > 0) isActive = true;
    }
    int equivFlag = m_Partition->getEquivFlag(equiv);
    if (equivFlag > 0) isActive = true;

    numMissing = std::max(numMissing, 0);
    if (isActive == true) {
      std::vector<LO> equivLocalDofs;
      std::vector<SM> xCoords, yCoords, zCoords;
      getLocalDofsAndCoordinates(equivNodes, equivLocalDofs, 
				 xCoords, yCoords, zCoords);
      const LO numRows = equivLocalDofs.size();
      std::vector<SX> nullSpace;
      getNullSpace(equivLocalDofs, xCoords, yCoords, zCoords, nullSpace); 
      std::vector<int> independentCols;
      determineIndependentColumns(numRows, nullSpace, numMissing, equivType,
				  independentCols);
      determineEquivConstraints(numRows, nullSpace, independentCols, 
				equivConstraints);
      equivConstraintsLocalDofs = independentCols;      
    }
  }

  void processFluxConstraints(const std::vector<LO> & equivNodes,
			      std::vector<SX> & equivConstraints)
  {
    LO numNode = equivNodes.size();
    LO numFlux(0), numDof(0);
    for (LO i=0; i<numNode; i++) {
      LO node = equivNodes[i];
      for (LO j=m_nodeBegin[node]; j<m_nodeBegin[node+1]; j++) {
	if (m_localDofs[j] == 7) numFlux++;
	numDof++;
      }
    }
    if (numFlux > 0) {
      LO numEquivConstraints = 1; // change later to account for irregular mesh decomps
      equivConstraints.resize(numDof);
      numDof = 0;
      for (LO i=0; i<numNode; i++) {
	LO node = equivNodes[i];
	for (LO j=m_nodeBegin[node]; j<m_nodeBegin[node+1]; j++) {
	  if (m_localDofs[j] == 7) equivConstraints[numDof] = 1;
	  numDof++;
	}
      }
    }
  }

};

} // namespace bddc

#endif // BDDC_CONSTRAINTS_H
  
