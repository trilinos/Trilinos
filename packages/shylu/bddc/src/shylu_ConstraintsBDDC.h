
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

#ifndef CONSTRAINTSBDDC_H
#define CONSTRAINTSBDDC_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>
#include "shylu_enumsBDDC.h"

#include "shylu_PartitionOfUnityBDDC.h"
#include "shylu_SubdomainBDDC.h"
#include "shylu_UtilBDDC.h"

using Teuchos::RCP;
using Teuchos::rcp;

// Author: Clark R. Dohrmann
namespace bddc {
  
template <class SX, class SM, class LO, class GO> 
  class ConstraintsBDDC
{
public:
  //
  // Convenience typedefs
  //
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;
  typedef Tpetra::Map<LO,GO,Node>                                 Map;
  typedef Tpetra::CrsGraph<LO,GO,Node>                            CrsGraph;
  typedef Tpetra::CrsMatrix<SX,LO,GO,Node>                        CrsMatrix;
  typedef Tpetra::CrsMatrix<GO,LO,GO,Node>                        CrsMatrixGO;
  typedef Tpetra::Export<LO,GO,Node>                              Export;
  typedef Tpetra::Import<LO,GO,Node>                              Import;
  typedef Tpetra::Vector<SX,LO,GO,Node>                           Vector;
  typedef Tpetra::Vector<GO,LO,GO,Node>                           VectorGO;
  typedef Tpetra::MultiVector<SX,LO,GO,Node>                      MV;
  typedef Tpetra::MultiVector<GO,LO,GO,Node>                      MVGO;

  ConstraintsBDDC()
  {
  }
  ConstraintsBDDC
    (LO numNodes,
     const LO* nodeBegin,
     const LO* localDofs,
     const std::vector< std::vector<LO> > & subNodes,
     std::vector< SubdomainBDDC<SX,SM,LO,GO>* > & Subdomain,
     RCP< PartitionOfUnity<SX,SM,LO,GO> > & Partition,
     RCP<Export> & exporterB,
     RCP<Import> & importerB,
     const std::vector<SM> & diagBoundary,
     RCP<Teuchos::ParameterList> & Parameters) :
  m_numNodes(numNodes),
    m_nodeBegin(nodeBegin),
    m_localDofs(localDofs),
    m_subNodes(subNodes),
    m_Subdomain(Subdomain),
    m_Partition(Partition),
    m_exporterB(exporterB),
    m_importerB(importerB),
    m_dofMapB(exporterB->getSourceMap()),
    m_dofMapB1to1(importerB->getSourceMap()),
    m_Comm(exporterB->getSourceMap()->getComm()),
    m_diagBoundary(diagBoundary),
    m_Parameters(Parameters),
    m_numSub(Subdomain.size()),
    m_numDofsB(m_dofMapB->getNodeNumElements()),
    m_numDofs(nodeBegin[numNodes]),
    m_spatialDim(Parameters->get("Spatial Dimension", 3)),
    m_problemType(Parameters->get("Problem Type BDDC", SCALARPDE)),
    m_IGO(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid())
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
      std::vector< std::vector<LO> > equivClassDofs(numEquiv);
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
	determineEquivConstraints(i, equivNodes, equiv, 
				  dofMapEquiv, numCorners, numEdges, numFaces,
				  numEquivConstraints[j], equivConstraints[j]);
	resetDofMap(equivNodes, dofMapEquiv);
      }
      m_Subdomain[i]->setEquivalenceClasses(equivClassDofs);
      m_Subdomain[i]->setInitialConstraints(numEquivConstraints,
					    equivConstraints);
      resetDofMap(m_subNodes[i], dofMap);
    }
  }

private:
  LO m_numNodes;
  const LO* m_nodeBegin;
  const LO* m_localDofs;
  const std::vector< std::vector<LO> > & m_subNodes;
  std::vector< SubdomainBDDC<SX,SM,LO,GO>* > & m_Subdomain;
  RCP< PartitionOfUnity<SX,SM,LO,GO> > & m_Partition;
  RCP<Export> m_exporterB;
  RCP<Import> m_importerB;
  RCP<const Map> m_dofMapB, m_dofMapB1to1;
  RCP<const Teuchos::Comm<int> > m_Comm;
  const std::vector<SM> m_diagBoundary;
  RCP<Teuchos::ParameterList> & m_Parameters;
  LO m_numSub, m_numDofsB, m_numDofs, m_spatialDim;
  enum ProblemType m_problemType;
  Tpetra::global_size_t m_IGO;

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
	assert (dof != -1);
	equivDofs.push_back(dof);
      }
    }
  }

  void determineEquivConstraints(LO sub,
				 const std::vector<LO> & equivNodes,
				 LO equiv,
				 std::vector<LO> & dofMapEquiv,
				 LO numCorners, 
				 LO numEdges, 
				 LO numFaces,
				 LO & numEquivConstraints, 
				 std::vector<SX> & equivConstraints)
  {
    enum EquivType equivType = m_Partition->getEquivType(equiv);
    // really only doing anything for corners right now
    if (equivType == CORNER) {
      LO node = equivNodes[0];
      LO numDof = m_nodeBegin[node+1] - m_nodeBegin[node];
      numEquivConstraints = numDof;
      equivConstraints.resize(numDof*numDof);
      for (LO i=0; i<numDof; i++) {
	LO dof = m_nodeBegin[node] + i;
	assert (dofMapEquiv[dof] != -1);
	equivConstraints[i+i*numDof] = 1;
      }
    }
    else if (((equivType == FACE) && (m_spatialDim == 3)) ||
	     ((equivType == EDGE) && (m_spatialDim == 2))) {
      if (m_Parameters->get("Use Flux Constraints", false) == true) {
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
	  numEquivConstraints = 1; // change later to account for irregular mesh decomps
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
	if (m_problemType == SCALARPDE) {
	  if (numCorners < 1) {
	  }
	}
	if (m_problemType == ELASTICITY) {
	  if (m_spatialDim == 3) {
	    if (numCorners <= 2) {
	      if (numCorners == 1) {
	      }
	      else {
	      }
	    }
	  }
	}
      }
    }
    else if (equivType == EDGE) {
      if (m_spatialDim == 2) {
	if (m_problemType == SCALARPDE) {
	  if (numCorners < 1) {
	  }
	}
	if (m_problemType == ELASTICITY) {
	  if (numCorners < 2) {
	    if (numCorners == 1) {
	    }
	    else {
	    }
	  }
	}
      }
    }
  }

};

} // namespace bddc

#endif // CONSTRAINTSBDDC_H
  
