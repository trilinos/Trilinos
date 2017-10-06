
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

#ifndef WEIGHTSBDDC_H
#define WEIGHTSBDDC_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
#include <vector>
#include "shylu_enumsBDDC.h"

#include "shylu_PartitionOfUnityBDDC.h"
#include "shylu_SubdomainBDDC.h"

using Teuchos::RCP;
using Teuchos::rcp;

// Author: Clark R. Dohrmann
namespace bddc {
  
template <class SX, class SM, class LO, class GO> 
  class WeightsBDDC
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
  typedef Tpetra::MultiVector<SX,LO,GO,Node>                      MV;

  WeightsBDDC()
  {
  }
  WeightsBDDC
    (std::vector< SubdomainBDDC<SX,SM,LO,GO>* > & Subdomain,
     RCP< PartitionOfUnity<SX,SM,LO,GO> > & Partition,
     RCP<Export> & exporterB,
     const std::vector< std::vector<LO> > & subBoundaryDofs,
     const std::vector<SM> & diagBoundary,
     RCP<Teuchos::ParameterList> & Parameters) :
  m_Subdomain(Subdomain),
    m_Partition(Partition),
    m_exporterB(exporterB),
    m_dofMapB(exporterB->getSourceMap()),
    m_dofMapB1to1(exporterB->getTargetMap()),
    m_Comm(exporterB->getSourceMap()->getComm()),
    m_subBoundaryDofs(subBoundaryDofs),
    m_diagBoundary(diagBoundary),
    m_Parameters(Parameters),
    m_numSub(Subdomain.size()),
    m_numDofsB(m_dofMapB->getNodeNumElements())
  {
  }

  void determineWeights()
  {
    determineEquivAndWeightTypes();
    std::vector< std::vector< std::vector<SX> > > deluxeWeights;
    determineDeluxeWeights(deluxeWeights);
    determineEquivWeights(deluxeWeights);
  }

  SM checkWeights()
  {
    return checkBDDCWeights();
  }

private:
  std::vector< SubdomainBDDC<SX,SM,LO,GO>* > & m_Subdomain;
  RCP< PartitionOfUnity<SX,SM,LO,GO> > & m_Partition;
  RCP<Export> m_exporterB;
  RCP<const Map> m_dofMapB, m_dofMapB1to1;
  RCP<const Teuchos::Comm<int> > m_Comm;
  const std::vector< std::vector<LO> > m_subBoundaryDofs;
  const std::vector<SM> m_diagBoundary;
  RCP<Teuchos::ParameterList> & m_Parameters;
  LO m_numSub, m_numDofsB;
  std::vector< std::vector<enum bddc::EquivType> > m_equivType;
  std::vector< std::vector<enum bddc::WeightType> > m_weightType;

  void determineEquivAndWeightTypes()
  {
    const std::vector< std::vector<LO> > & subdomainEquivClasses = 
      m_Partition->getSubdomainEquivClasses();
    m_equivType.resize(m_numSub);
    m_weightType.resize(m_numSub);
    for (LO i=0; i<m_numSub; i++) {
      LO numEquiv = m_Subdomain[i]->getNumEquiv();
      m_equivType[i].resize(numEquiv);
      m_weightType[i].resize(numEquiv);
      for (LO j=0; j<numEquiv; j++) {
	LO equiv = subdomainEquivClasses[i][j];
	m_equivType[i][j] = m_Partition->getEquivType(equiv);
	m_weightType[i][j] = getWeightType(m_equivType[i][j]);
      }
    }
  }

  enum WeightType getWeightType(enum EquivType equivType)
  {
    enum WeightType weightType;
    if (equivType == CORNER) {
      weightType = m_Parameters->get("Weight Type Corner", STIFFNESS);
    }
    else if (equivType == EDGE) {
      weightType = m_Parameters->get("Weight Type Edge", STIFFNESS);
    }
    else {
      weightType = m_Parameters->get("Weight Type Face", STIFFNESS);
    }
    return weightType;
  }

  void determineDeluxeWeights
    (std::vector< std::vector< std::vector<SX> > > & deluxeWeights)
  {
    // quick return if no deluxe scaling on any equivalance classes
    LO deluxeFlag(0);
    for (LO i=0; i<m_numSub; i++) {
      LO numEquiv = m_Subdomain[i]->getNumEquiv();
      for (LO j=0; j<numEquiv; j++) {
	if (m_weightType[i][j] == DELUXE) deluxeFlag = 1;
      }
    }
    LO deluxeFlagMax;
    Teuchos::reduceAll<int, LO> (*m_Comm, Teuchos::REDUCE_MAX, 1, 
				 &deluxeFlag, &deluxeFlagMax);
    if (deluxeFlagMax == 0) return;
    // calculate equivalence class Schur complements
    deluxeWeights.resize(m_numSub);
    std::vector<bool> rowFlag(m_numDofsB, false);
    std::vector<LO> columns;
    RCP<CrsGraph> ScGraph = rcp( new CrsGraph(m_dofMapB, m_dofMapB, 0) );
    for (LO i=0; i<m_numSub; i++) {
      LO numEquiv = m_Subdomain[i]->getNumEquiv();
      deluxeWeights[i].resize(numEquiv);
      LO localDofB = 0;
      for (LO j=0; j<numEquiv; j++) {
	const LO* equivDofsGlobal = &m_subBoundaryDofs[i][localDofB];
	LO numEquivDofs = m_Subdomain[i]->getNumEquivDofs(j);
	if (m_weightType[i][j] == bddc::DELUXE) {
	  std::vector<SX> & Sc = deluxeWeights[i][j];
	  m_Subdomain[i]->calculateSchurComplement(j, Sc);
	  columns.resize(numEquivDofs);
	  for (LO k=0; k<numEquivDofs; k++) {
	    columns[k] = equivDofsGlobal[k];
	  }
	  for (LO k=0; k<numEquivDofs; k++) {
	    LO row = equivDofsGlobal[k];
	    if (rowFlag[row] == false) {
	      rowFlag[row] = true;
	      ScGraph->insertLocalIndices(row, Teuchos::ArrayView<LO>(columns));
	    }
	  }
	}
	localDofB += numEquivDofs;
      }
    }
    ScGraph->fillComplete(m_dofMapB1to1, m_dofMapB1to1);
    // assemble on-processor Schur complements
    CrsMatrix Sc(ScGraph);
    std::vector<SX> values;
    std::vector<LO> indices;
    for (LO i=0; i<m_numSub; i++) {
      LO numEquiv = m_Subdomain[i]->getNumEquiv();
      LO localDofB = 0;
      for (LO j=0; j<numEquiv; j++) {
	const LO* equivDofsGlobal = &m_subBoundaryDofs[i][localDofB];
	LO numEquivDofs = m_Subdomain[i]->getNumEquivDofs(j);
	std::vector<SX> & ScSub = deluxeWeights[i][j];
	if (ScSub.size() > 0) {
	  values.resize(numEquivDofs);
	  indices.resize(numEquivDofs);
	  for (LO k=0; k<numEquivDofs; k++) {
	    for (LO m=0; m<numEquivDofs; m++) {
	      values[m] = ScSub[k+m*numEquivDofs];
	      indices[m] = equivDofsGlobal[m];
	    }
	    LO row = equivDofsGlobal[k];
	    Sc.sumIntoLocalValues(row, Teuchos::ArrayView<LO>(indices),
				  Teuchos::ArrayView<SX>(values));
	  }
	}
	localDofB += numEquivDofs;
      }
    }
    Sc.fillComplete(m_dofMapB1to1, m_dofMapB1to1);
    // assemble Schur complements across processors
    CrsMatrix Sc1to1(m_dofMapB1to1, 0);
    Sc1to1.doExport(Sc, *m_exporterB, Tpetra::ADD);
    Sc1to1.fillComplete(m_dofMapB1to1, m_dofMapB1to1);
    CrsMatrix ScSum(ScGraph);
    ScSum.doImport(Sc1to1, *m_exporterB, Tpetra::ADD);
    ScSum.fillComplete(m_dofMapB1to1, m_dofMapB1to1);
    // calculate deluxe weight matrices
    std::vector<SX> ScEquiv;
    std::vector<LO> globalToLocalMap(m_numDofsB, -1);
    for (LO i=0; i<m_numSub; i++) {
      LO numEquiv = m_Subdomain[i]->getNumEquiv();
      LO localDofB = 0;
      for (LO j=0; j<numEquiv; j++) {
	const LO* equivDofsGlobal = &m_subBoundaryDofs[i][localDofB];
	LO numEquivDofs = m_Subdomain[i]->getNumEquivDofs(j);
	std::vector<SX> & ScSub = deluxeWeights[i][j];
	if (ScSub.size() > 0) {
	  Teuchos::ArrayView<const LO> Indices;
	  Teuchos::ArrayView<const SX> Values;
	  ScEquiv.resize(numEquivDofs*numEquivDofs);
	  for (LO k=0; k<numEquivDofs; k++) {
	    globalToLocalMap[equivDofsGlobal[k]] = k;
	  }
	  for (LO k=0; k<numEquivDofs; k++) {
	    LO row = equivDofsGlobal[k];
	    ScSum.getLocalRowView(row, Indices, Values);
	    assert (Indices.size() == numEquivDofs);
	    for (LO m=0; m<numEquivDofs; m++) {
	      LO localCol = globalToLocalMap[Indices[m]];
	      assert (localCol != -1);
	      ScEquiv[k+localCol*numEquivDofs] = Values[m];
	    }
	  }
	  for (LO k=0; k<numEquivDofs; k++) {
	    globalToLocalMap[equivDofsGlobal[k]] = -1;
	  }
	  Teuchos::LAPACK<int, SX> LAPACK;
	  int INFO(0), N(numEquivDofs);
	  std::vector<int> IPIV(N);
	  LAPACK.GESV(N, N, &ScEquiv[0], N, &IPIV[0], &deluxeWeights[i][j][0], 
		      N, &INFO);
	}
	localDofB += numEquivDofs;
      }
    }
  }

  void determineEquivWeights
    (std::vector< std::vector< std::vector<SX> > > & deluxeWeights)
  {
    const std::vector< std::vector<LO> > & subdomainEquivClasses = 
      m_Partition->getSubdomainEquivClasses();
    const std::vector<LO> & equivCard = m_Partition->getEquivCardinality();
    for (LO i=0; i<m_numSub; i++) {
      LO numEquiv = m_Subdomain[i]->getNumEquiv();
      std::vector< std::vector<LO> > rowBeginWeight(numEquiv), 
	columnsWeight(numEquiv);
      std::vector< std::vector<SX> > valuesWeight(numEquiv);
      const LO* equivDofsGlobal = &m_subBoundaryDofs[i][0];
      for (LO j=0; j<numEquiv; j++) {
	LO equiv = subdomainEquivClasses[i][j];
	LO numEquivDofs = m_Subdomain[i]->getNumEquivDofs(j);
	if (m_weightType[i][j] == CARDINALITY) {
	  setCardinalityWeights(i, j, rowBeginWeight[j], columnsWeight[j],
				valuesWeight[j], equivCard[equiv]);
	}
	else if (m_weightType[i][j] == STIFFNESS) {
	  setStiffnessWeights(i, j, rowBeginWeight[j], columnsWeight[j],
			      valuesWeight[j], equivDofsGlobal);
	}
	else {
	  setDeluxeWeights(i, j, rowBeginWeight[j], columnsWeight[j],
			   valuesWeight[j], deluxeWeights[i][j]);
	}
	equivDofsGlobal += numEquivDofs;
      }
      m_Subdomain[i]->setEquivalenceClassWeightMatrices
	(rowBeginWeight, columnsWeight, valuesWeight);
    }
  }

  void setCardinalityWeights(LO sub, 
			     LO equiv, 
			     std::vector<LO> & rowBeginWeight, 
			     std::vector<LO> & columnsWeight,
			     std::vector<SX> & valuesWeight, 
			     LO equivCard)
  {
    LO numEquivDofs = m_Subdomain[sub]->getNumEquivDofs(equiv);
    rowBeginWeight.resize(numEquivDofs+1, 0);
    columnsWeight.resize(numEquivDofs);
    valuesWeight.resize(numEquivDofs);
    for (LO i=0; i<numEquivDofs; i++) {
      columnsWeight[i] = i;
      valuesWeight[i] = SX(1)/SM(equivCard);
      rowBeginWeight[i+1] = rowBeginWeight[i] + 1;
    }
  }

  void setStiffnessWeights(LO sub, 
			   LO equiv, 
			   std::vector<LO> & rowBeginWeight, 
			   std::vector<LO> & columnsWeight,
			   std::vector<SX> & valuesWeight,
			   const LO* equivDofsGlobal)
  {
    LO numEquivDofs = m_Subdomain[sub]->getNumEquivDofs(equiv);
    rowBeginWeight.resize(numEquivDofs+1, 0);
    columnsWeight.resize(numEquivDofs);
    valuesWeight.resize(numEquivDofs);
    std::vector<SM> diagValues;
    m_Subdomain[sub]->getBoundaryDiagValues(equiv, diagValues);
    assert (LO(diagValues.size()) == numEquivDofs);
    for (LO i=0; i<numEquivDofs; i++) {
      columnsWeight[i] = i;
      LO dofB = equivDofsGlobal[i];
      valuesWeight[i] = diagValues[i]/m_diagBoundary[dofB];
      rowBeginWeight[i+1] = rowBeginWeight[i] + 1;
    }
  }

  void setDeluxeWeights(LO sub, 
			LO equiv, 
			std::vector<LO> & rowBeginWeight, 
			std::vector<LO> & columnsWeight,
			std::vector<SX> & valuesWeight,
			std::vector<SX> & deluxeWeights)
  {
    LO numEquivDofs = m_Subdomain[sub]->getNumEquivDofs(equiv);
    LO numTerms = numEquivDofs*numEquivDofs;
    assert (LO(deluxeWeights.size()) == numTerms);
    rowBeginWeight.resize(numEquivDofs+1, 0);
    columnsWeight.resize(numTerms);
    valuesWeight.resize(numTerms);
    numTerms = 0;
    for (LO i=0; i<numEquivDofs; i++) {
      for (LO j=0; j<numEquivDofs; j++) {
	columnsWeight[numTerms] = j;
	valuesWeight[numTerms] = deluxeWeights[i+j*numEquivDofs];
	numTerms++;
      }
      rowBeginWeight[i+1] = numTerms;
    }
  }

  SM checkBDDCWeights()
  {
    size_t numProcDofB = m_dofMapB->getNodeNumElements();
    std::vector<SX> sumWeightsB(numProcDofB);
    for (LO i=0; i<m_numSub; i++) {
      LO numSubDofB = m_Subdomain[i]->getNumBoundaryDofs();
      std::vector<SX> x(numSubDofB, 1), Wx(numSubDofB);
      bool applyTranspose = false;
      m_Subdomain[i]->applyWeights(&x[0], applyTranspose, &Wx[0]);
      assert (m_subBoundaryDofs[i].size() == numSubDofB);
      for (LO j=0; j<numSubDofB; j++) {
	LO dofB = m_subBoundaryDofs[i][j];
	sumWeightsB[dofB] += Wx[j];
      }
    }
    MV SumWeightsB(m_dofMapB, 1);
    for (LO i=0; i<numProcDofB; i++) {
      SumWeightsB.replaceLocalValue(i, 0, sumWeightsB[i]);
    }
    MV SumWeightsB1to1(m_dofMapB1to1, 1);
    SumWeightsB1to1.doExport(SumWeightsB, *m_exporterB, Tpetra::ADD);
    Teuchos::ArrayRCP<const SX> values = SumWeightsB1to1.getData(0);
    SM tol(1e-10);
    SM maxError = 0;
    for (size_t i=0; i<m_dofMapB1to1->getNodeNumElements(); i++) {
      SM error = std::abs(values[i]-1);
      if (error > maxError) maxError = error;
    }
    SM maxErrorAll;
    Teuchos::reduceAll<int, SM>(*m_Comm, Teuchos::REDUCE_MAX, 1,
				&maxError, &maxErrorAll);
    return maxErrorAll;
  }

};

} // namespace bddc

#endif // WEIGHTSBDDC_H
  
