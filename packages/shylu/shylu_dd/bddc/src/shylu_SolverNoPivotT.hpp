
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

#ifndef BDDC_SOLVERNOPIVOTT_H
#define BDDC_SOLVERNOPIVOTT_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "shylu_SolverBaseBDDC.hpp"
#include "NoPivotT.h"
#include <metis.h>

#if defined(_OPENMP) && defined(HAVE_SHYLU_DDBDDC_SHYLU_NODEHTS)
#include "shylu_hts.hpp"
#endif

namespace bddc {
  
template <class SX, class SM, class LO, class GO> class SolverNoPivotT : 
  public SolverBase<SX,SM,LO,GO>
{
public:
  SolverNoPivotT(LO numRows,
		LO* rowBegin,
		LO* columns,
		SX* values,
		Teuchos::ParameterList & Parameters) :
  SolverBase<SX,SM,LO,GO>(numRows, rowBegin, columns, values, Parameters),
    m_Solver(0)
  {
    m_useHtsForSolves = useHtsForSolves();
  }

  ~SolverNoPivotT()
  {
    delete m_Solver;
#if defined(_OPENMP) && defined(HAVE_SHYLU_DDBDDC_SHYLU_NODEHTS)
    if (m_useHtsForSolves) {
      HTST::delete_Impl(m_Uimpl);
      HTST::delete_Impl(m_Limpl);
    }
#endif
  }

  const std::vector<double> & getTimingsSuperFactorSerial() const
  {
    return m_Solver->getTimingsSuperFactorSerial();
  }

  const int* getPermutation() const
  {
    return m_Solver->getPermutation();
  }

  int Initialize()
  {
    LO numRows = this->m_numRows;
    if (numRows == 0) return 0;
    // get fill-reducing reordering
    LO* nullPtr(0);
    LO* perm = this->m_Parameters.get("Permutation Vector", nullPtr);
    if (perm == nullPtr) {
      Reorder(numRows, this->m_rowBegin, this->m_columns, m_permInput);
      perm = &m_permInput[0];
    }
    // set solver parameters
    bool useThreadsForFactorization = 
      this->m_Parameters.get("Use Threads For Factorization", false);
    bool useThreadsForTriangularSolves = 
      this->m_Parameters.get("Use Threads For Triangular Solves", false);
    int blockSizeForDenseTriangularSolves = 
      this->m_Parameters.get("Block Size For Dense Triangular Solves", 64);
    double minNumberTasksMultiplier =
      this->m_Parameters.get("Minimum Number Task Multiplier", double(1));
    bool useLargerMemoryUpdates = 
      this->m_Parameters.get("Use Larger Memory Updates", true);
    int threadOptionFactor =
      this->m_Parameters.get("Thread Option Factor", 0);
    int threadOptionLower =
      this->m_Parameters.get("Thread Option Lower", 0);
    int threadOptionUpper =
      this->m_Parameters.get("Thread Option Upper", 0);
    int numAdditionalTaskLevels = 
      this->m_Parameters.get("Number of Additional Task Levels", 0);
    m_Solver = new nopivot::NoPivotT<SX>();
    m_Solver->setThreadsModel
      (useThreadsForFactorization, useThreadsForTriangularSolves,
       threadOptionFactor, threadOptionLower, threadOptionUpper,
       minNumberTasksMultiplier, numAdditionalTaskLevels,
       blockSizeForDenseTriangularSolves, useLargerMemoryUpdates);
    m_Solver->SymbolicFactorization(this->m_numRows,
				    this->m_rowBegin,
				    this->m_columns,
				    perm);
    m_Solver->NumericalFactorization(this->m_values);
    initializeHts();
    return 0;
  }

  bool IsDirectSolver()
  {
    return true;
  }

  void MySolve(int NRHS,
	       SX* Rhs, 
	       SX* Sol)
  {
    if (this->m_numRows == 0) return;
    if (m_useHtsForSolves) {
#if defined(_OPENMP) && defined(HAVE_SHYLU_DDBDDC_SHYLU_NODEHTS)
      HTST::solve_omp(m_Limpl, Rhs, NRHS, Sol);
      HTST::solve_omp(m_Uimpl, Sol, NRHS);
#endif
    }
    else {
      int solveOption = 0; // not used
      m_Solver->Solve(NRHS, Rhs, Sol, solveOption);
    }
  }

  bool MyExactSolver() 
  {
    return true;
  }

private:
  std::vector<idx_t> m_permInput;
  nopivot::NoPivotT<SX>* m_Solver;
  bool m_useHtsForSolves;
#if defined(_OPENMP) && defined(HAVE_SHYLU_DDBDDC_SHYLU_NODEHTS)
  typedef Experimental::HTS<LO, LO, SM> HTST;
  typename HTST::Impl *m_Limpl, *m_Uimpl;
#endif

  bool useHtsForSolves() 
  {
#if defined(_OPENMP) && defined(HAVE_SHYLU_DDBDDC_SHYLU_NODEHTS)
    if (this->m_Parameters.get("Use HTS For Solves", false)) {
      return true;
    }
    else {
      return false;
    }
#else
    return false;
#endif
  }

  void initializeHts()
  {
    if ( ! m_useHtsForSolves) return;
#if defined(_OPENMP) && defined(HAVE_SHYLU_DDBDDC_SHYLU_NODEHTS)
    std::vector<int> ir, jc;
    std::vector<SX> d;
    std::vector<SX> diag;
    int in_nnz, out_nnz;
    m_Solver->extractLAsCcsMatrix(jc, ir, d, in_nnz, out_nnz);
    m_Solver->extractDiag(diag);
    const int nthreads = omp_get_max_threads();
    {
      typename HTST::CrsMatrix* L = HTST::make_CrsMatrix
        (ir.size() - 1, ir.data(), jc.data(), d.data(), true);
      m_Limpl = HTST::preprocess
        (L, 1, nthreads, true, m_Solver->getPermutation(), NULL, NULL);
      HTST::delete_CrsMatrix(L);
    }
    {
      typename HTST::CrsMatrix* U = HTST::make_CrsMatrix
        (ir.size() - 1, ir.data(), jc.data(), d.data());
      m_Uimpl = HTST::preprocess
        (U, 1, nthreads, true, NULL, m_Solver->getPermutation(), diag.data());
      HTST::delete_CrsMatrix(U);
    }
    delete m_Solver; m_Solver = 0;
#endif
  }
  
  void determineComponents(LO *A1, 
			   LO *A2, 
			   LO N, 
			   std::vector<LO> & component)
  {
    LO i;
    component.resize(N);
    component.assign(N, -1);
    componentsFunction(A1, A2, N, &component[0]);
  }

  void componentsFunction(LO *A1, 
			  LO *A2, 
			  LO N, 
			  LO* component)
  {
    LO i, comp_num;
    comp_num = 0;
    for (i=0; i<N; i++) {
      if (component[i] == -1) {
	depthFirstSearch(i, comp_num, component, A1, A2);
	comp_num++;
      }
    }
  }

  void depthFirstSearch(const LO v, 
			const LO comp_num, 
			LO* component,
			LO* A1, 
			LO* A2)
  {
    LO i, adj_vertex;
    component[v] = comp_num;
    for (i=A2[v]; i<A2[v+1]; i++) {
      adj_vertex = A1[i];
      if (component[adj_vertex] == -1) 
	depthFirstSearch(adj_vertex, comp_num, component, A1, A2);
    }
  }

  LO determineNumComponents(std::vector<LO> & component)
  {
    LO maxComp(0);
    for (size_t i=0; i<component.size(); i++) {
      if (component[i] > maxComp) maxComp = component[i];
    }
    return maxComp+1;
  }

  void Reorder(LO numRows,
	       LO* rowBegin,
	       LO* columns,
	       std::vector<idx_t> & perm)
  {
    // Metis parameters
    std::vector<idx_t> options(METIS_NOPTIONS, 0);
    METIS_SetDefaultOptions(&options[0]);
    //  options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_EDGE;
    options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_NODE;
    const int cStyleNumbering = 0;
    options[METIS_OPTION_NUMBERING] = cStyleNumbering;
    // determine if graph has disconnected parts
    std::vector<LO> component;
    determineComponents(columns, rowBegin, numRows, component);
    LO numComponents = determineNumComponents(component);
    std::vector<idx_t> rowBeginGraph, columnsGraph;
    perm.resize(numRows);
    if (numComponents == 1) {
      // graph used by Metis should not include diagonal entries
      determineMetisGraph(numRows, rowBegin, columns, rowBeginGraph,
			  columnsGraph);
      std::vector<idx_t> invPerm(numRows);
      METIS_NodeND(&numRows, &rowBeginGraph[0], &columnsGraph[0], 0, 
		   &options[0], &perm[0], &invPerm[0]);
    }
    else {
      std::vector<LO> imap(numRows, -1), activeRows(numRows);
      LO numRowsGraph, sumNumRowsGraph(0);
      for (LO i=0; i<numComponents; i++) {
	determineSubGraph(i, rowBegin, columns, component, imap, activeRows,
			  numRowsGraph, rowBeginGraph, columnsGraph);
	std::vector<idx_t> permGraph(numRowsGraph), invPermGraph(numRowsGraph);
	METIS_NodeND(&numRowsGraph, &rowBeginGraph[0], &columnsGraph[0], 0, 
		     &options[0], &permGraph[0], &invPermGraph[0]);
	for (LO j=0; j<numRowsGraph; j++) {
	  perm[activeRows[j]] = permGraph[j] + sumNumRowsGraph;
	}
	sumNumRowsGraph += numRowsGraph;
      }
    }
  }

  void determineSubGraph(LO comp, 
			 LO* rowBegin, 
			 LO* columns, 
			 std::vector<LO> & component, 
			 std::vector<LO> & imap, 
			 std::vector<LO> & activeRows,
			 LO & numRowsGraph, 
			 std::vector<idx_t> & rowBeginGraph, 
			 std::vector<idx_t> & columnsGraph)
  {
    LO numRows = component.size();
    numRowsGraph = 0;
    for (LO i=0; i<numRows; i++) {
      if (component[i] == comp) {
	activeRows[numRowsGraph] = i;
	imap[i] = numRowsGraph++;
      }
    }
    rowBeginGraph.resize(numRowsGraph+1, 0);
    LO numTerms = 0;
    for (LO i=0; i<numRowsGraph; i++) {
      LO row = activeRows[i];
      for (LO j=rowBegin[row]; j<rowBegin[row+1]; j++) {
	LO col = imap[columns[j]];
	if (col != -1) {
	  if (col != i) numTerms++;
	}
      }
      rowBeginGraph[i+1] = numTerms;
    }
    columnsGraph.resize(numTerms);
    numTerms = 0;
    for (LO i=0; i<numRowsGraph; i++) {
      LO row = activeRows[i];
      for (LO j=rowBegin[row]; j<rowBegin[row+1]; j++) {
	LO col = imap[columns[j]];
	if (col != -1) {
	  if (col != i) columnsGraph[numTerms++] = col;
	}
      }
    }
    for (LO i=0; i<numRowsGraph; i++) imap[activeRows[i]] = -1;
  }
    
  void determineMetisGraph(LO numRows, 
			   const LO* rowBegin, 
			   const LO* columns, 
			   std::vector<idx_t> & rowBeginGraph,
			   std::vector<idx_t> & columnsGraph)
  {
    rowBeginGraph.resize(numRows+1, 0);
    LO numTerms = rowBegin[numRows];
    columnsGraph.resize(numTerms);
    numTerms = 0;
    for (LO i=0; i<numRows; i++) {
      for (LO j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	LO col = columns[j];
	if (col != i) {
	  columnsGraph[numTerms] = col;
	  numTerms++;
	}
      }
      rowBeginGraph[i+1] = numTerms;
    }
  }

 };
  
} // namespace bddc

#endif // BDDC_SOLVERNOPIVOTT_H
  
