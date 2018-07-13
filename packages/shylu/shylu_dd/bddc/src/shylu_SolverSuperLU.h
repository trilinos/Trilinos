
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

#ifndef SOLVERSUPERLUBDDC_H
#define SOLVERSUPERLUBDDC_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "shylu_SolverBaseBDDC.h"

typedef int int_t;
#include "supermatrix.h"
#include "slu_util.h"

extern "C"{
extern void
dgssvx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, int, SuperMatrix *, SuperMatrix *,
       double *, double *, double *, double *,
       mem_usage_t *, SuperLUStat_t *, int *);
extern void
dCreate_Dense_Matrix(SuperMatrix *, int, int, double *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_CompCol_Matrix(SuperMatrix *, int, int, int, double *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
sgssvx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, float *, float *, SuperMatrix *, SuperMatrix *,
       void *, int, SuperMatrix *, SuperMatrix *,
       float *, float *, float *, float *,
       mem_usage_t *, SuperLUStat_t *, int *);
extern void
sCreate_Dense_Matrix(SuperMatrix *, int, int, float *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
sCreate_CompCol_Matrix(SuperMatrix *, int, int, int, float *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
}

namespace bddc {
  
template <class SX> class SolverSuperLU : 
  public SolverBase<SX>
{
public:
  SolverSuperLU(int numRows,
		int* rowBegin,
		int* columns,
		SX* values,
		Teuchos::ParameterList & Parameters) :
  SolverBase<SX>(numRows, rowBegin, columns, values, Parameters)
  {
  }

  ~SolverSuperLU()
  {
    Destroy_SuperMatrix_Store(&m_A);
    Destroy_SuperNode_Matrix(&m_L);
    Destroy_CompCol_Matrix(&m_U);
  }

  int Initialize()
  {
    int numRows = this->m_numRows;
    if (numRows == 0) return 0;
    set_default_options(&m_options);
    m_options.Fact = DOFACT;
    m_options.Equil = NO;
    m_options.ColPerm = MMD_AT_PLUS_A;
    m_options.Trans = NOTRANS;
    m_options.IterRefine = NOREFINE;
    m_options.DiagPivotThresh = 0.001;
    m_options.SymmetricMode = YES;
    m_options.PivotGrowth = NO;
    m_options.ConditionNumber = NO;
    m_options.PrintStat = NO;
    int* rowBegin = this->m_rowBegin;
    int* columns = this->m_columns;
    SX* values =  this->m_values;
    std::vector<SX> rhsx(numRows), rhsb(numRows);
    int nrhs(1);
    generateMatrices(numRows, nrhs, rowBegin, columns, values);
    m_etree.resize(numRows);
    m_perm_r.resize(numRows);
    m_perm_c.resize(numRows);
    m_R.resize(numRows);
    m_C.resize(numRows);
    m_ferr.resize(nrhs);
    m_berr.resize(nrhs);
    SuperLUStat_t stat;
    StatInit(&stat);
    m_B.ncol = 0; // signal to only do the factorization
    mem_usage_t mem_usage;
    factorMatrices(stat, mem_usage, values[0]);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&m_X);
    Destroy_SuperMatrix_Store(&m_B);
    /*
    std::cout << "SuperLU mem_usage.for_lu MB = " << mem_usage.for_lu/1e6 
	      << std::endl;
    std::cout << "SuperLU mem_usage.total_needed MB = "
	      << mem_usage.total_needed/1e6 << std::endl;
    */
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
    int numRows = this->m_numRows;
    if (numRows == 0) return;
    m_ferr.resize(NRHS);
    m_berr.resize(NRHS);
    m_options.Fact = FACTORED; // already factored matrix
    SuperLUStat_t stat;
    StatInit(&stat);
    int info(0);
    SX rpg(0), rcond(0);
    mem_usage_t mem_usage;
    solveEquations(numRows, NRHS, Rhs, Sol, stat, mem_usage, rpg, rcond, info);
    assert (info == 0);
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&m_X);
    Destroy_SuperMatrix_Store(&m_B);
  }

  bool MyExactSolver() 
  {
    return true;
  }

private:
  SuperMatrix m_A, m_L, m_U, m_X, m_B;
  superlu_options_t m_options;
  std::vector<int> m_etree, m_perm_r, m_perm_c;
  std::vector<SX> m_R, m_C, m_ferr, m_berr;
  char m_equed[1];

  void solveEquations(int numRows, 
		      int NRHS, 
		      double* Rhs, 
		      double* Sol, 
		      SuperLUStat_t & stat,
		      mem_usage_t & mem_usage,
		      double & rpg, 
		      double & rcond, 
		      int & info)
  {
    dCreate_Dense_Matrix(&m_B, numRows, NRHS, Rhs, numRows, 
			 SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&m_X, numRows, NRHS, Sol, numRows, 
			 SLU_DN, SLU_D, SLU_GE);
    dgssvx(&m_options, &m_A, m_perm_c.data(), m_perm_r.data(), m_etree.data(), 
	   m_equed, m_R.data(), m_C.data(), &m_L, &m_U, NULL, 0, &m_B, &m_X, 
	   &rpg, &rcond, m_ferr.data(), m_berr.data(), &mem_usage, &stat, 
	   &info);
  }

  void solveEquations(int numRows, 
		      int NRHS, 
		      float* Rhs, 
		      float* Sol, 
		      SuperLUStat_t & stat,
		      mem_usage_t & mem_usage,
		      float & rpg, 
		      float & rcond, 
		      int & info)
  {
    sCreate_Dense_Matrix(&m_B, numRows, NRHS, Rhs, numRows, 
			 SLU_DN, SLU_S, SLU_GE);
    sCreate_Dense_Matrix(&m_X, numRows, NRHS, Sol, numRows, 
			 SLU_DN, SLU_S, SLU_GE);
    sgssvx(&m_options, &m_A, m_perm_c.data(), m_perm_r.data(), m_etree.data(), 
	   m_equed, m_R.data(), m_C.data(), &m_L, &m_U, NULL, 0, &m_B, &m_X, 
	   &rpg, &rcond, m_ferr.data(), m_berr.data(), &mem_usage, &stat, 
	   &info);
  }

  void generateMatrices(int numRows, 
			int nrhs,
			int* rowBegin, 
			int* columns, 
			double* values)
  {
    int nnz = rowBegin[numRows];
    double *nullVec(0);
    dCreate_CompCol_Matrix(&m_A, numRows, numRows, nnz, values, columns,
			   rowBegin, SLU_NR, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&m_B, numRows, nrhs, nullVec, numRows, 
			 SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&m_X, numRows, nrhs, nullVec, numRows, 
			 SLU_DN, SLU_D, SLU_GE);
  }

  void generateMatrices(int numRows, 
			int nrhs,
			int* rowBegin, 
			int* columns, 
			float* values)
  {
    int nnz = rowBegin[numRows];
    float *nullVec(0);
    sCreate_CompCol_Matrix(&m_A, numRows, numRows, nnz, values, columns,
			   rowBegin, SLU_NR, SLU_S, SLU_GE);
    sCreate_Dense_Matrix(&m_B, numRows, nrhs, nullVec, numRows, 
			 SLU_DN, SLU_S, SLU_GE);
    sCreate_Dense_Matrix(&m_X, numRows, nrhs, nullVec, numRows, 
			 SLU_DN, SLU_S, SLU_GE);
  }

  void factorMatrices(SuperLUStat_t & stat,
		      mem_usage_t & mem_usage,
		      double value)
  {
    int info(0);
    double rpg(0), rcond(0);
    dgssvx(&m_options, &m_A, m_perm_c.data(), m_perm_r.data(), m_etree.data(), 
	   m_equed, m_R.data(), m_C.data(), &m_L, &m_U, NULL, 0, &m_B, &m_X, 
	   &rpg, &rcond, m_ferr.data(), m_berr.data(), &mem_usage, &stat, 
	   &info);
    assert (info == 0);
  }

  void factorMatrices(SuperLUStat_t & stat,
		      mem_usage_t & mem_usage,
		      float value)
  {
    int info(0);
    float rpg(0), rcond(0);
    sgssvx(&m_options, &m_A, m_perm_c.data(), m_perm_r.data(), m_etree.data(), 
	   m_equed, m_R.data(), m_C.data(), &m_L, &m_U, NULL, 0, &m_B, &m_X, 
	   &rpg, &rcond, m_ferr.data(), m_berr.data(), &mem_usage, &stat, 
	   &info);
    assert (info == 0);
  }

 };
  
} // namespace bddc

#endif // SOLVERSUPERLUBDDC_H
  
