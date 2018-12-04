
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

#include "ShyLU_DDBDDC_config.h"

#if defined(HAVE_SHYLU_DDBDDC_SUPERLU)

#include "shylu_errorBDDC.hpp"
#include "shylu_SuperLU.hpp"

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
  
  SuperLU::SuperLU()
  {
  }
  
  SuperLU::~SuperLU()
  {
    if (m_numRows == 0) return;
    Destroy_SuperMatrix_Store(&m_A);
    Destroy_SuperNode_Matrix(&m_L);
    Destroy_CompCol_Matrix(&m_U);
  }
  
  void SuperLU::initialize(const int numRows)
  {
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
    m_etree.resize(numRows);
    m_perm_r.resize(numRows);
    m_perm_c.resize(numRows);
  }
  
  int SuperLU::initialize(int numRows,
			  int* rowBegin,
			  int* columns,
			  double* values)
  {
    m_numRows = numRows;
    if (numRows == 0) return 0;
    initialize(numRows);
    generateMatrices(numRows, rowBegin, columns, values);
    mem_usage_t mem_usage;
    factorMatrices(values[0], mem_usage);
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

  int SuperLU::initialize(int numRows,
			  int* rowBegin,
			  int* columns,
			  float* values)
  {
    m_numRows = numRows;
    if (numRows == 0) return 0;
    initialize(numRows);
    generateMatrices(numRows, rowBegin, columns, values);
    mem_usage_t mem_usage;
    factorMatrices(values[0], mem_usage);
    Destroy_SuperMatrix_Store(&m_X);
    Destroy_SuperMatrix_Store(&m_B);
    return 0;
}

  void SuperLU::solve(int NRHS,
		      double* Rhs, 
		      double* Sol)
  {
    int numRows = m_numRows;
    if (numRows == 0) return;
    m_ferr_d.resize(NRHS);
    m_berr_d.resize(NRHS);
    m_options.Fact = FACTORED; // already factored matrix
    SuperLUStat_t stat;
    StatInit(&stat);
    int info(0);
    double rpg(0), rcond(0);
    mem_usage_t mem_usage;
    solveEquations(numRows, NRHS, Rhs, Sol, stat, mem_usage, rpg, rcond, info);
    BDDC_TEST_FOR_EXCEPTION(info != 0, std::runtime_error, 
			    "SuperLU solveEquations error");
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&m_X);
    Destroy_SuperMatrix_Store(&m_B);
  }

  void SuperLU::solve(int NRHS,
		      float* Rhs, 
		      float* Sol)
  {
    int numRows = m_numRows;
    if (numRows == 0) return;
    m_ferr_f.resize(NRHS);
    m_berr_f.resize(NRHS);
    m_options.Fact = FACTORED; // already factored matrix
    SuperLUStat_t stat;
    StatInit(&stat);
    int info(0);
    float rpg(0), rcond(0);
    mem_usage_t mem_usage;
    solveEquations(numRows, NRHS, Rhs, Sol, stat, mem_usage, rpg, rcond, info);
    BDDC_TEST_FOR_EXCEPTION(info != 0, std::runtime_error, 
			    "SuperLU solveEquations error");
    StatFree(&stat);
    Destroy_SuperMatrix_Store(&m_X);
    Destroy_SuperMatrix_Store(&m_B);
  }

  void SuperLU::SuperLU::solveEquations(int numRows, 
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
	   m_equed, m_R_d.data(), m_C_d.data(), &m_L, &m_U, NULL, 0, &m_B, &m_X, 
	   &rpg, &rcond, m_ferr_d.data(), m_berr_d.data(), &mem_usage, &stat, 
	   &info);
  }

  void SuperLU::solveEquations(int numRows, 
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
	   m_equed, m_R_f.data(), m_C_f.data(), &m_L, &m_U, NULL, 0, &m_B, &m_X, 
	   &rpg, &rcond, m_ferr_f.data(), m_berr_f.data(), &mem_usage, &stat, 
	   &info);
  }

  void SuperLU::generateMatrices(int numRows, 
				 int* rowBegin, 
				 int* columns, 
				 double* values)
  {
    int nrhs = 1;
    int nnz = rowBegin[numRows];
    double *nullVec(0);
    dCreate_CompCol_Matrix(&m_A, numRows, numRows, nnz, values, columns,
			   rowBegin, SLU_NR, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&m_B, numRows, nrhs, nullVec, numRows, 
			 SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix(&m_X, numRows, nrhs, nullVec, numRows, 
			 SLU_DN, SLU_D, SLU_GE);
    m_R_d.resize(numRows);
    m_C_d.resize(numRows);
    m_ferr_d.resize(nrhs);
    m_berr_d.resize(nrhs);
  }

  void SuperLU::generateMatrices(int numRows, 
				 int* rowBegin, 
				 int* columns, 
				 float* values)
  {
    int nrhs = 1;
    int nnz = rowBegin[numRows];
    float *nullVec(0);
    sCreate_CompCol_Matrix(&m_A, numRows, numRows, nnz, values, columns,
			   rowBegin, SLU_NR, SLU_S, SLU_GE);
    sCreate_Dense_Matrix(&m_B, numRows, nrhs, nullVec, numRows, 
			 SLU_DN, SLU_S, SLU_GE);
    sCreate_Dense_Matrix(&m_X, numRows, nrhs, nullVec, numRows, 
			 SLU_DN, SLU_S, SLU_GE);
    m_R_f.resize(numRows);
    m_C_f.resize(numRows);
    m_ferr_f.resize(nrhs);
    m_berr_f.resize(nrhs);
  }

  void SuperLU::factorMatrices(double value,
			       mem_usage_t & mem_usage)
  {
    SuperLUStat_t stat;
    StatInit(&stat);
    m_B.ncol = 0; // signal to only do the factorization
    int info(0);
    double rpg(0), rcond(0);
    dgssvx(&m_options, &m_A, m_perm_c.data(), m_perm_r.data(), m_etree.data(), 
	   m_equed, m_R_d.data(), m_C_d.data(), &m_L, &m_U, NULL, 0, &m_B, &m_X, 
	   &rpg, &rcond, m_ferr_d.data(), m_berr_d.data(), &mem_usage, &stat, 
	   &info);
    BDDC_TEST_FOR_EXCEPTION(info != 0, std::runtime_error, 
			    "SuperLU dgssvx error");
    StatFree(&stat);
  }

  void SuperLU::factorMatrices(float value,
			       mem_usage_t & mem_usage)
  {
    SuperLUStat_t stat;
    StatInit(&stat);
    m_B.ncol = 0; // signal to only do the factorization
    int info(0);
    float rpg(0), rcond(0);
    sgssvx(&m_options, &m_A, m_perm_c.data(), m_perm_r.data(), m_etree.data(), 
	   m_equed, m_R_f.data(), m_C_f.data(), &m_L, &m_U, NULL, 0, &m_B, &m_X, 
	   &rpg, &rcond, m_ferr_f.data(), m_berr_f.data(), &mem_usage, &stat, 
	   &info);
    BDDC_TEST_FOR_EXCEPTION(info != 0, std::runtime_error, 
			    "SuperLU sgssvx error");
    StatFree(&stat);
  }
  
} // namespace bddc

#endif // HAVE_SHYLU_DDBDDC_SUPERLU
  
