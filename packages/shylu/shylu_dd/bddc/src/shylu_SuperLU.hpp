
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

#ifndef BDDC_SUPERLU_HPP
#define BDDC_SUPERLU_HPP
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <vector>

#include "slu_ddefs.h"
#include "supermatrix.h"
#include "slu_util.h"

namespace bddc {
  
class SuperLU
{
public:
  SuperLU();
  ~SuperLU();
  int initialize(int numRows,
		 int* rowBegin,
		 int* columns,
		 double* values);
  int initialize(int numRows,
		 int* rowBegin,
		 int* columns,
		 float* values);
  void solve(int NRHS,
	     double* Rhs, 
	     double* Sol);

  void solve(int NRHS,
	     float* Rhs, 
	     float* Sol);
private:
  void initialize(const int numRows);
  void solveEquations(int numRows, 
		      int NRHS, 
		      double* Rhs, 
		      double* Sol, 
		      SuperLUStat_t & stat,
		      mem_usage_t & mem_usage,
		      double & rpg, 
		      double & rcond, 
		      int & info);
  void solveEquations(int numRows, 
		      int NRHS, 
		      float* Rhs, 
		      float* Sol, 
		      SuperLUStat_t & stat,
		      mem_usage_t & mem_usage,
		      float & rpg, 
		      float & rcond, 
		      int & info);
  void generateMatrices(int numRows, 
			int* rowBegin, 
			int* columns, 
			double* values);
  void generateMatrices(int numRows, 
			int* rowBegin, 
			int* columns, 
			float* values);
  void factorMatrices(double value,
		      mem_usage_t & mem_usage);
  void factorMatrices(float value,
		      mem_usage_t & mem_usage);

private:
  int m_numRows{0};
  SuperMatrix m_A, m_L, m_U, m_X, m_B;
  superlu_options_t m_options;
  std::vector<int> m_etree, m_perm_r, m_perm_c;
  std::vector<double> m_R_d, m_C_d, m_ferr_d, m_berr_d;
  std::vector<float> m_R_f, m_C_f, m_ferr_f, m_berr_f;
  char m_equed[1];

};
  
} // namespace bddc
  
#endif
