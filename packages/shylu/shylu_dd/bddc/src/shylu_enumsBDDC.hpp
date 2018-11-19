
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

#ifndef BDDC_ENUMS_H
#define BDDC_ENUMS_H
#include <stdio.h>
#include <iostream>
#include <fstream>

namespace bddc {

enum ProblemType {SCALARPDE, ELASTICITY, HELMHOLTZ};

enum EquivType{CORNER, EDGE, FACE};

enum WeightType{STIFFNESS, CARDINALITY, DELUXE};

enum MatrixType{OPERATOR, PRECONDITIONER};

enum AnalysisType {STANDARD = 0, HELMHOLTZA = 1};

enum Timings{
  TIME_INIT_STATIC_COND, TIME_STATIC_EXPANSION, TIME_SUB_CORR, 
  TIME_COARSE_CORR, TIME_APPLY_OPER, TIME_INITIALIZATION, 
  TIME_APPLY_PRECONDITIONER, TIME_INIT_STEP1, TIME_DOF_MANAGER,
  TIME_PARTITION_OF_UNITY, TIME_COARSE_SPACE_PREP, TIME_VCS_REDUCTION,
  TIME_VCS_EXPORT, TIME_VCS_SOLVE, TIME_VCS_IMPORT, TIME_VCS_EXPANSION,
  TIME_INTERIOR_FACTORIZATIONS, TIME_BASE_CONSTRAINTS, 
  TIME_INIT_CORRECTION, TIME_LENGTH
};

} // namespace bddc

#endif // BDDC_ENUMS_H
