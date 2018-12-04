
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

#ifndef BDDC_SOLVERFACTORY_H
#define BDDC_SOLVERFACTORY_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ShyLU_DDBDDC_config.h"
#include "shylu_SolverBaseBDDC.hpp"
#include "shylu_SolverLAPACK.hpp"
#include "shylu_errorBDDC.hpp"

#if defined(HAVE_SHYLU_DDBDDC_AMESOS2)
#include "shylu_SolverKLU2.hpp"
#endif

#if defined(HAVE_SHYLU_DDBDDC_SHYLU_NODETACHO)
#include "shylu_SolverTacho.hpp"
#endif

#if defined(HAVE_SHYLU_DDBDDC_SUPERLU)
#include "shylu_SolverSuperLU.hpp"
#endif

#if defined(HAVE_SHYLU_DDBDDC_PARDISO_MKL)
#include "shylu_SolverPardisoBDDC.hpp"
#endif

#if defined(USE_INTEL_CLUSTER_PARDISO)
#include "shylu_SolverClusterPardiso.hpp"
#endif

#if defined(HAVE_SHYLU_DDBDDC_MUELU)
#include "shylu_SolverMueLu.hpp"
#endif

//#define NODALAMGBDDC
#ifdef NODALAMGBDDC
#include "shylu_SolverNodalAMG.hpp"
#endif

namespace bddc {

template <class SX> class SolverFactory 
{
 public: // functions
  SolverFactory() { };

  ~SolverFactory() { };

  SolverBase<SX>* Generate(int numRows,
			   int* rowBegin,
			   int* columns,
			   SX* values,
			   Teuchos::ParameterList & Parameters,
			   MPI_Comm* pComm = 0)
  {
    SolverBase<SX>* SolverPtr = NULL;
    std::string solverString = Parameters.get("Solver", "SuperLU");
    if (solverString == "SuperLU") {
#if defined(HAVE_SHYLU_DDBDDC_SUPERLU)
      SolverPtr = new SolverSuperLU<SX>(numRows,
					rowBegin,
					columns,
					values,
					Parameters);
#else
      BDDC_TEST_FOR_EXCEPTION(1, std::runtime_error, 
			      "SuperLU solver is not available");
#endif
    }
    else if (solverString == "MueLu") {
#if defined(HAVE_SHYLU_DDBDDC_MUELU)
      SolverPtr = new SolverMueLu<SX>(numRows,
				      rowBegin,
				      columns,
				      values,
				      Parameters);
#else
      BDDC_TEST_FOR_EXCEPTION(1, std::runtime_error, 
			      "MueLu preconditioner is not available");
#endif
    }
    else if (solverString == "KLU2") {
#if defined(HAVE_SHYLU_DDBDDC_AMESOS2)
      SolverPtr = new SolverKLU2<SX>(numRows,
				     rowBegin,
				     columns,
				     values,
				     Parameters);
#else
      BDDC_TEST_FOR_EXCEPTION(1, std::runtime_error, 
			      "KLU2 solver is not available");
#endif
    }
    else if (solverString == "Tacho") {
#if defined(HAVE_SHYLU_DDBDDC_SHYLU_NODETACHO)
      SolverPtr = new SolverTacho<SX>(numRows,
				      rowBegin,
				      columns,
				      values,
				      Parameters);
#else
      BDDC_TEST_FOR_EXCEPTION(1, std::runtime_error, 
			      "Tacho solver is not available");
#endif
    }
    else if (solverString == "Pardiso") {
#if defined(HAVE_SHYLU_DDBDDC_PARDISO_MKL)
      SolverPtr = new SolverPardiso<SX>(numRows, 
					rowBegin,
					columns,
					values,
					Parameters);
#else
      BDDC_TEST_FOR_EXCEPTION(1, std::runtime_error, 
			      "Pardiso solver is not available");
#endif
    }
    else if (solverString == "Cluster Pardiso") {
#if defined(USE_INTEL_CLUSTER_PARDISO)
      SolverPtr = new SolverClusterPardiso<SX>(numRows, 
					       rowBegin,
					       columns,
					       values,
					       Parameters,
					       pComm);
#else
      BDDC_TEST_FOR_EXCEPTION(1, std::runtime_error, 
			      "Cluster Sparse solver is not available");
#endif
    }
    else if (solverString == "LAPACK") {
      SolverPtr = new SolverLAPACK<SX>(numRows,
				       rowBegin,
				       columns,
				       values,
				       Parameters);
    }
    else if (solverString == "NodalAMG") {
#if defined(NODALAMGBDDC)
      SolverPtr = new SolverNodalAMG<SX>(numRows,
					 rowBegin,
					 columns,
					 values,
					 Parameters);
#else
      BDDC_TEST_FOR_EXCEPTION(1, std::runtime_error, 
			      "NodalAMG preconditioner is not available");
#endif
    }
    else {
      std::string text = "Requested solver is " + solverString;
      std::cout << text << std::endl;
      BDDC_TEST_FOR_EXCEPTION(1, std::runtime_error, 
			      "solver not found");
    }
    return SolverPtr;
  }
 private:
 protected:
};
 
} // namespace bddc

#endif // BDDC_SOLVERFACTORY_H
