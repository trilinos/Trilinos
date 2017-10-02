
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

#ifndef SOLVERFACTORYBDDC_H
#define SOLVERFACTORYBDDC_H

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>
#include <assert.h>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ShyLUBDDC_config.h"
#include "shylu_SolverBaseBDDC.h"
#include "shylu_SolverLAPACK.h"

/*
#if defined(HAVE_SHYLUBDDC_TRILINOSSS)
#include "shylu_SolverKLU2.h"
#endif
*/

#if defined(HAVE_SHYLUBDDC_SHYLUTACHO)
#include "shylu_SolverTacho.h"
#endif

#if defined(HAVE_SHYLUBDDC_SUPERLU)
#include "shylu_SolverSuperLU.h"
#endif

#if defined(HAVE_SHYLUBDDC_PARDISO_MKL)
#include "shylu_SolverPardisoBDDC.h"
#endif

#if defined(USE_INTEL_CLUSTER_PARDISO)
#include "shylu_SolverClusterPardiso.h"
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
      SolverPtr = new SolverSuperLU<SX>(numRows,
					rowBegin,
					columns,
					values,
					Parameters);
    }
    else if (solverString == "Tacho") {
#if defined(HAVE_SHYLUBDDC_SHYLUTACHO)
      SolverPtr = new SolverTacho<SX>(numRows,
				      rowBegin,
				      columns,
				      values,
				      Parameters);
#else
      std::cout << "Error: Tacho solver is not available\n";
#endif
    }
    else if (solverString == "Pardiso") {
#if defined(HAVE_SHYLUBDDC_PARDISO_MKL)
      SolverPtr = new SolverPardiso<SX>(numRows, 
					rowBegin,
					columns,
					values,
					Parameters);
#else
      std::cout << "Error: Pardiso solver is not available\n";
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
      std::cout << "Error: Cluster Pardiso solver is not available\n";
#endif
    }
    else if (solverString == "LAPACK") {
      SolverPtr = new SolverLAPACK<SX>(numRows,
				       rowBegin,
				       columns,
				       values,
				       Parameters);
    }
    else {
      std::string msg("Error: no acceptable direct bddc solver.");
      msg += " Requested solver is ";
      msg += solverString;
      throw msg;
    }
    return SolverPtr;
  }
 private:
 protected:
};
 
} // namespace bddc

#endif // SOLVERFACTORYBDDC_H
