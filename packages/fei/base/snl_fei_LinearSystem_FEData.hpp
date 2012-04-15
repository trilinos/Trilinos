/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _snl_fei_LinearSystem_FEData_hpp_
#define _snl_fei_LinearSystem_FEData_hpp_

#include <fei_macros.hpp>
#include <fei_mpi.h>
#include <fei_utils.hpp>
#include <fei_LinearSystem.hpp>
#include <fei_Vector.hpp>
#include <fei_Matrix.hpp>
#include <fei_fwd.hpp>

namespace fei {
  class DirichletBCManager;
}

namespace snl_fei {
  /** implementation of fei::LinearSystem specialized for
     FiniteElementData */
  class LinearSystem_FEData : public fei::LinearSystem {
  public:
    /** constructor */
    LinearSystem_FEData(fei::SharedPtr<FiniteElementData>& fedata,
			fei::SharedPtr<fei::MatrixGraph>& matrixGraph);

    /** destructor */
    virtual ~LinearSystem_FEData();

    /** implementation of loadLagrangeConstraint */
    int loadLagrangeConstraint(int constraintID,
			       const double *weights,
			       double rhsValue);

    /** implementation of loadPenaltyConstraint */
    int loadPenaltyConstraint(int constraintID,
			      const double *weights,
			      double penaltyValue,
			      double rhsValue);

    /** Signal that all boundary-conditions and constraint coefficients have
	been loaded, and they may now be applied to the linear system.
    */
    int loadComplete(bool applyBCs=true,
                     bool globalAssemble=true);

    /** Retrieve FiniteElementData object */
    fei::SharedPtr<FiniteElementData> getFiniteElementData() { return( feData_ ); }

    /** Set parameters on this object. Currently two parameters are recognized:
	"debugOutput 'path'" where 'path' is the path to the location where
	debug-log files will be produced.<br>
	"name 'string'" where 'string' is an identifier that will be used in
	debug-log file-names.
    */
    int parameters(int numParams,
		   const char* const* paramStrings)
      { return( feData_->parameters(numParams, (char**)paramStrings) ); }

    /** implementation of parameters */
    int parameters(const fei::ParameterSet& params)
      {
	int numParams = 0;
	const char** paramStrings = NULL;
	std::vector<std::string> stdstrings;
	fei::utils::convert_ParameterSet_to_strings(&params, stdstrings);
	fei::utils::strings_to_char_ptrs(stdstrings, numParams, paramStrings);

	int err = parameters(numParams, paramStrings);

	delete [] paramStrings;

	return(err);
      }

    /** set previously specified BC values on given vector */
    int setBCValuesOnVector(fei::Vector* vector);

    /** set lookup object */
    void setLookup(Lookup* lookup)
      { lookup_ = lookup; }

    /** Query whether specified eqn has prescribed BC value. */
    bool eqnIsEssentialBC(int globalEqnIndex) const;

    /** Retrieve BC eqn indices. */
    void getEssentialBCs(std::vector<int>& bcEqns,
                         std::vector<double>& bcVals) const;

    /** Retrieve constrained eqn indices */
    void getConstrainedEqns(std::vector<int>& crEqns) const;

  private:
    int implementBCs(bool applyBCs);

    MPI_Comm comm_;
    int localProc_;
    int numProcs_;
    fei::SharedPtr<fei::Matrix> matrix_;
    fei::SharedPtr<fei::Vector> soln_;
    fei::SharedPtr<fei::Vector> rhs_;
    fei::SharedPtr<FiniteElementData> feData_;
    Lookup* lookup_;

    std::vector<char*> attributeNames_;
    std::vector<void*> attributes_;
  };//class LinearSystem_FEData
}//namespace snl_fei

#endif // _snl_fei_LinearSystem_FEData_hpp_
