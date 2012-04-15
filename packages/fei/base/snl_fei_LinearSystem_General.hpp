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


#ifndef _snl_fei_LinearSystem_General_hpp_
#define _snl_fei_LinearSystem_General_hpp_

#include <fei_macros.hpp>
#include <fei_mpi.h>
#include <fei_CSVec.hpp>
#include <fei_LinearSystem.hpp>
#include <fei_Matrix.hpp>
#include <fei_Vector.hpp>
#include <fei_fwd.hpp>
#include <fei_Logger.hpp>

namespace fei {
  class DirichletBCManager;
}

namespace snl_fei {
  /** implementation of fei::LinearSystem interface */
  class LinearSystem_General : public fei::LinearSystem,
                               private fei::Logger {
  public:
    /** constructor */
    LinearSystem_General(fei::SharedPtr<fei::MatrixGraph>& matrixGraph);

    /** denstructor */
    virtual ~LinearSystem_General();

    /** Essential (dirichlet) boundary-condition function.
    */
    int loadEssentialBCs(int numIDs,
                         const int* IDs,
                         int idType,
                         int fieldID,
                         int offsetIntoField,
                         const double* prescribedValues);

    /** Essential (dirichlet) boundary-condition function.
    */
    int loadEssentialBCs(int numIDs,
                         const int* IDs,
                         int idType,
                         int fieldID,
                         const int* offsetIntoField,
                         const double* prescribedValues);

    /** load lagrange-multiplier constraint coefficients */
    int loadLagrangeConstraint(int constraintID,
			       const double *weights,
			       double rhsValue);

    /** load penalty constraint coefficients */
    int loadPenaltyConstraint(int constraintID,
			      const double *weights,
			      double penaltyValue,
			      double rhsValue);

    /** Signal that all boundary-conditions and constraint coefficients have
	been loaded, and they may now be applied to the linear system.
    */
    int loadComplete(bool applyBCs=true,
                     bool globalAssemble=true);

    /** Set parameters on this object. Currently two parameters are recognized:
	"debugOutput 'path'" where 'path' is the path to the location where
	debug-log files will be produced.<br>
	"name 'string'" where 'string' is an identifier that will be used in
	debug-log file-names.
    */
    int parameters(int numParams,
		   const char* const* paramStrings);

    /** parameters implementation */
    int parameters(const fei::ParameterSet& params);

    /** use stored BC values to modify specified vector */
    int setBCValuesOnVector(fei::Vector* vector);

    /** query whether specified eqn index has prescribed BC value */
    bool eqnIsEssentialBC(int globalEqnIndex) const;

    /** Retrieve eqn-indices and values for BC equations */
    void getEssentialBCs(std::vector<int>& bcEqns,
                         std::vector<double>& bcVals) const;

    /** Retrieve eqn-indices for constraints */
    void getConstrainedEqns(std::vector<int>& crEqns) const;

  private:
    void setName(const char* name);

    int fill_EssBCValues();

    int implementBCs(bool applyBCs);

    int enforceEssentialBC_LinSysCore();

    void enforceEssentialBC_step_1(fei::CSVec& essBCs);

    void enforceEssentialBC_step_2(fei::CSVec& essBCs);

    int getMatrixRow(fei::Matrix* matrix, int row,
		     std::vector<double>& coefs,
		     std::vector<int>& indices);

    MPI_Comm comm_;

    fei::CSVec* essBCvalues_;
    fei::CSVec* allEssBCs_;

    bool resolveConflictRequested_;
    bool bcs_trump_slaves_;
    bool explicitBCenforcement_;
    bool BCenforcement_no_column_mod_;

    int localProc_;
    int numProcs_;

    int firstLocalOffset_;
    int lastLocalOffset_;

    std::string name_;
    std::map<std::string, unsigned> named_loadcomplete_counter_;

    std::vector<int> iwork_;
    std::vector<double> dwork_;
    std::string dbgprefix_;
  };//class LinearSystem_General
}//namespace snl_fei

#endif // _snl_fei_LinearSystem_General_hpp_
