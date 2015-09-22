/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
