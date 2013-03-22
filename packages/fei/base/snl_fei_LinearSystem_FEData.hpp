/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
