/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_LinearSystem_hpp_
#define _fei_LinearSystem_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Matrix.hpp>
#include <fei_Vector.hpp>
#include <fei_DirichletBCManager.hpp>

namespace fei {
  class Factory;
  class ParameterSet;

  /** A simple container to bind a matrix and two vectors together as the
      matrix, rhs and solution of a linear system.
  */
  class LinearSystem {
  public:
    /** LinearSystem Factory interface */
    class Factory {
    public:
      /** Usual virtual destructor */
      virtual ~Factory(){}

      /** Produce an instance of a LinearSystem. */
      virtual fei::SharedPtr<fei::LinearSystem>
	createLinearSystem(fei::SharedPtr<fei::MatrixGraph>& matrixGraph);
    };

    /** Constructor */
    LinearSystem(fei::SharedPtr<fei::MatrixGraph>& matrixGraph);

    /** Destructor */
    virtual ~LinearSystem();

    /** Set parameters on this object. Currently two parameters are recognized:
	"debugOutput 'path'" where 'path' is the path to the location where
	debug-log files will be produced.<br>
	"name 'string'" where 'string' is an identifier that will be used in
	debug-log file-names.
    */
    virtual int parameters(int numParams,
                           const char* const* paramStrings) = 0;

    /** Set parameters on this object.
     */
    virtual int parameters(const fei::ParameterSet& params) = 0;

    /** Set the matrix for this linear system. */
    virtual void setMatrix(fei::SharedPtr<fei::Matrix>& matrix);

    /** Get the matrix for this linear system. */
    virtual fei::SharedPtr<fei::Matrix> getMatrix()
      { return(matrix_); }

    /** Get the matrix for this linear system. */
    virtual fei::SharedPtr<const fei::Matrix> getMatrix() const
      { fei::SharedPtr<const fei::Matrix> const_mat(matrix_);
        return(const_mat); }

    /** Set the right-hand-side for this linear system. */
    virtual void setRHS(fei::SharedPtr<fei::Vector>& rhs)
      { rhs_ = rhs; }

    /** Get the right-hand-side for this linear system. */
    virtual fei::SharedPtr<fei::Vector> getRHS()
      { return(rhs_); }

    /** Get the right-hand-side for this linear system. */
    virtual fei::SharedPtr<const fei::Vector> getRHS() const
      { fei::SharedPtr<const fei::Vector> const_rhs(rhs_);
        return(const_rhs); }

    /** Set the solution for this linear system. */
    virtual void setSolutionVector(fei::SharedPtr<fei::Vector>& soln)
      { soln_ = soln; }

    /** Get the solution for this linear system. */
    virtual fei::SharedPtr<fei::Vector> getSolutionVector()
      { return(soln_); }

    /** Get the solution for this linear system. */
    virtual fei::SharedPtr<const fei::Vector> getSolutionVector() const
      { fei::SharedPtr<const fei::Vector> const_soln(soln_);
      return(const_soln); }

    /** Store an attribute on this LinearSystem which can be retrieved
        later by name. */
    virtual int putAttribute(const char* name,
                             void* attribute);

    /** Retrieve an attribute which was previously stored. */
    virtual int getAttribute(const char* name,
                             void*& attribute);

    /** Essential boundary-condition function that simply accepts a list
        of prescribed values, rather than the 'old' FEI's confusing approach
        of accepting arrays of alpha, beta and gamma values that nobody every
        really understood.

        For each specified ID, a value is being prescribed for a specified
        fieldID and a specified offset into that field.

        @param numIDs
        @param IDs
        @param idType
        @param fieldID
        @param offsetIntoField
        @param prescribedValues Input. List of values. Has length numIDs.
    */
    virtual int loadEssentialBCs(int numIDs,
                                 const int* IDs,
                                 int idType,
                                 int fieldID,
                                 int offsetIntoField,
                                 const double* prescribedValues);

    /** Essential boundary-condition function that simply accepts a list
        of prescribed values, rather than the 'old' FEI's confusing approach
        of accepting arrays of alpha, beta and gamma values that nobody every
        really understood.

        For each specified ID, a value is being prescribed for a specified
        fieldID and a specified offset into that field. The offset into the
        field can be different for each prescribed value.

        @param numIDs
        @param IDs
        @param idType
        @param fieldID
        @param offsetsIntoField Input. List of values, length numIDs.
        @param prescribedValues Input. List of values. Has length numIDs.
    */
    virtual int loadEssentialBCs(int numIDs,
                                 const int* IDs,
                                 int idType,
                                 int fieldID,
                                 const int* offsetsIntoField,
                                 const double* prescribedValues);

    /** Lagrange constraint coefficient loading function.
	@param constraintID Input. Must be an identifier of a lagrange 
	constraint that was initialized on the fei::MatrixGraph object which
	was used to construct the matrix for this linear system.
	@param weights Input. List, with length given by the sum of the sizes
	of the constrained fields.
	@param rhsValue
    */
    virtual int loadLagrangeConstraint(int constraintID,
				       const double *weights,
				       double rhsValue) = 0;

    /** Penalty constraint coefficient loading function.
	@param constraintID Input. Must be an identifier of a lagrange 
	constraint that was initialized on the fei::MatrixGraph object which
	was used to construct the matrix for this linear system.
	@param weights Input. List, with length given by the sum of the sizes
	of the constrained fields.
	@param penaltyValue
	@param rhsValue
    */
    virtual int loadPenaltyConstraint(int constraintID,
				      const double *weights,
				      double penaltyValue,
				      double rhsValue) = 0;

    /** Signal that all boundary-conditions and constraint coefficients have been
	loaded, and they may now be applied to the linear system.
    */
    virtual int loadComplete(bool applyBCs=true,
                             bool globalAssemble=true) = 0;

    /** Request that any boundary-condition values that have been provided via
	loadEssentialBCs() be set in the specified vector.
    */
    virtual int setBCValuesOnVector(fei::Vector* vector) = 0;

    /** Query whether a specified equation-index has a prescribed
	essential boundary-condition.
    */
    virtual bool eqnIsEssentialBC(int globalEqnIndex) const = 0;

    /** Fill caller-supplied vectors with the global equation-indices (which
	reside on the local processor) that have essential boundary-conditions
	prescribed, and fill a second vector with the prescribed values.
    */
    virtual void getEssentialBCs(std::vector<int>& bcEqns,
                                 std::vector<double>& bcVals) const = 0;

    /** Fill a caller-supplied vector with the global equation-indices (which
       reside on the local processor) that are involved in constraint-relations.
    */
    virtual void getConstrainedEqns(std::vector<int>& crEqns) const = 0;

   protected:
    fei::SharedPtr<fei::Matrix> matrix_;
    fei::SharedPtr<fei::Vector> soln_;
    fei::SharedPtr<fei::Vector> rhs_;

    fei::SharedPtr<fei::MatrixGraph> matrixGraph_;
    fei::DirichletBCManager* dbcManager_;

    std::vector<char*> attributeNames_;
    std::vector<void*> attributes_;
  };//class LinearSystem
}//namespace fei

#endif // _fei_LinearSystem_hpp_
