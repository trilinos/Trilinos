/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_Solver_hpp_
#define _snl_fei_Solver_hpp_

#include <fei_macros.hpp>
#include <fei_Solver.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Vector.hpp>
#include <fei_Matrix.hpp>
#include <fei_LinearSystem.hpp>
#include <fei_LinearSystemCore.hpp>
#include <fei_FiniteElementData.hpp>

#include <fei_ErrMacros.hpp>

namespace snl_fei {

  /** snl_fei:: implementation of the fei::Solver interface.
   */
  class Solver : public virtual fei::Solver {
  public:
    /** Constructor */
    Solver(){}

    /** Destructor */
    virtual ~Solver(){}

    int solve(fei::LinearSystem* linearSystem,
	      fei::Matrix* preconditioningMatrix,
	      const fei::ParameterSet& parameterSet,
	      int& iterationsTaken,
	      int& status);

  private:
    int solve(fei::LinearSystem* linearSystem,
	      fei::Matrix* preconditioningMatrix,
	      int numParams,
	      const char* const* solverParams,
	      int& iterationsTaken,
	      int& status);

  };//class Solver
}//namespace snl_fei

#endif // _snl_fei_Solver_hpp_
