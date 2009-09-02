/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Solver_hpp_
#define _fei_Solver_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>

namespace fei {
  class ParameterSet;
  class Matrix;
  class LinearSystem;

  /** Interface for requesting that a linear-system be solved.
   */
  class Solver {
  public:
    /** Solver Factory interface */
    class Factory {
    public:
      /** Usual virtual destructor */
      virtual ~Factory(){}

      /** Produce an instance of a Solver */
      virtual fei::SharedPtr<fei::Solver> createSolver(const char* name=0) = 0;
    };

    /** virtual destructor */
    virtual ~Solver(){}

    /** Solve a linear system
     */
    virtual int solve(fei::LinearSystem* linearSystem,
		      fei::Matrix* preconditioningMatrix,
		      const fei::ParameterSet& parameterSet,
		      int& iterationsTaken,
		      int& status);
  };//class Solver
}//namespace fei

#endif // _fei_Factory_hpp_
