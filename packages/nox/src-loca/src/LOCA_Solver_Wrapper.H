// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_SOLVER_WRAPPER_H
#define LOCA_SOLVER_WRAPPER_H

#include "NOX_Solver_Generic.H"	// base class

namespace LOCA {

  //! %LOCA %Solver classes
  namespace Solver{

    /*!
     * \brief A %wrapper class for wrapping a %NOX solver.
     */
    /*!
     * The %LOCA::Solver::Wrapper class provides a wrapper for %NOX solvers
     * to change the group data returned by getSolutionGroup() and 
     * getPreviousSolutionGroup() so that status tests can operate correctly.
     * %LOCA continuation and bifurcation algorithms are implemented by
     * extended groups which augment the nonlinear equations defining
     * equilibrium solutions with appropriate continuation/bifurcation
     * equations.  Therefore the groups returned by getSolutionGroup() and
     * getPreviousSolutionGroup() will be groups corresponding to the 
     * continuation/bifurcation equations, not the underlying group.
     * Status tests that are designed to use concrete data from the original
     * underlying group (e.g., NOX::StatusTest::NormWRMS) will then fail
     * (usually via a segmentation fault or raised exception) when applied
     * to the extended continuation or bifurcation groups.  
     *
     * This solver wrapper class fixes this problem by reimplementing the
     * solution group accessor methods to return the underlying groups of 
     * the solution groups 
     * (via LOCA::Extended::MultiAbstractGroup::getUnderlyingGroup()) if they 
     * are extended groups (derived from LOCA::Extended::MultiAbstractGroup).  
     * If the groups are not extended groups, the original solution groups are 
     * returned.  All other NOX::Solver::Generic methods are passed to the
     * wrapped solver.
     *
     * The LOCA::StatusTest::Wrapper class uses this wrapper class to wrap
     * the solver supplied via the checkStatus method which is then forwarded
     * to the original status test.  When used properly, the group "seen" by
     * the status test is of the appropriate type for the status test.
     */
    class Wrapper : public NOX::Solver::Generic {

    public:

      //! Constructor with a non-const ref-count pointer to a %NOX solver
      /*!
       * The constructor calls resetWrapper() to grab the proper solution
       * groups.
       */
      Wrapper(const Teuchos::RefCountPtr<NOX::Solver::Generic>& solver);

      //! Constructor with a const ref-count pointer to a %NOX solver
      /*!
       * The constructor calls resetWrapper() to grab the proper solution
       * groups.
       */
      Wrapper(const Teuchos::RefCountPtr<const NOX::Solver::Generic>& solver);

      //! Destructor
      ~Wrapper();

      //! Implementation of reset method (forwarded to wrapped solver)
      bool reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
		 const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests, 
		 const Teuchos::RefCountPtr<Teuchos::ParameterList>& params);

      //! Implementation of reset method (forwarded to wrapped solver)
      bool reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp);
  
      //! Implementation of reset method (forwarded to wrapped solver)
      bool reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
		 const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests);
      
      //! Implementation of getStatus method (forwarded to wrapped solver)
      NOX::StatusTest::StatusType getStatus();
      
      //! Implementation of step method (forwarded to wrapped solver)
      NOX::StatusTest::StatusType step();
      
      //! Implementation of solve method (forwarded to wrapped solver)
      NOX::StatusTest::StatusType solve();
      
      /*!
       * \brief Returns underlying group if solution group is extended, 
       * solution group otherwise.
       */
      const NOX::Abstract::Group& getSolutionGroup() const;
      
      /*!
       * \brief Returns underlying group if previous solution group is 
       * extended, previous solution group otherwise.
       */
      const NOX::Abstract::Group& getPreviousSolutionGroup() const;
      
      /*! 
       * \brief Implementation of getNumIterations method (forwarded to 
       * wrapped solver)
       */
      int getNumIterations() const;
      
      /*! 
       * \brief Implementation of getList method (forwarded to 
       * wrapped solver)
       */
      const Teuchos::ParameterList& getList() const;

    protected:

      //! Grabs appropriate solution groups from solver
      void resetWrapper();

    protected:
      
      //! Pointer to the solver object
      Teuchos::RefCountPtr<NOX::Solver::Generic> solverPtr;

      //! Pointer to the const solver object
      Teuchos::RefCountPtr<const NOX::Solver::Generic> constSolverPtr;

      //! Pointer to solution group
      Teuchos::RefCountPtr<const NOX::Abstract::Group> solnGrpPtr;

      //! Pointer to previous solution group
      Teuchos::RefCountPtr<const NOX::Abstract::Group> oldSolnGrpPtr;

    };
  }
}
#endif
