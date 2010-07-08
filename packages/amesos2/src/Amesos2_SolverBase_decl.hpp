// @HEADER
// ***********************************************************************
//
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

/**
  \file   Amesos2_SolverBase_decl.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Thu Dec 31 10:32:08 2009
  
  \brief  Pure virtual base class used by all other concrete solver
          classes.
*/

#ifndef AMESOS2_SOLVER_BASE_DECL_HPP
#define AMESOS2_SOLVER_BASE_DECL_HPP

#include <Teuchos_ParameterListAcceptor.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

#include "Amesos2_Control.hpp"
#include "Amesos2_Status.hpp"
#include "Amesos2_Timers.hpp"

namespace Amesos {

  
class SolverBase : Teuchos::ParameterListAcceptor {
  
public:
  
  //@{ \name Solver-dependent methods
  virtual SolverBase& preOrdering( void ) = 0;
    

  virtual SolverBase& symbolicFactorization( void ) = 0;


  virtual SolverBase& numericFactorization( void ) = 0;


  virtual void solve( void ) = 0;


  virtual bool matrixShapeOK( void ) = 0;


  virtual SolverBase& setParameters(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList ) = 0;

    
  //@} End solver-dependent methods

  virtual const Teuchos::RCP<const Teuchos::Comm<int> >& getComm( void ) const = 0;


  virtual int getNumSymbolicFact( void ) const = 0;


  virtual int getNumNumericFact( void ) const = 0;


  virtual int getNumSolve( void ) const = 0;


  virtual void describe( Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;


  virtual void printTiming( Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default ) const = 0;

  /**
     Redefined from Teuchos::ParameterListAcceptor

     \param parameterList
  */
  virtual void setParameterList( Teuchos::RCP<Teuchos::ParameterList> const& parameterlist ) = 0;


  /// This is an empty stub
  /** 
      This is an empty stub
      
      \return 
  */
  virtual Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList( void ) = 0;


  /// This is an empty stub
  /** 
      This is an empty stub
      
      \return 
  */
  virtual Teuchos::RCP<Teuchos::ParameterList> unsetParameterList( void ) = 0;


  /** \brief Extracts timing information from the current solver.
   *
   * Results are placed in the parameter list \c timingParameterList
   * 
   * \param timingParameterList Accepts timing information from the
   * current solver
   */
  virtual void getTiming( Teuchos::ParameterList& timingParameterList ) const = 0;


  // Class members
protected:

  // Hold status information about a solver
  Status status_;

  // Parameters for solving
  Control control_;

  // Various timing statistics
  Timers timers_;
    
};                              // End class SolverBase

  
} // end namespace Amesos

#endif	// AMESOS2_SOLVER_BASE_DECL_HPP
