// @HEADER
// ***********************************************************************
//
//                Amesos2: Direct Sparse Solver Package
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
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

#include "Amesos2_Control.hpp"
#include "Amesos2_Status.hpp"
#include "Amesos2_Timers.hpp"

namespace Amesos {


/**
 * \brief Base class used by all other concrete solver classes.
 *
 * Specifies a uniform interface for interaction with Amesos2 solver wrappers
 * to third-party libraries.
 *
 * This class holds no definitions itself; it is pure virtual.
 */
class SolverBase : Teuchos::ParameterListAcceptor {

public:

  /// \name Solver-dependent methods
  //@{

  /// Pre-orders the matrix for minimal fill-in
  virtual SolverBase& preOrdering( void ) = 0;


  /// Performs symbolic factorization on the matrix
  virtual SolverBase& symbolicFactorization( void ) = 0;


  ///Performs numeric factorization on the matrix
  virtual SolverBase& numericFactorization( void ) = 0;


  /// Solves \f$ A X = B\f$ (or \f$ A^T X = B\f$ )
  virtual void solve( void ) = 0;


  /// Returns \c true if the solver can handle the matrix shape
  virtual bool matrixShapeOK( void ) = 0;


  /// Set/update internal variables and solver options.
  virtual SolverBase& setParameters(
    const Teuchos::RCP<Teuchos::ParameterList> & parameterList ) = 0;


  /**
   * \brief Return a const parameter list of all of the valid parameters that
   * this->setParameterList(...)  will accept.
   */
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const = 0;

  //@} End solver-dependent methods


  /// Returns a pointer to the Tpetra::Comm communicator with this matrix
  virtual Teuchos::RCP<const Teuchos::Comm<int> > getComm( void ) const = 0;


  /// Returns the number of symbolic factorizations performed by this object.
  virtual int getNumSymbolicFact( void ) const = 0;


  /// Returns the number of numeric factorizations performed by this object.
  virtual int getNumNumericFact( void ) const = 0;


  /// Returns the number of solves performed by this object.
  virtual int getNumSolve( void ) const = 0;


  /// Returns a short description of this Solver
  virtual std::string description() const = 0;


  /// Prints the status information about the current solver with some level
  /// of verbosity.
  virtual void describe( Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;


  /// Prints timing information about the current solver.
  virtual void printTiming( Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default ) const = 0;


  /**
   * Redefined from Teuchos::ParameterListAcceptor
   *
   * \param parameterList
   */
  virtual void setParameterList( Teuchos::RCP<Teuchos::ParameterList> const& parameterlist ) = 0;


  /**
   * \brief This is an empty stub
   */
  virtual Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList( void ) = 0;


  /**
   * \brief This is an empty stub
   */
  virtual Teuchos::RCP<Teuchos::ParameterList> unsetParameterList( void ) = 0;


  /**
   * \brief Extracts timing information from the current solver.
   *
   * Results are placed in the parameter list \c timingParameterList
   *
   * \param timingParameterList Accepts timing information from the
   * current solver
   */
  virtual void getTiming( Teuchos::ParameterList& timingParameterList ) const = 0;


  /// Return the name of this solver.
  virtual std::string name() const = 0;


};                              // End class SolverBase


} // end namespace Amesos

#endif	// AMESOS2_SOLVER_BASE_DECL_HPP
