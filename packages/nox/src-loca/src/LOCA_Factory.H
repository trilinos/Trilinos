// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#ifndef LOCA_FACTORY_H
#define LOCA_FACTORY_H

#include "Teuchos_RefCountPtr.hpp"

#include "LOCA_Abstract_Factory.H"
#include "LOCA_Parameter_SublistParser.H"

// Factories
#include "LOCA_Eigensolver_Factory.H"

namespace LOCA {

  //! %Factory class for creating strategies
  /*!
   * The %Factory class provides a single location for instantiating
   * various strategies based on parameter list choices.  It provides a 
   * create() method for each type of strategy which instantiates strategy
   * objects for that type.  The argument \c topLevelParams to the constructors
   * and reset() method should provide the name of the requested strategy
   * as well as any parameters for a particular strategy.  A user-supplied
   * factory may also be provided for instantiating user-defined strategies.
   * If a user-defined factory is supplied, each create method will first
   * attempt to instantiate the strategy using it, and then instantiate
   * strategies itself if necessary.
   */
  class Factory {
    
  public:

    //! Constructor
    /*!
     * \param topLevelparams [in] Top-level parameter list.  This list will
     * be parsed to determine which strategies to instantiate.
     */
    Factory(const Teuchos::RefCountPtr<NOX::Parameter::List>& topLevelParams);
    
    //! Constructor with user-supplied factory
    /*!
     * \param topLevelparams [in] Top-level parameter list.  This list will
     * be parsed to determine which strategies to instantiate.
     * \param userFactory [in] A user-supplied factory for instantiating
     * user-defined strategies.
     */
    Factory(const Teuchos::RefCountPtr<NOX::Parameter::List>& topLevelParams,
	    const Teuchos::RefCountPtr<LOCA::Abstract::Factory>& userFactory);

    //! Destructor
    ~Factory();

    //! Reset
    /*!
     * \param topLevelParams [in]  Re-parses parameter list 
     * This method should be called whenever a different strategy is required.
     */
    NOX::Abstract::Group::ReturnType
    reset(const Teuchos::RefCountPtr<NOX::Parameter::List>& topLevelParams);

    /*! 
     * @name Strategy create methods
     */
    //@{

    //! Create eigensolver strategy
    /*!
     * Instantiates an eigensolver strategy based on the "Method" parameter
     * of the "Eigensolver" sublist.  See LOCA::Eigensolver::Factory for a 
     * description of available strategies.
     */
    Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy>
    createEigensolver();

    //@}

  protected:

    //! Parse parameter list
    void parseParameters();

  protected:

    //! User provided factory
    Teuchos::RefCountPtr<LOCA::Abstract::Factory> factory;

    //! Flag indicating whether we have a factory or not
    bool haveFactory;

    //! Parser to find sublists
    LOCA::Parameter::SublistParser sublistParser;

    //! Eigensolver factory
    LOCA::Eigensolver::Factory eigensolverFactory;
    
  }; // Class Factory

} // Namespace LOCA

#endif
