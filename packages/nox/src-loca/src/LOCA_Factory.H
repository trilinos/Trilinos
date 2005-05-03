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

#include "NOX_Abstract_Group.H"

#include "LOCA_MultiPredictor_Factory.H"
#include "LOCA_MultiContinuation_Factory.H"
#include "LOCA_BorderedSystem_Factory.H"
#include "LOCA_Eigensolver_Factory.H"
#include "LOCA_EigenvalueSort_Factory.H"

// Forward declarations
namespace NOX {
  namespace Parameter {
    class List;
  }
}
namespace LOCA {
  class GlobalData;
  namespace Abstract {
    class Factory;
  }
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  //! %Factory class for creating strategies
  /*!
   * The %Factory class provides a single location for instantiating
   * various strategies based on parameter list choices.  It provides a 
   * create() method for each type of strategy which instantiates strategy
   * objects for that type.  Each create method takes as arguments a 
   * ref-count pointer to a LOCA::Parameter::SublistParser and a parameter
   * list.  The parameter list determines which strategy to choose and also
   * should provide any parameters the strategy requires.  The sublist parser
   * provides a parsed version of the top-level parameter list and allows
   * strategies to easily obtain other sublists from the top-level list.
   * A user-supplied factory may also be provided for instantiating 
   * user-defined strategies.  If a user-defined factory is supplied, 
   * each create method will first attempt to instantiate the strategy using 
   * it, and then instantiate strategies itself if necessary.
   */
  class Factory {
    
  public:

    //! Constructor
    /*!
     * \param global_data [in] Global data object.  The constructor sets the
     * factory member of the global data to this.
     */
    Factory(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data);
    
    //! Constructor with user-supplied factory
    /*!
     * \param global_data [in] Global data object.  The constructor sets the
     * factory member of the global data to this.
     * \param userFactory [in] A user-supplied factory for instantiating
     * user-defined strategies.
     */
    Factory(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	    const Teuchos::RefCountPtr<LOCA::Abstract::Factory>& userFactory);

    //! Destructor
    virtual ~Factory();

    /*! 
     * @name Strategy create methods
     */
    //@{

    //! Create predictor strategy
    /*!
     * Instantiates a predictor strategy based on the 
     * "Method" parameter of the "Predictor" sublist.  See 
     * LOCA::MultiPredictor::Factory for a description of available strategies.
     */
    Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>
    createPredictorStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& predictorParams);

    //! Create continuation strategy
    /*!
     * Instantiates a continuation strategy based on the "Continuation Method" 
     * parameter of the "Stepper" sublist.  See 
     * LOCA::MultiContinuation::Factory for a description of available 
     * strategies.
     */
    Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractStrategy>
    createContinuationStrategy(
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<NOX::Parameter::List>& stepperParams,
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const vector<int>& paramIDs);

    //! Create bordered system solver strategy
    /*!
     * Instantiates an bordered system solver strategy based on the 
     * "Bordered Solver Method" parameter
     * of the "Linear Solver" sublist.  See LOCA::BorderedSystem::Factory for 
     * a description of available strategies.
     */
    Teuchos::RefCountPtr<LOCA::BorderedSystem::AbstractStrategy>
    createBorderedSystemStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& solverParams);

    //! Create eigensolver strategy
    /*!
     * Instantiates an eigensolver strategy based on the "Method" parameter
     * of the "Eigensolver" sublist.  See LOCA::Eigensolver::Factory for a 
     * description of available strategies.
     */
    Teuchos::RefCountPtr<LOCA::Eigensolver::AbstractStrategy>
    createEigensolverStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& eigenParams);

    //! Create eigenvalue sort strategy
    /*!
     * Instantiates an eigenvalue sorting strategy based on the 
     * "Sorting Method" parameter of the "Eigensolver" sublist.  See 
     * LOCA::EigenvalueSort::Factory for a description of available strategies.
     */
    Teuchos::RefCountPtr<LOCA::EigenvalueSort::AbstractStrategy>
    createEigenvalueSortStrategy(
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<NOX::Parameter::List>& eigenParams);

    //@}

  private:

    //! Private to prohibit copying
    Factory(const Factory&);

    //! Private to prohibit copying
    Factory& operator = (const Factory&);

  protected:

    //! Global data
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

    //! User provided factory
    Teuchos::RefCountPtr<LOCA::Abstract::Factory> factory;

    //! Flag indicating whether we have a factory or not
    bool haveFactory;

    //! Predictor factory
    LOCA::MultiPredictor::Factory predictorFactory;

    //! Continuation factory
    LOCA::MultiContinuation::Factory continuationFactory;

    //! Bordered system factory
    LOCA::BorderedSystem::Factory borderedFactory;

    //! Eigensolver factory
    LOCA::Eigensolver::Factory eigensolverFactory;

    //! EigenvalueSort factory
    LOCA::EigenvalueSort::Factory eigenvalueSortFactory;
    
  }; // Class Factory

} // Namespace LOCA

#endif
