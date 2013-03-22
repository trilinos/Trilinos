// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#ifndef LOCA_FACTORY_H
#define LOCA_FACTORY_H

#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H"

#include "LOCA_MultiPredictor_Factory.H"
#include "LOCA_MultiContinuation_Factory.H"
#include "LOCA_Bifurcation_Factory.H"
#include "LOCA_StepSize_Factory.H"
#include "LOCA_BorderedSolver_Factory.H"
#include "LOCA_Eigensolver_Factory.H"
#include "LOCA_EigenvalueSort_Factory.H"
#include "LOCA_SaveEigenData_Factory.H"
#include "LOCA_AnasaziOperator_Factory.H"
#include "LOCA_TurningPoint_MooreSpence_SolverFactory.H"
#include "LOCA_Pitchfork_MooreSpence_SolverFactory.H"
#include "LOCA_Hopf_MooreSpence_SolverFactory.H"

// Forward declarations
namespace Teuchos {
  class ParameterList;
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
    Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);
    
    //! Constructor with user-supplied factory
    /*!
     * \param global_data [in] Global data object.  The constructor sets the
     * factory member of the global data to this.
     * \param userFactory [in] A user-supplied factory for instantiating
     * user-defined strategies.
     */
    Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data,
	    const Teuchos::RCP<LOCA::Abstract::Factory>& userFactory);

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
    Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>
    createPredictorStrategy(
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& predictorParams);

    //! Create continuation strategy
    /*!
     * Instantiates a continuation strategy based on the "Continuation Method" 
     * parameter of the "Stepper" sublist.  See 
     * LOCA::MultiContinuation::Factory for a description of available 
     * strategies.
     */
    Teuchos::RCP<LOCA::MultiContinuation::AbstractStrategy>
    createContinuationStrategy(
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& stepperParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const std::vector<int>& paramIDs);

    //! Create bifurcation strategy
    /*!
     * Instantiates a bifurcation strategy based on the "Method" 
     * parameter of the "Bifurcation" sublist.  See 
     * LOCA::Bifurcation::Factory for a description of available 
     * strategies.
     */
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
    createBifurcationStrategy(
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& bifurcationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp);

    //! Create step size control strategy
    /*!
     * Instantiates a step size control strategy based on the 
     * "Method" parameter of the "Step Size" sublist.  See 
     * LOCA::StepSize::Factory for a description of available strategies.
     */
    Teuchos::RCP<LOCA::StepSize::AbstractStrategy>
    createStepSizeStrategy(
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& stepsizeParams);

    //! Create bordered system solver strategy
    /*!
     * Instantiates an bordered system solver strategy based on the 
     * "Bordered Solver Method" parameter
     * of the "Linear Solver" sublist.  See LOCA::BorderedSolver::Factory for 
     * a description of available strategies.
     */
    Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy>
    createBorderedSolverStrategy(
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

    //! Create eigensolver strategy
    /*!
     * Instantiates an eigensolver strategy based on the "Method" parameter
     * of the "Eigensolver" sublist.  See LOCA::Eigensolver::Factory for a 
     * description of available strategies.
     */
    Teuchos::RCP<LOCA::Eigensolver::AbstractStrategy>
    createEigensolverStrategy(
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);

    //! Create eigenvalue sort strategy
    /*!
     * Instantiates an eigenvalue sorting strategy based on the 
     * "Sorting Method" parameter of the "Eigensolver" sublist.  See 
     * LOCA::EigenvalueSort::Factory for a description of available strategies.
     */
    Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy>
    createEigenvalueSortStrategy(
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);

    //! Create strategy to save eigenvector/value data
    /*!
     * Instantiates a strategy to save eigenvector/value data based on the 
     * "Save Eigen Data Method" parameter of the "Eigensolver" sublist.  
     * See LOCA::SaveEigenData::Factory for a description of available 
     * strategies.
     */
    Teuchos::RCP<LOCA::SaveEigenData::AbstractStrategy>
    createSaveEigenDataStrategy(
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);

    //! Create Anasazi operator strategy
    /*!
     * Instantiates an Anasazi operator strategy based on the "Operator" 
     * parameter of the "Eigensolver" sublist.  See 
     * LOCA::AnasaziOperator::Factory for a description of available 
     * strategies.
     */
    Teuchos::RCP<LOCA::AnasaziOperator::AbstractStrategy>
    createAnasaziOperatorStrategy(
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& eigenParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& solverParams,
	 const Teuchos::RCP<NOX::Abstract::Group>& grp);

    //! Create Moore-Spence turning point solver strategy
    /*!
     * Instantiates a solver strategy based on the 
     * "Solver Method" parameter of the "Bifurcation" sublist.  See
     * LOCA::TurningPoint::MooreSpence::SolverFactory for a description of 
     * available strategies.
     */
    Teuchos::RCP<LOCA::TurningPoint::MooreSpence::SolverStrategy>
    createMooreSpenceTurningPointSolverStrategy(
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

    //! Create Moore-Spence pitchfork solver strategy
    /*!
     * Instantiates a solver strategy based on the 
     * "Solver Method" parameter of the "Bifurcation" sublist.  See
     * LOCA::Pitchfork::MooreSpence::SolverFactory for a description of 
     * available strategies.
     */
    Teuchos::RCP<LOCA::Pitchfork::MooreSpence::SolverStrategy>
    createMooreSpencePitchforkSolverStrategy(
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

    //! Create Moore-Spence Hopf solver strategy
    /*!
     * Instantiates a solver strategy based on the 
     * "Solver Method" parameter of the "Bifurcation" sublist.  See
     * LOCA::Hopf::MooreSpence::SolverFactory for a description of 
     * available strategies.
     */
    Teuchos::RCP<LOCA::Hopf::MooreSpence::SolverStrategy>
    createMooreSpenceHopfSolverStrategy(
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

    //@}

  private:

    //! Private to prohibit copying
    Factory(const Factory&);

    //! Private to prohibit copying
    Factory& operator = (const Factory&);

  protected:

    //! Global data
    Teuchos::RCP<LOCA::GlobalData> globalData;

    //! User provided factory
    Teuchos::RCP<LOCA::Abstract::Factory> factory;

    //! Flag indicating whether we have a factory or not
    bool haveFactory;

    //! Predictor factory
    LOCA::MultiPredictor::Factory predictorFactory;

    //! Continuation factory
    LOCA::MultiContinuation::Factory continuationFactory;

    //! Bifurcation factory
    LOCA::Bifurcation::Factory bifurcationFactory;

    //! Step size factory
    LOCA::StepSize::Factory stepsizeFactory;

    //! Bordered system factory
    LOCA::BorderedSolver::Factory borderedFactory;

    //! Eigensolver factory
    LOCA::Eigensolver::Factory eigensolverFactory;

    //! EigenvalueSort factory
    LOCA::EigenvalueSort::Factory eigenvalueSortFactory;

    //! Save eigen data factory
    LOCA::SaveEigenData::Factory saveEigenFactory;

    //! Anasazi operator factory
    LOCA::AnasaziOperator::Factory anasaziOperatorFactory;

    //! Moore-Spence turning point solver factory
    LOCA::TurningPoint::MooreSpence::SolverFactory mooreSpenceTurningPointSolverFactory;

    //! Moore-Spence pitchfork point solver factory
    LOCA::Pitchfork::MooreSpence::SolverFactory mooreSpencePitchforkSolverFactory;

    //! Moore-Spence Hopf point solver factory
    LOCA::Hopf::MooreSpence::SolverFactory mooreSpenceHopfSolverFactory;
    
  }; // Class Factory

} // Namespace LOCA

#endif
