#ifndef PIKE_BLACK_BOX_MODEL_EVALUATOR_HPP
#define PIKE_BLACK_BOX_MODEL_EVALUATOR_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include <string>

namespace pike {

  /** \brief Pure virtual interface to a user implemented physics model. */
  class BlackBoxModelEvaluator : public Teuchos::Describable,
				 public Teuchos::VerboseObject<pike::BlackBoxModelEvaluator> {

  public:

    virtual ~BlackBoxModelEvaluator();

    //! Unique name for this model evaluator.
    virtual std::string name() const = 0;

    /** \brief Perform a complete solve of this physics. */
    virtual void solve() = 0;

    /** \brief Returns true if the last call to solve() was
     *  successful.  Note that this does NOT imply convergence of the
     *  global coupled system but only convergence of the individual
     *  (local) application solve.
     */
    virtual bool isLocallyConverged() const = 0;

    /** \brief Returns true if the metrics local to this application
     *	consider the global coupled problem to be converged.  NOTE:
     *	This function is deprecated.

     *  This method allows individual applications to track responses
     *	internally and determine convergence of the global coupled
     *	system.  For example, an application might monitor the
     *	relative error of local responses between iterations without
     *	having to expose the responses through the model evaluator
     *	interface.  This is primarily for user convenience, codes
     *	should really be using the StatusTest objects for checking
     *	global convergence.

     *  IMPORTANT NOTE: It is dangerous to use this function for
     *	global convergence as opposed to pike::StatusTest objects
     *	because it assumes that each call to solve is for a new
     *	iteration (i.e. assumes implicit knowledge of the solver
     *	algorithm).  If a solve is called for some other purpose
     *	(perturbing for sensitivities), then relative change of
     *	internal responses is not a good indicator!  We strongly
     *	recommend that you use pike::StatusTest objects for
     *	determining global convergence!
     */
    virtual bool isGloballyConverged() const;

    /**@{ \name Optional Support for Parameters.
  
       This group of methods is optional: default methods are
       implemented.  Users can override the default implementations if
       the particular application model evaluator supports setting
       parameters.
     */
    
    //! Returns true if the parameter is supported by this model evaluator.
    virtual bool supportsParameter(const std::string& pName) const;
    
    //! Get the number of parameters, Np.
    virtual int getNumberOfParameters() const;
    
    //! Returns the parameter name for index l where 0 <= l < Np.
    virtual std::string getParameterName(const int l) const;
    
    //! Returns the parameter index for the parameter name.
    virtual int getParameterIndex(const std::string& pName) const;
    
    //! Sets the parameter, p, for index l where 0 <= l < Np. 
    virtual void setParameter(const int l, const Teuchos::ArrayView<const double>& p);
    
    /**@} */
    
    /**@{ \name Optional Support for Responses.
       
       This group of methods is optional: default methods are
       implemented.  Users can override the default implementations if
       the particular application model evaluator supports returning
       responses.
    */
    
    //! Returns true if the response is supported by this model evaluator.
    virtual bool supportsResponse(const std::string& rName) const;
    
    //! Get the number of responses, Ng.
    virtual int getNumberOfResponses() const;
    
    //! Returns the response name for index j where 0 <= j < Ng.
    virtual std::string getResponseName(const int j) const;
    
    //! Returns the response index for the string name.
    virtual int getResponseIndex(const std::string& rName) const;
    
    //! Returns the response for index j where 0 <= j < Ng. 
    virtual Teuchos::ArrayView<const double> getResponse(const int j) const;
    
    /**@} */
    
    /**@{ \name Optional Support for Transient Applications.
       
       This group of methods is optional: default methods are
       implemented.  Users can override the default implementations if
       the particular application model evaluator supports a transient
       solve.
       
       This interface allows pike to query physics models for a
       possible next time step, then coordinate with the coupled codes
       to set a consistent step across the physics models and then
       solve() from the pike::BlackBoxModelEvaluator base class
       method.
       
       NOTE: During a time step, the solver will typically use picard
       iteration to converge the coupled system for the specified time
       step.  This means that each physics model will be solved with a
       call the the pike::BlackBoxModelEvaluator::solve() multiple
       times FOR THE SAME TIME STEP.  When the simulation is converged
       and ready to move on to a new time step, each physics model is
       notified by a call to acceptTimeStep().  If any of the codes
       failed for a specified time step, the solver can set a new
       reduced time step and attempt a new solve of the coupled
       system.
    */
    
    /** \brief Returns true if transient support is enabled for this
	model evaluator.
    */
    virtual bool isTransient() const;
    
    /** \brief Returns the last time value that was "accepted" via a
	call to
	pike::TransientBlackBoxModelEvaluator::acceptTimeStep().  If
	no steps have been accepted, this value is the initial start
	time of the transient simulation time.
    */
    virtual double getCurrentTime() const;
    
    /** \brief Returns the current time that was used in the last call
	to pike::BlackBoxModelEvaluator::solve().  If this value is
	different than the accepted value, then a call to solve() has
	been performed, but the time step has not been accepted yet.
    */
    virtual double getTentativeTime() const;
    
    /** \brief Returns true if a tentative time step has been solved
	but not accepted yet.
    */
    virtual bool solvedTentativeStep() const;
    
    /** \brief Returns the currently set time step size. */
    virtual double getCurrentTimeStepSize() const;
    
    /** \brief Returns the time step size that the physics model wants
	to take next (e.g. to satisfy its stability or local
	truncation error targets).
    */ 
    virtual double getDesiredTimeStepSize() const;
    
    /** \brief Returns the maximum time step size that a physics model
	can handle (e.g. to satisfy a CFL limit for an explicit time
	integration method).
    */
    virtual double getMaxTimeStepSize() const;

    /** \brief Sets the time step size for the next solve. */
    virtual void setNextTimeStepSize(const double& dt);

    /** \brief Accepts a tentative time step solve. */
    virtual void acceptTimeStep();

    /**@} */

  };

}

#endif
