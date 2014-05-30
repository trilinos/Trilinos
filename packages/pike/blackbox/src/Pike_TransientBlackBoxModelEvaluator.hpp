#ifndef PIKE_TRANSIENT_BLACK_BOX_MODEL_EVALUATOR_HPP
#define PIKE_TRANSIENT_BLACK_BOX_MODEL_EVALUATOR_HPP

#include "Pike_BlackBoxModelEvaluator.hpp"

namespace pike {

  /** A pure virtual base class that extends the pike::BlackBoxModelEvaluator for transient problems.

      This interface allows pike to query physics models for a
      possible next time step, then coordinate with the coupled codes
      to set a consistent step across the physics models and then
      solve from the pike::BlackBoxModelEvaluator.

      NOTE: During a time step, the solver will typically use picard
      iteration to converge the coupled system for the specified time
      step.  This means that each physics model will be solved with a
      call the the pike::BlackBoxModelEvaluator::solve() multiple
      times FOR THE SAME TIME STEP.  When the simulation is converged
      and ready to move on to a new time step, each physics model is
      notified by a call to acceptTimeStep().  If any of the codes
      failed for a specified time step, the solver can set a new
      reduced time step and attempt a new solve of the coupled system.
   */
  template<typename Scalar>
  class TransientBlackBoxModelEvaluator {

  public:

    virtual ~TransientBlackBoxModelEvaluator() {}

    /** \brief Returns the last time value that was "accepted" via a call to pike::TransientBlackBoxModelEvaluator::acceptTimeStep().  If no steps have been accepted, this value is the initial start time of the transient simulation time. */
    virtual Scalar getCurrentTime() const = 0;

    /** \brief Returns the current time that was used in the last call to pike::BlackBoxModelEvaluator::solve().  If this value is different than the accepted value, then a call to solve() has been performed, but the time step has not been accepted yet. **/
    virtual Scalar getTentativeTime() const = 0;

    /** \brief Returns true if a tentative time step has been solved but not accepted yet */
    virtual bool solvedTentativeStep() const = 0; 

    /** \brief Returns the currently set time step size. */
    virtual Scalar getCurrentTimeStepSize() const = 0;

    /** \brief Returns the time step size that the physics model wants to take next (e.g. to satisfy its stability or local truncation error targets). */ 
    virtual Scalar getDesiredTimeStepSize() const = 0;

    /** \brief Returns the maximum time step size that a physics model can handle (e.g. to satisfy a CFL limit for an explicit time integration method). */
    virtual Scalar getMaxTimeStepSize() const = 0;

    /** \brief Sets the time step size for the next solve. */
    virtual void setNextTimeStepSize(const Scalar& dt) = 0;

    /** \brief Accepts a tentative time step solve. */
    virtual void acceptTimeStep() = 0;

  };

}

#endif
