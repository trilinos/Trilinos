#ifndef PIKE_BLACK_BOX_MODEL_EVALUATOR_HPP
#define PIKE_BLACK_BOX_MODEL_EVALUATOR_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Pike_BlackBoxModelEvaluator_ParameterMixIn.hpp"
#include "Pike_BlackBoxModelEvaluator_ResponseMixIn.hpp"
#include <string>

namespace pike {

  /** \brief Pure virtual interface to a user implemented physics model. */
  class BlackBoxModelEvaluator : public Teuchos::Describable,
				 public Teuchos::VerboseObject<pike::BlackBoxModelEvaluator>,
				 public pike::ParameterMixIn,
				 public pike::ResponseMixIn {

  public:

    virtual ~BlackBoxModelEvaluator() {}

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
     *	consider the global coupled problem to be converged.

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
     *	because it assume that each call to solve is for a new
     *	iteration (i.e. assumes implicit knowledge of the solver
     *	algorithm).  If a solve is called for some other purpose
     *	(perturbing for sensitivities), then relative change of
     *	internal responses is not a good indicator!  We strongly
     *	recommend that you use pike::StatusTest objects for
     *	determining global convergence!
     */
    virtual bool isGloballyConverged() const = 0;

  };

}

#endif
