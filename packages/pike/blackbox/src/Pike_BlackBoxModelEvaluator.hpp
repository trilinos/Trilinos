#ifndef PIKE_BLACK_BOX_MODEL_EVALUATOR_HPP
#define PIKE_BLACK_BOX_MODEL_EVALUATOR_HPP

#include "Teuchos_ArrayView.hpp"
#include <string>

namespace pike {

  /** \brief Pure virtual interface to a user implemented physics model. */
  class BlackBoxModelEvaluator {

  public:

    virtual ~BlackBoxModelEvaluator() {}

    //! Unique name for this model evaluator.
    virtual std::string name() const = 0;

    /** \brief Perform a complete solve of this physics.

	@return Returns true if the local solve is successful.
     */
    virtual bool solve() = 0;

    /** \brief Returns true if the last call to solve() was
     *  successful.  Note that this does NOT imply convergence of the
     *  global coupled system but only convergence of the individual
     *  (local) application solve.
     */
    virtual bool isLocallyConverged() const = 0;

    /** \brief Optional function for assessing convergence of the
     *	globally coupled problem.  Returns true if the metrics local
     *	to this application are converged.  This allows for users to
     *	track responses internally and determine convergence of the
     *	global coupled system based on relative error of local
     *	responses without having to expose the responses through a
     *	pike::Response.  This is primarily for backwards
     *	compatibility, new codes should use the StatusTest objects for
     *	global convergence (and for checking the local convergence).
     */
    virtual bool isGloballyConverged() const = 0;

    //! Returns the response for index i.
    virtual Teuchos::ArrayView<const double> getResponse(const int i) const = 0;

    //! Returns the response index for the string name.
    virtual int getResponseIndex(const std::string& name) const = 0;

    //! Returns the response name for the index.
    virtual std::string getResponseName(const int i) const = 0;

    //! Returns true if the response is supported by this model evaluator.
    virtual bool supportsResponse(const std::string& name) const = 0;

    //! Get the number of responses, Ng.
    virtual int getNumberOfResponses() const = 0;

  };

}

#endif
