#ifndef PIKE_BLACK_BOX_MODEL_EVALUATOR_HPP
#define PIKE_BLACK_BOX_MODEL_EVALUATOR_HPP

#include "Pike_Response.hpp"
#include <string>

namespace pike {

  class BlackBoxModelEvaluator {

  public:

    virtual ~BlackBoxModelEvaluator() {}

    //! Unique name for this model evaluator.
    virtual std::string name() const = 0;

    /** \brief Perform a complete solve of this physics.

	@return Returns true if the solve is successful.
     */
    virtual bool solve() = 0;

    /** \brief Returns true if the last call to solve() was
     *  successful.  Note that this does NOT imply convergence of the
     *  global coupled system but only convergence of the local
     *  application solve.
    */
    virtual bool isConverged() const = 0;

    //! Returns the response for index i.
    virtual Teuchos::RCP<pike::Response> getResponse(const int i) const = 0;

    //! Returns the response index for the string name.
    virtual int getResponseIndex(const std::string name) const = 0;

    //! Returns true if the response is supported by this model evaluator
    virtual bool supportsResponse(const std::string name) const = 0;

  };

}

#endif
