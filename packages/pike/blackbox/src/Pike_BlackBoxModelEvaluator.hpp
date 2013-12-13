#ifndef PIKE_BLACK_BOX_MODEL_EVALUATOR_HPP
#define PIKE_BLACK_BOX_MODEL_EVALUATOR_HPP

#include "Pike_Response.hpp"
#include <string>

namespace pike {

  class BlackBoxModelEvaluator {

    virtual ~BlackBoxModelEvaluator() {}

    virtual const std::string name() const = 0;

    virtual void solve() = 0;

    virtual bool isConverged() const = 0;

    virtual Teuchos::RCP<pike::Response> getResponse(const int i) const = 0;

    virtual int getResponseIndex(const std::string name) const = 0;

    virtual bool supportsResponse(const std::string name) const = 0;

  };

}

#endif
