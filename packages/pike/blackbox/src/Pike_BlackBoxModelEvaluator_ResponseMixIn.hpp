#ifndef PIKE_BLACK_BOX_MODEL_EVALUATOR_RESPONSE_MIXIN_HPP
#define PIKE_BLACK_BOX_MODEL_EVALUATOR_RESPONSE_MIXIN_HPP

#include "Teuchos_ArrayView.hpp"
#include <string>

namespace pike {

  /** \brief MixIn pure virtual class for supporting responses. */
  class ResponseMixIn {

  public:

    virtual ~ResponseMixIn() {}

    //! Returns true if the response is supported by this model evaluator.
    virtual bool supportsResponse(const std::string& rName) const = 0;

    //! Get the number of responses, Ng.
    virtual int getNumberOfResponses() const = 0;

    //! Returns the response name for index j where 0 <= j < Ng.
    virtual std::string getResponseName(const int j) const = 0;

    //! Returns the response index for the string name.
    virtual int getResponseIndex(const std::string& rName) const = 0;

    //! Returns the response for index j where 0 <= j < Ng. 
    virtual Teuchos::ArrayView<const double> getResponse(const int j) const = 0;

  };

}

#endif
