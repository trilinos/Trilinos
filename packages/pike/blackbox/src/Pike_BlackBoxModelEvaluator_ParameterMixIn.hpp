#ifndef PIKE_BLACK_BOX_MODEL_EVALUATOR_PARAMETER_MIXIN_HPP
#define PIKE_BLACK_BOX_MODEL_EVALUATOR_PARAMETER_MIXIN_HPP

#include "Teuchos_ArrayView.hpp"
#include <string>

namespace pike {

  /** \brief MixIn pure virtual class for supporting responses. */
  class ParameterMixIn {

  public:

    virtual ~ParameterMixIn() {}

    //! Returns true if the response is supported by this model evaluator.
    virtual bool supportsParameter(const std::string& name) const = 0;

    //! Get the number of parameters, Np.
    virtual int getNumberOfParameters() const = 0;

    //! Returns the parameter name for index l where 0 <= l < Np.
    virtual std::string getParameterName(const int l) const = 0;

    //! Returns the parameter index for the string name.
    virtual int getParameterIndex(const std::string& name) const = 0;

    //! Sets the parameter, p, for index l where 0 <= l < Np. 
    virtual void setParameter(const int l, const Teuchos::ArrayView<const double>& p) = 0;

  };

}

#endif
