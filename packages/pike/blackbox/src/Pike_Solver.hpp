#ifndef PIKE_SOLVER_HPP
#define PIKE_SOLVER_HPP

#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"

namespace pike {

  class Solver : public Teuchos::ParameterListAcceptorDefaultBase,
		 public Teuchos::Describable,
		 public Teuchos::VerboseObject<pike::Solver> {

  public:

    virtual ~Solver() {}

    virtual void step() = 0;

    virtual void solve() = 0;

    virtual int getNumberOfIterations() const = 0;

    //@{ \name Methods to override from Teuchos::ParameterListAcceptor.

    virtual void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList) = 0;

    virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const = 0;

    //@}

  };

}

#endif
