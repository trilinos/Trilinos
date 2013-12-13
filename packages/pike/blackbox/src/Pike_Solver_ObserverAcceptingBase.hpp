#ifndef PIKE_SOLVER_OBSERVER_ACCEPTING_BASE_HPP
#define PIKE_SOLVER_OBSERVER_ACCEPTING_BASE_HPP

namespace pike {

  class ObserverAcceptingBase {

    virtual ~ObserverAcceptingBase() {}

    void addObserver(const Teuchos::RCP<pike::Observer>& observer) = 0;

    void getObserver(const std::string& name) const = 0;

  };

}

#endif
