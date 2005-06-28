#ifndef Capo_INTEGRATOR_H
#define Capo_INTEGRATOR_H

#include "Thyra_VectorBase.hpp"
//#include "Thyra_EpetraThyraWrappers.hpp"

namespace CAPO {

  class Integrator {
  public:
    Integrator() {};
    
    virtual ~Integrator() {};
    
    //! Perform picard iteration: given current guess y, param value param,
    //! and previously-converged solution in ynew (needed for converging
    //! implicit time steps), return next iterate ynew.
    
    
    virtual bool Integrate(Teuchos::RefCountPtr<Thyra::VectorBase<double> >& y, Teuchos::RefCountPtr<Thyra::VectorBase<double> >& ynew, const double T, const double lambda) = 0;
  };

}

#endif
