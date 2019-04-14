// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_ODEParameters_hpp
#define Tempus_ODEParameters_hpp


#include "Tempus_TimeDerivative.hpp"
#include "Tempus_WrapperModelEvaluator.hpp"


namespace Tempus {


template<class Scalar>
class ODEParameters
{
  public:
    /// Constructor
    ODEParameters()
      : timeDer_(Teuchos::null), timeStepSize_(Scalar(0.0)), stageNumber_(0),
        alpha_(Scalar(0.0)), beta_(Scalar(0.0)),
        solutionHistory_(Teuchos::null), evaluationType_(SOLVE_FOR_X)
    {}
    /// Constructor
    ODEParameters(Teuchos::RCP<TimeDerivative<Scalar> > timeDer,
                  Scalar timeStepSize, Scalar alpha, Scalar beta,
                  const Teuchos::RCP<SolutionHistory<Scalar> >& sh,
                  EVALUATION_TYPE evaluationType = SOLVE_FOR_X)
      : timeDer_(timeDer), timeStepSize_(timeStepSize), stageNumber_(0),
        alpha_(alpha), beta_(beta), solutionHistory_(sh),
        evaluationType_(evaluationType)
    {}

    Teuchos::RCP<TimeDerivative<Scalar> > timeDer_;
    Scalar                                timeStepSize_;
    int                                   stageNumber_;
    Scalar                                alpha_;  // Only valid for Implicit
    Scalar                                beta_;   // Only valid for Implicit
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory_;
    EVALUATION_TYPE                       evaluationType_;
};

} // namespace Tempus
#endif // Tempus_ODEParameters_hpp
