// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperObserverComposite_decl_hpp
#define Tempus_StepperObserverComposite_decl_hpp

#include "Tempus_StepperObserver.hpp"
#include "Tempus_TimeStepControl.hpp"
#include <vector>

namespace Tempus {

/** \brief This observer is a composite observer,
 *
 *  which takes other StepperObservers and sequentially calls each
 *  individual observer function.
 */
template<class Scalar>
class StepperObserverComposite
  : virtual public Tempus::StepperObserver<Scalar>
{
public:

  /// Default constructor
  StepperObserverComposite();

  /// Destructor
  virtual ~StepperObserverComposite();

  /// \name Override StepperObserver basic methods
  //@{
    /// Observe the beginning of the time integrator.
    virtual void observeBeginTakeStep(Teuchos::RCP<SolutionHistory<Scalar> > sh, Stepper<Scalar> & stepper) override;

    /// Observe the beginning of the time step loop.
    virtual void observeEndTakeStep(Teuchos::RCP<SolutionHistory<Scalar> > sh, Stepper<Scalar> & stepper) override;

    // add observer to the composite observer list
    void addObserver(const Teuchos::RCP<StepperObserver<Scalar> > &observer);

    // clear all observer from the composite observer list
    void clearObservers();

    // get the number of RK stepper observers present in the composite
    std::size_t getSize() const { return observers_.size(); }
  //@}

private:

  std::vector<Teuchos::RCP<StepperObserver<Scalar > > > observers_;

};

} // namespace Tempus
#endif // Tempus_StepperObserverComposite_decl_hpp
