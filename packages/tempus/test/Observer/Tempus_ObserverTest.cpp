// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorObserverLogging.hpp"
#include "Tempus_IntegratorObserverComposite.hpp"

#include "../TestModels/SinCosModel.hpp"
#include "../TestModels/VanDerPolModel.hpp"
#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <vector>

namespace Tempus_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

using Tempus::IntegratorBasic;
using Tempus::SolutionHistory;
using Tempus::SolutionState;


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(Observer, IntegratorObserverLogging)
{
  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_Observer_SinCos.xml");

  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  RCP<SinCosModel<double> > model =
    Teuchos::rcp(new SinCosModel<double> (scm_pl));

  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  RCP<Tempus::IntegratorBasic<double> > integrator =
    Tempus::integratorBasic<double>(pl, model);

  RCP<Tempus::SolutionHistory<double> > sh  = integrator->getSolutionHistory();
  RCP<Tempus::TimeStepControl<double> > tsc = integrator->getTimeStepControl();
  RCP<Tempus::IntegratorObserverLogging<double> > loggingObs =
    Teuchos::rcp(new Tempus::IntegratorObserverLogging<double>(sh,tsc));
  integrator->setObserver(loggingObs);

  // Integrate to timeMax
  bool integratorStatus = integrator->advanceTime();
  TEST_ASSERT(integratorStatus)

  // Test if at 'Final Time'
  double time = integrator->getTime();
  double timeFinal = pl->sublist("Demo Integrator")
    .sublist("Time Step Control").get<double>("Final Time");
  TEST_FLOATING_EQUALITY(time, timeFinal, 1.0e-14);

  // Construct the reference counter and order for comparison.
  std::map<std::string,int> refCounters;
  std::list<std::string> refOrder;

  refCounters[loggingObs->nameObserveStartIntegrator_ ] = 1;
  refCounters[loggingObs->nameObserveStartTimeStep_   ] = 10;
  refCounters[loggingObs->nameObserveNextTimeStep_    ] = 10;
  refCounters[loggingObs->nameObserveBeforeTakeStep_  ] = 10;
  refCounters[loggingObs->nameObserveAfterTakeStep_   ] = 10;
  refCounters[loggingObs->nameObserveAcceptedTimeStep_] = 10;
  refCounters[loggingObs->nameObserveEndIntegrator_   ] = 1;

  refOrder.push_back(loggingObs->nameObserveStartIntegrator_ );
  for (int i=0 ; i<10; ++i) {
    refOrder.push_back(loggingObs->nameObserveStartTimeStep_   );
    refOrder.push_back(loggingObs->nameObserveNextTimeStep_    );
    refOrder.push_back(loggingObs->nameObserveBeforeTakeStep_  );
    refOrder.push_back(loggingObs->nameObserveAfterTakeStep_   );
    refOrder.push_back(loggingObs->nameObserveAcceptedTimeStep_);
  }
  refOrder.push_back(loggingObs->nameObserveEndIntegrator_   );

  const std::map<std::string,int>& counters = *(loggingObs->getCounters());
  const std::list<std::string>&    order    = *(loggingObs->getOrder());

  // Compare against reference.
  TEST_EQUALITY(
       counters.find(loggingObs->nameObserveStartIntegrator_ )->second,
    refCounters.find(loggingObs->nameObserveStartIntegrator_ )->second);
  TEST_EQUALITY(
       counters.find(loggingObs->nameObserveStartTimeStep_   )->second,
    refCounters.find(loggingObs->nameObserveStartTimeStep_   )->second);
  TEST_EQUALITY(
       counters.find(loggingObs->nameObserveNextTimeStep_    )->second,
    refCounters.find(loggingObs->nameObserveNextTimeStep_    )->second);
  TEST_EQUALITY(
       counters.find(loggingObs->nameObserveBeforeTakeStep_  )->second,
    refCounters.find(loggingObs->nameObserveBeforeTakeStep_  )->second);
  TEST_EQUALITY(
       counters.find(loggingObs->nameObserveAfterTakeStep_   )->second,
    refCounters.find(loggingObs->nameObserveAfterTakeStep_   )->second);
  TEST_EQUALITY(
       counters.find(loggingObs->nameObserveAcceptedTimeStep_)->second,
    refCounters.find(loggingObs->nameObserveAcceptedTimeStep_)->second);
  TEST_EQUALITY(
       counters.find(loggingObs->nameObserveEndIntegrator_   )->second,
    refCounters.find(loggingObs->nameObserveEndIntegrator_   )->second);

  TEUCHOS_ASSERT(order.size() == refOrder.size());
  std::list<std::string>::const_iterator orderIter = order.begin();
  std::list<std::string>::const_iterator refOrderIter = refOrder.begin();
  for ( ; orderIter != order.end(); ++orderIter,++refOrderIter) {
    //std::cout << *orderIter << std::endl;
    TEST_EQUALITY(*orderIter, *refOrderIter);
  }

  // Test the reset.
  loggingObs->resetLogCounters();
  TEST_EQUALITY(
    counters.find(loggingObs->nameObserveStartIntegrator_ )->second, 0);
  TEST_EQUALITY(
    counters.find(loggingObs->nameObserveStartTimeStep_   )->second, 0);
  TEST_EQUALITY(
    counters.find(loggingObs->nameObserveNextTimeStep_    )->second, 0);
  TEST_EQUALITY(
    counters.find(loggingObs->nameObserveBeforeTakeStep_  )->second, 0);
  TEST_EQUALITY(
    counters.find(loggingObs->nameObserveAfterTakeStep_   )->second, 0);
  TEST_EQUALITY(
    counters.find(loggingObs->nameObserveAcceptedTimeStep_)->second, 0);
  TEST_EQUALITY(
    counters.find(loggingObs->nameObserveEndIntegrator_   )->second, 0);
  TEST_EQUALITY(order.size(), 0);

  Teuchos::TimeMonitor::summarize();
}

/* This class is used to test that the composite observer sets the
   solution history and time step control on all the underlying
   observers */
template<typename Scalar>
class MockObserver : public Tempus::IntegratorObserver<Scalar> {

  Teuchos::RCP<Tempus::SolutionHistory<Scalar> > sh_;
  Teuchos::RCP<Tempus::TimeStepControl<Scalar> > tsc_;

public:
  void observeStartIntegrator() 
  {
    TEUCHOS_ASSERT(nonnull(sh_));
    TEUCHOS_ASSERT(nonnull(tsc_));
  }

  void observeStartTimeStep() override
  {
    TEUCHOS_ASSERT(nonnull(sh_));
    TEUCHOS_ASSERT(nonnull(tsc_));
  }

  void observeNextTimeStep(Tempus::Status & integratorStatus) override
  {
    TEUCHOS_ASSERT(nonnull(sh_));
    TEUCHOS_ASSERT(nonnull(tsc_));
  }

  void observeBeforeTakeStep() override
  {
    TEUCHOS_ASSERT(nonnull(sh_));
    TEUCHOS_ASSERT(nonnull(tsc_));
  }

  void observeAfterTakeStep() override
  {
    TEUCHOS_ASSERT(nonnull(sh_));
    TEUCHOS_ASSERT(nonnull(tsc_));
  }

  void observeAcceptedTimeStep(Tempus::Status & integratorStatus) override
  {
    TEUCHOS_ASSERT(nonnull(sh_));
    TEUCHOS_ASSERT(nonnull(tsc_));
  }

  void observeEndIntegrator(const Tempus::Status integratorStatus) override
  {
    TEUCHOS_ASSERT(nonnull(sh_));
    TEUCHOS_ASSERT(nonnull(tsc_));
  }

  void setSolutionHistory(Teuchos::RCP<Tempus::SolutionHistory<Scalar> > sh) override
  { sh_ = sh; }

  void setTimeStepControl(Teuchos::RCP<Tempus::TimeStepControl<Scalar> > tsc) override
  { tsc_ = tsc; }
};

TEUCHOS_UNIT_TEST( Observer, IntegratorObserverComposite) {

  // Read params from .xml file
  RCP<ParameterList> pList =
    getParametersFromXmlFile("Tempus_Observer_SinCos.xml");
  //
  // Setup the SinCosModel
  RCP<ParameterList> scm_pl = sublist(pList, "SinCosModel", true);
  RCP<SinCosModel<double> > model =
    Teuchos::rcp(new SinCosModel<double> (scm_pl));

  Tempus::Status status = Tempus::Status::PASSED;

  // Setup the Integrator and reset initial time step
  RCP<ParameterList> pl = sublist(pList, "Tempus", true);
  RCP<Tempus::IntegratorBasic<double> > integrator =
    Tempus::integratorBasic<double>(pl, model);

  RCP<Tempus::SolutionHistory<double> > sh  = integrator->getSolutionHistory();
  RCP<Tempus::TimeStepControl<double> > tsc = integrator->getTimeStepControl();

  RCP<Tempus::IntegratorObserverLogging<double> > loggingObs =
    Teuchos::rcp(new Tempus::IntegratorObserverLogging<double>(sh,tsc));

  // creating another logging observer
  RCP<Tempus::IntegratorObserverLogging<double> > loggingObs2 =
    Teuchos::rcp(new Tempus::IntegratorObserverLogging<double>(sh,tsc));

  // Create the Mock observer. This will test that the sh and tsc
  // parameters are set on the sub-observers correctly. We have cases
  // where the observers are created without access to sh and tsc and
  // need to be registered after time integrator construction.
  RCP<Tempus::IntegratorObserver<double>> mockObs = rcp(new MockObserver<double>);

  RCP<Tempus::IntegratorObserverComposite<double> > compObs = 
      Teuchos::rcp(new Tempus::IntegratorObserverComposite<double>(sh, tsc));

  compObs->addObserver(loggingObs);
  compObs->addObserver(loggingObs2);
  compObs->addObserver(mockObs);

  compObs->observeStartIntegrator();
  compObs->observeStartTimeStep();
  compObs->observeBeforeTakeStep();
  compObs->observeAfterTakeStep();
  compObs->observeAcceptedTimeStep(status);


  for (int i=0 ; i<10; ++i) {
      compObs->observeStartTimeStep();
      compObs->observeBeforeTakeStep();
      compObs->observeAfterTakeStep();
      compObs->observeAcceptedTimeStep(status);
  }
  compObs->observeEndIntegrator(status);

 const std::map<std::string,int>& counters = *(loggingObs->getCounters());

  TEST_EQUALITY(
        counters.find(loggingObs->nameObserveStartIntegrator_  )->second, 1);
  TEST_EQUALITY(
        counters.find(loggingObs->nameObserveStartTimeStep_    )->second,11);
  TEST_EQUALITY(
        counters.find(loggingObs->nameObserveBeforeTakeStep_   )->second,11);
  TEST_EQUALITY(
        counters.find(loggingObs->nameObserveAfterTakeStep_    )->second,11);
  TEST_EQUALITY(
        counters.find(loggingObs->nameObserveAcceptedTimeStep_ )->second,11);
  TEST_EQUALITY(
        counters.find(loggingObs->nameObserveEndIntegrator_    )->second, 1);
  TEST_EQUALITY(
        counters.find(loggingObs2->nameObserveStartIntegrator_ )->second, 1);
  TEST_EQUALITY(
        counters.find(loggingObs2->nameObserveStartTimeStep_   )->second,11);
  TEST_EQUALITY(
        counters.find(loggingObs2->nameObserveBeforeTakeStep_  )->second,11);
  TEST_EQUALITY(
        counters.find(loggingObs2->nameObserveAfterTakeStep_   )->second,11);
  TEST_EQUALITY(
        counters.find(loggingObs2->nameObserveAcceptedTimeStep_)->second,11);
  TEST_EQUALITY(
        counters.find(loggingObs2->nameObserveEndIntegrator_   )->second, 1);

  // Test that setSolutionHistory and setTimeStepControl pass on new
  // values to underlying objects
  compObs->setSolutionHistory(Teuchos::null); // break history
  TEST_THROW(compObs->observeStartIntegrator(),std::logic_error);
  TEST_THROW(compObs->observeStartTimeStep(),std::logic_error);
  TEST_THROW(compObs->observeBeforeTakeStep(),std::logic_error);
  TEST_THROW(compObs->observeAfterTakeStep(),std::logic_error);
  TEST_THROW(compObs->observeAcceptedTimeStep(status),std::logic_error);
  TEST_THROW(compObs->observeEndIntegrator(status),std::logic_error);

  compObs->setSolutionHistory(sh); // fix history
  compObs->setTimeStepControl(Teuchos::null); // break ts control
  TEST_THROW(compObs->observeStartIntegrator(),std::logic_error);
  TEST_THROW(compObs->observeStartTimeStep(),std::logic_error);
  TEST_THROW(compObs->observeBeforeTakeStep(),std::logic_error);
  TEST_THROW(compObs->observeAfterTakeStep(),std::logic_error);
  TEST_THROW(compObs->observeAcceptedTimeStep(status),std::logic_error);
  TEST_THROW(compObs->observeEndIntegrator(status),std::logic_error);
  
  compObs->setTimeStepControl(tsc); // fix ts control
  TEST_NOTHROW(compObs->observeStartIntegrator());
  TEST_NOTHROW(compObs->observeStartTimeStep());
  TEST_NOTHROW(compObs->observeBeforeTakeStep());
  TEST_NOTHROW(compObs->observeAfterTakeStep());
  TEST_NOTHROW(compObs->observeAcceptedTimeStep(status));
  TEST_NOTHROW(compObs->observeEndIntegrator(status));
}


} // namespace Tempus_Test
