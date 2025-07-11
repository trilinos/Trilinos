//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorFactory_impl_hpp
#define Tempus_PhiEvaluatorFactory_impl_hpp

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluatorPFD.hpp"
#include "Tempus_PhiEvaluatorLeja.hpp"
//#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Tempus_PhiEvaluatorFactory.hpp"

namespace Tempus {

template <class Scalar>
Teuchos::RCP<PhiEvaluator<Scalar> > PhiEvaluatorFactory<Scalar>::createPhiEvaluator(
    std::string phiEvaluatorType,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  if (phiEvaluatorType == "") phiEvaluatorType = "PFD";
  return this->createPhiEvaluator(phiEvaluatorType, Teuchos::null, model);
}

template <class Scalar>
Teuchos::RCP<PhiEvaluator<Scalar> > PhiEvaluatorFactory<Scalar>::createPhiEvaluator(
    Teuchos::RCP<Teuchos::ParameterList> phiEvaluatorPL,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  std::string phiEvaluatorType = "PFD";
  if (phiEvaluatorPL != Teuchos::null)
    phiEvaluatorType = phiEvaluatorPL->get<std::string>("PhiEvaluator Type", "PFD");
  return this->createPhiEvaluator(phiEvaluatorType, phiEvaluatorPL, model);
}

template <class Scalar>
Teuchos::RCP<PhiEvaluator<Scalar> > PhiEvaluatorFactory<Scalar>::createPhiEvaluator(
    std::string phiEvaluatorType, Teuchos::RCP<Teuchos::ParameterList> phiEvaluatorPL,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  if (phiEvaluatorType == "PFD") {
    return createPhiEvaluatorPFD<Scalar>(phiEvaluatorPL);
  }
  else if (phiEvaluatorType == "Leja") {
    return createPhiEvaluatorLeja<Scalar>(phiEvaluatorPL);
  }
  else {
    Teuchos::RCP<Teuchos::FancyOStream> out =
        Teuchos::VerboseObjectBase::getDefaultOStream();
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 1, "PhiFactoryFactory::createPhiEvaluator");
    *out << "Unknown PhiEvaluator Type!  ('" + phiEvaluatorType + "').\n"
         << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Unknown 'PhiEvaluator Type' = " << phiEvaluatorType);
  }
}

}  // namespace Tempus
#endif  // Tempus_PhiEvaluatorFactory_impl_hpp
