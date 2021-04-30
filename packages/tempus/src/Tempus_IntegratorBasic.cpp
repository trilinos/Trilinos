// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_IntegratorBasic.hpp"
#include "Tempus_IntegratorBasic_impl.hpp"

namespace Tempus {

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(IntegratorBasic)

  // Nonmember ctor
  template Teuchos::RCP<IntegratorBasic<double> > createIntegratorBasic(
    Teuchos::RCP<Teuchos::ParameterList>        parameterList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model);

  // Nonmember ctor
  template Teuchos::RCP<IntegratorBasic<double> > createIntegratorBasic(
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    std::string stepperType);

  // Nonmember ctor
  template Teuchos::RCP<IntegratorBasic<double> > createIntegratorBasic();

  // Nonmember ctor
  template Teuchos::RCP<IntegratorBasic<double> > createIntegratorBasic(
    Teuchos::RCP<Teuchos::ParameterList>                     pList,
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<double> > > models);

} // namespace Tempus

#endif
