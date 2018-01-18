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

  // non-member ctor
  template Teuchos::RCP<IntegratorBasic<double> > integratorBasic(
    Teuchos::RCP<Teuchos::ParameterList>        parameterList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model);

  // non-member ctor
  template Teuchos::RCP<IntegratorBasic<double> > integratorBasic(
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    std::string stepperType);

  // non-member ctor
  template Teuchos::RCP<IntegratorBasic<double> > integratorBasic();

  // non-member ctor
  template Teuchos::RCP<IntegratorBasic<double> > integratorBasic(
    Teuchos::RCP<Teuchos::ParameterList>                     pList,
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<double> > > models);

} // namespace Tempus

#endif
