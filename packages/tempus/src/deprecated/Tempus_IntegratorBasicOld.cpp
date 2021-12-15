// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifdef __GNUC__
#  warning "This file Tempus_IntegratorBasicOld.cpp is deprecated!  Use Tempus_IntegratorBasic.cpp instead!"
#endif

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_IntegratorBasicOld.hpp"
#include "Tempus_IntegratorBasicOld_impl.hpp"

namespace Tempus {

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(IntegratorBasicOld)

  // Nonmember ctor
  template Teuchos::RCP<IntegratorBasicOld<double> > integratorBasic(
    Teuchos::RCP<Teuchos::ParameterList>        parameterList,
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model);

  // Nonmember ctor
  template Teuchos::RCP<IntegratorBasicOld<double> > integratorBasic(
    const Teuchos::RCP<Thyra::ModelEvaluator<double> >& model,
    std::string stepperType);

  // Nonmember ctor
  template Teuchos::RCP<IntegratorBasicOld<double> > integratorBasic();

  // Nonmember ctor
  template Teuchos::RCP<IntegratorBasicOld<double> > integratorBasic(
    Teuchos::RCP<Teuchos::ParameterList>                     pList,
    std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<double> > > models);

} // namespace Tempus

#endif
