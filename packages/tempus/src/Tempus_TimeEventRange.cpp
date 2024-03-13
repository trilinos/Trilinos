// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "Tempus_TimeEventRange.hpp"
#include "Tempus_TimeEventRange_impl.hpp"

namespace Tempus {

TEMPUS_INSTANTIATE_TEMPLATE_CLASS(TimeEventRange)

// Nonmember constructor
template Teuchos::RCP<TimeEventRange<double> > createTimeEventRange(
    Teuchos::RCP<Teuchos::ParameterList> pl);

}  // namespace Tempus

#endif
