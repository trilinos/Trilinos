// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_LINEAR_OP_DEFAULT_BASE_HPP
#define THYRA_LINEAR_OP_DEFAULT_BASE_HPP

#include "Thyra_LinearOpDefaultBase_decl.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_describeLinearOp.hpp"

namespace Thyra {


// Overridden from Teuchos::Describable


template<class Scalar>
std::string LinearOpDefaultBase<Scalar>::description() const
{
  std::ostringstream oss;
  const Teuchos::RCP<const VectorSpaceBase<Scalar> >
    l_range = this->range();
  const Teuchos::RCP<const VectorSpaceBase<Scalar> >
    l_domain = this->domain();
  oss << Teuchos::Describable::description();
  if(!l_range.get()) {
    oss << "{range=NULL,domain=NULL}"; 
  }
  else {
    const Ordinal dimDomain = l_domain->dim(), dimRange = l_range->dim();
    oss
      << "{rangeDim=" << dimRange
      << ",domainDim=" << dimDomain << "}";
  }
  return oss.str();
}


template<class Scalar>
void LinearOpDefaultBase<Scalar>::describe(
  Teuchos::FancyOStream &out,
  const Teuchos::EVerbosityLevel verbLevel
  ) const
{
  describeLinearOp(*this, out, verbLevel);
}


}	// end namespace Thyra

#endif // THYRA_LINEAR_OP_DEFAULT_BASE_HPP
