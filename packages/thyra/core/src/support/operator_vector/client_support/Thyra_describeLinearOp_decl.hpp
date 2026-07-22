// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DESCRIBE_LINEAR_OP_DECL_HPP
#define THYRA_DESCRIBE_LINEAR_OP_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerbosityLevel.hpp"


namespace Thyra {


/** \brief Basic implementation of a describe function for a linear operator.
 *
 * \ingroup Thyra_Op_Vec_ANA_Developmnet_support_code_grp
 */
template<class Scalar>
void describeLinearOp(
  const LinearOpBase<Scalar> &A,
  Teuchos::FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  );


}	// end namespace Thyra


#endif // THYRA_DESCRIBE_LINEAR_OP_DECL_HPP
