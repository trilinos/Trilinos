// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_SPMD_apply_op_decl.hpp"
#include "RTOpPack_SPMD_apply_op_def.hpp"


Teuchos::RCP<Teuchos::FancyOStream>& RTOpPack::spmdApplyOpDumpOut()
{
  static Teuchos::RCP<Teuchos::FancyOStream> dumpOut;
  return dumpOut;
}


void RTOpPack::set_SPMD_apply_op_dump_out(const RCP<FancyOStream> &dumpOut)
{
  spmdApplyOpDumpOut() = dumpOut;
}


#ifdef HAVE_RTOP_EXPLICIT_INSTANTIATION


#include "Teuchos_ExplicitInstantiationHelpers.hpp"


namespace RTOpPack {


TEUCHOS_MACRO_TEMPLATE_INSTANT_SCALAR_TYPES(
  RTOPPACK_SPMD_APPLY_OP_INSTANT_SCALAR)


} // namespace RTOpPack


#endif // HAVE_TEUCHOS_EXCPLICIT_INSTANTIATION
