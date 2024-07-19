// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_EReductionType.hpp"
#include "Teuchos_TestForException.hpp"

namespace Teuchos {

#ifdef HAVE_TEUCHOS_MPI
namespace Details {

MPI_Op
getMpiOpForEReductionType (const enum EReductionType reductionType)
{
  switch (reductionType) {
  case REDUCE_SUM: return MPI_SUM;
  case REDUCE_MIN: return MPI_MIN;
  case REDUCE_MAX: return MPI_MAX;
  case REDUCE_AND: return MPI_LAND; // logical AND, not bitwise AND
  case REDUCE_BOR: return MPI_BOR; // bitwise OR
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      "The given EReductionType value is invalid.");
  }
}

} // namespace Details
#endif // HAVE_TEUCHOS_MPI

const char*
toString (const EReductionType reductType)
{
  switch (reductType) {
  case REDUCE_SUM: return "REDUCE_SUM";
  case REDUCE_MIN: return "REDUCE_MIN";
  case REDUCE_MAX: return "REDUCE_MAX";
  case REDUCE_AND: return "REDUCE_AND";
  case REDUCE_BOR: return "REDUCE_BOR";
  default:
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::invalid_argument, "Teuchos::toString(EReductionType): "
       "Invalid EReductionType value " << reductType << ".  Valid values "
       "are REDUCE_SUM = " << REDUCE_SUM << ", REDUCE_MIN = " << REDUCE_MIN
       << ", REDUCE_MAX = " << REDUCE_MIN << ", REDUCE_AND = " << REDUCE_AND
       << ", and REDUCE_BOR = " << REDUCE_BOR
       << ".");
  }
}

} // namespace Teuchos
