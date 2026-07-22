// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Range1D.hpp"

namespace Teuchos {

const Range1D Range1D::Invalid(Range1D::INVALID);

} // end namespace Teuchos

std::ostream& Teuchos::operator<<(std::ostream &out, const Range1D& rng)
{
  out << "Range1D{";
  if (rng == Range1D::Invalid) {
    out << "Invalid";
  }
  else {
    out << rng.lbound() << "," << rng.ubound();
  }
  out << "}";
  return out;
}
