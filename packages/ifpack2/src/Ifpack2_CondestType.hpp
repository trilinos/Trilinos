// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_CONDESTTYPE_HPP
#define IFPACK2_CONDESTTYPE_HPP

namespace Ifpack2 {

//! Ifpack2::CondestType: enum to define the type of condition number estimate.

enum CondestType {
  Cheap,  //!< cheap estimate
  CG,     //!< Uses AztecOO's CG
  GMRES   //!< Uses AztecOO's GMRES
};

}//namespace Ifpack2

#endif
