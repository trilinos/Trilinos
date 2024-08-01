// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_FieldType_hpp__
#define __Panzer_FieldType_hpp__

namespace panzer {

  //! The type of discretization to use for a field pattern
  enum class FieldType
  {
    CG,       //! Continuous Galerkin Formulation
    DG,       //! Discontinuous Galerkin Formulation
  };

}

#endif
