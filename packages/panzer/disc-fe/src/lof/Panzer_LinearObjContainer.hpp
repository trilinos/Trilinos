// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_LinearObjContainer_hpp__
#define __Panzer_LinearObjContainer_hpp__

#include "PanzerDiscFE_config.hpp"

#include "Panzer_GlobalEvaluationData.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace panzer {

// **********************************************************************************************
// *********************** LINEAR OBJ CONTAINER *************************************************
// **********************************************************************************************

/**
 * \brief Abstract, linear-algebra-library-agnostic container for the vectors and matrix used in a nonlinear/transient solve.
 *
 * Concrete subclasses (e.g. panzer::TpetraLinearObjContainer) hold the
 * actual solution vector `x`, its time derivative `dxdt`, the residual
 * vector `f`, and the Jacobian matrix, in whatever linear algebra
 * library they wrap, and are produced/consumed by a matching
 * LinearObjFactory. As a GlobalEvaluationData_Default, it participates
 * in the gather/scatter preEvaluate lookup mechanism but performs no
 * ghost/global communication itself (that is handled by the owning
 * LinearObjFactory).
 */
class LinearObjContainer : public GlobalEvaluationData_Default {
public:
   virtual ~LinearObjContainer() {}

   /// \brief Bitmask flags identifying which of the linear objects (solution, time derivative, residual, matrix) are present/requested.
   typedef enum { X=0x1, DxDt=0x2, F=0x4, Mat=0x8} Members;

   /// \brief Zeroes out all linear objects currently held by this container.
   virtual void initialize() = 0;
};

}

#endif
