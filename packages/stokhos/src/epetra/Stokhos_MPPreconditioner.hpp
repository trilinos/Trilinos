// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MP_PRECONDITIONER_HPP
#define STOKHOS_MP_PRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Epetra_Operator.h"
#include "Stokhos_BlockDiagonalOperator.hpp"
#include "Epetra_Vector.h"

namespace Stokhos {

  /*! 
   * \brief An abstract class to represent a generic stochastic Galerkin 
   * preconditioner as an Epetra_Operator.
   */
  class MPPreconditioner : public virtual Epetra_Operator {
  public:

    //! Constructor
    MPPreconditioner() {}

    //! Destructor
    virtual ~MPPreconditioner() {}

    //! Setup preconditioner
    virtual void 
    setupPreconditioner(
      const Teuchos::RCP<Stokhos::BlockDiagonalOperator>& mp_op, 
      const Epetra_Vector& x) = 0;

  }; // class MPPreconditioner

} // namespace Stokhos

#endif // STOKHOS_MP_PRECONDITIONER_HPP
