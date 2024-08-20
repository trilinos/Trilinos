// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SG_PRECONDITIONER_HPP
#define STOKHOS_SG_PRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Epetra_Operator.h"
#include "Stokhos_SGOperator.hpp"
#include "Epetra_Vector.h"

namespace Stokhos {

  /*! 
   * \brief An abstract class to represent a generic stochastic Galerkin 
   * preconditioner as an Epetra_Operator.
   */
  class SGPreconditioner : public virtual Epetra_Operator {
  public:

    //! Constructor
    SGPreconditioner() {}

    //! Destructor
    virtual ~SGPreconditioner() {}

    //! Setup preconditioner
    virtual void 
    setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op, 
			const Epetra_Vector& x) = 0;

  }; // class SGPreconditioner

} // namespace Stokhos

#endif // STOKHOS_SG_PRECONDITIONER_HPP
