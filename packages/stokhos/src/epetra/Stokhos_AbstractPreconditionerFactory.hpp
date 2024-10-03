// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_ABSTRACT_PRECONDITIONER_FACTORY_HPP
#define STOKHOS_ABSTRACT_PRECONDITIONER_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Epetra_Operator.h"

namespace Stokhos {

  //! An abstract class to represent a generic preconditioner factory.
  class AbstractPreconditionerFactory {
  public:

    //! Constructor
    AbstractPreconditionerFactory() {}

    //! Destructor
    virtual ~AbstractPreconditionerFactory() {}

    //! Compute preconditioner operator
    virtual Teuchos::RCP<Epetra_Operator> 
    compute(const Teuchos::RCP<Epetra_Operator>& mat,
	    bool compute_prec = true) = 0;

    //! Recompute preconditioner operator for a new matrix
    virtual void
    recompute(const Teuchos::RCP<Epetra_Operator>& mat,
	      const Teuchos::RCP<Epetra_Operator>& prec) = 0;

  }; // class PreconditionerFactory

} // namespace Stokhos

#endif // STOKHOS_ABSTRACT_PRECONDITIONER_FACTORY_HPP
