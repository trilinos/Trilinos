// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_PRECONDITIONER_FACTORY_HPP
#define STOKHOS_PRECONDITIONER_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Operator.h"
#include "Stokhos_AbstractPreconditionerFactory.hpp"

namespace Stokhos {

  //! An class for building preconditioners
  class PreconditionerFactory : public Stokhos::AbstractPreconditionerFactory {
  public:

    //! Constructor
    PreconditionerFactory(
      const std::string& prec_name,
      const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Destructor
    virtual ~PreconditionerFactory() {}

    //! Compute preconditioner operator
    virtual Teuchos::RCP<Epetra_Operator> 
    compute(const Teuchos::RCP<Epetra_Operator>& mat,
	    bool compute_prec = true);

    //! Recompute preconditioner operator for a new matrix
    virtual void
    recompute(const Teuchos::RCP<Epetra_Operator>& mat,
	      const Teuchos::RCP<Epetra_Operator>& prec);

  protected:

    //! Preconditioner factory
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> factory;

  }; // class PreconditionerFactory

} // namespace Stokhos

#endif // STOKHOS_PRECONDITIONER_FACTORY_HPP
