// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_ML_PRECONDITIONER_FACTORY_HPP
#define STOKHOS_ML_PRECONDITIONER_FACTORY_HPP

#include "Teuchos_ParameterList.hpp"

#include "Stokhos_AbstractPreconditionerFactory.hpp"

namespace Stokhos {

  //! A factory for building ML preconditioners
  class MLPreconditionerFactory : 
    public Stokhos::AbstractPreconditionerFactory {
  public:

    //! Constructor
    MLPreconditionerFactory(const Teuchos::RCP<Teuchos::ParameterList>& p);

    //! Destructor
    virtual ~MLPreconditionerFactory() {}

    //! Compute preconditioner
    virtual Teuchos::RCP<Epetra_Operator> 
    compute(const Teuchos::RCP<Epetra_Operator>& op,
	    bool compute_prec = true);

    //! Recompute preconditioner operator for a new matrix
    virtual void
    recompute(const Teuchos::RCP<Epetra_Operator>& op,
	      const Teuchos::RCP<Epetra_Operator>& prec);

  protected:

    //! Preconditioner parameters
    Teuchos::RCP<Teuchos::ParameterList> precParams;

  }; // class MLPreconditionerFactory

} // namespace Stokhos

#endif // STOKHOS_ML_PRECONDITIONER_FACTORY_HPP
