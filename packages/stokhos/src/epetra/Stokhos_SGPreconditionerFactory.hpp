// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SG_PRECONDITIONER_FACTORY_HPP
#define STOKHOS_SG_PRECONDITIONER_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_SGPreconditioner.hpp"
#include "EpetraExt_MultiComm.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Epetra_Map.h"
#include "Stokhos_AbstractPreconditionerFactory.hpp"

namespace Stokhos {

  //! Factory for generating stochastic Galerkin preconditioners
  class SGPreconditionerFactory {
  public:

    //! Constructor
    SGPreconditionerFactory(
      const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Destructor
    virtual ~SGPreconditionerFactory() {}

    //! Return whether a preconditioner will be supported
    virtual bool isPrecSupported() const;

    //! Build preconditioner operator
    virtual Teuchos::RCP<Stokhos::SGPreconditioner> 
    build(
      const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk,
      const Teuchos::RCP<const Epetra_Map>& base_map,
      const Teuchos::RCP<const Epetra_Map>& sg_map);

  protected:

    //! Build preconditioner factory for mean
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> 
    buildMeanPreconditionerFactory();

  private:
    
    //! Private to prohibit copying
    SGPreconditionerFactory(const SGPreconditionerFactory&);
    
    //! Private to prohibit copying
    SGPreconditionerFactory& operator=(const SGPreconditionerFactory&);

  protected:

    //! Preconditioner parameters
    Teuchos::RCP<Teuchos::ParameterList> params;

    //! Preconditioner method
    std::string prec_method;

  }; // class SGPreconditionerFactory

} // namespace Stokhos

#endif // STOKHOS_SG_PRECONDITIONER_FACTORY_HPP
