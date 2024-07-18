// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MP_PRECONDITIONER_FACTORY_HPP
#define STOKHOS_MP_PRECONDITIONER_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_MPPreconditioner.hpp"
#include "EpetraExt_MultiComm.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Epetra_Map.h"
#include "Stokhos_AbstractPreconditionerFactory.hpp"

namespace Stokhos {

  //! Factory for generating stochastic Galerkin preconditioners
  class MPPreconditionerFactory {
  public:

    //! Constructor
    MPPreconditionerFactory(
      const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Destructor
    virtual ~MPPreconditionerFactory() {}

    //! Build preconditioner operator
    virtual Teuchos::RCP<Stokhos::MPPreconditioner> 
    build(
      const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm,
      int num_mp_blocks,
      const Teuchos::RCP<const Epetra_Map>& base_map,
      const Teuchos::RCP<const Epetra_Map>& mp_map);

  protected:

    //! Build preconditioner factory for each point
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> 
    buildPointPreconditionerFactory();

  private:
    
    //! Private to prohibit copying
    MPPreconditionerFactory(const MPPreconditionerFactory&);
    
    //! Private to prohibit copying
    MPPreconditionerFactory& operator=(const MPPreconditionerFactory&);

  protected:

    //! Preconditioner parameters
    Teuchos::RCP<Teuchos::ParameterList> params;

  }; // class MPPreconditionerFactory

} // namespace Stokhos

#endif // STOKHOS_MP_PRECONDITIONER_FACTORY_HPP
