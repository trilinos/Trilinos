// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SG_OPERATOR_FACTORY_HPP
#define STOKHOS_SG_OPERATOR_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_SGOperator.hpp"
#include "EpetraExt_MultiComm.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Epetra_Map.h"

namespace Stokhos {

  //! Factory for generating stochastic Galerkin preconditioners
  class SGOperatorFactory {
  public:

    //! Constructor
    SGOperatorFactory(
      const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Destructor
    virtual ~SGOperatorFactory() {}

    //! Build preconditioner operator
    virtual Teuchos::RCP<Stokhos::SGOperator> 
    build(
      const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk,
      const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map,
      const Teuchos::RCP<const Epetra_Map>& domain_sg_map,
      const Teuchos::RCP<const Epetra_Map>& range_sg_map);

  private:
    
    //! Private to prohibit copying
    SGOperatorFactory(const SGOperatorFactory&);
    
    //! Private to prohibit copying
    SGOperatorFactory& operator=(const SGOperatorFactory&);

  protected:

    //! Operator parameters
    Teuchos::RCP<Teuchos::ParameterList> params;

  }; // class SGOperatorFactory

} // namespace Stokhos

#endif // STOKHOS_SG_OPERATOR_FACTORY_HPP
