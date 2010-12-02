// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
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
