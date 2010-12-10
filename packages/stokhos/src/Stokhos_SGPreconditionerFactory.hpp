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

#ifndef STOKHOS_SG_PRECONDITIONER_FACTORY_HPP
#define STOKHOS_SG_PRECONDITIONER_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_SGPreconditioner.hpp"
#include "EpetraExt_MultiComm.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Epetra_Map.h"
#include "Stokhos_PreconditionerFactory.hpp"

namespace Stokhos {

  //! Factory for generating stochastic Galerkin preconditioners
  class SGPreconditionerFactory {
  public:

    //! Constructor
    SGPreconditionerFactory(
      const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Destructor
    virtual ~SGPreconditionerFactory() {}

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
    Teuchos::RCP<Stokhos::PreconditionerFactory> 
    buildMeanPreconditionerFactory();

    //! Build preconditioner factory
    Teuchos::RCP<Stokhos::PreconditionerFactory> 
    buildPreconditionerFactory(
      const std::string& prec_name,
      const Teuchos::RCP<Teuchos::ParameterList>& precParams);

  private:
    
    //! Private to prohibit copying
    SGPreconditionerFactory(const SGPreconditionerFactory&);
    
    //! Private to prohibit copying
    SGPreconditionerFactory& operator=(const SGPreconditionerFactory&);

  protected:

    //! Preconditioner parameters
    Teuchos::RCP<Teuchos::ParameterList> params;

  }; // class SGPreconditionerFactory

} // namespace Stokhos

#endif // STOKHOS_SG_PRECONDITIONER_FACTORY_HPP
