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

#include "Stokhos_MPPreconditionerFactory.hpp"
#include "Stokhos_MPBlockDiagonalPreconditioner.hpp"
#include "Stokhos_MPMeanBasedPreconditioner.hpp"
#include "Stokhos_PreconditionerFactory.hpp"
#include "Teuchos_TestForException.hpp"

Stokhos::MPPreconditionerFactory::
MPPreconditionerFactory(const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  params(params_)
{
}

Teuchos::RCP<Stokhos::MPPreconditioner> 
Stokhos::MPPreconditionerFactory::
build(const Teuchos::RCP<const EpetraExt::MultiComm>& mp_comm,
      int num_mp_blocks,
      const Teuchos::RCP<const Epetra_Map>& base_map,
      const Teuchos::RCP<const Epetra_Map>& mp_map)
{
  Teuchos::RCP<Stokhos::MPPreconditioner> mp_prec;
  std::string prec_method = params->get("Preconditioner Method", 
					"Block Diagonal");
  if (prec_method == "Block Diagonal") {
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory = 
      buildPointPreconditionerFactory();
    mp_prec = Teuchos::rcp(new Stokhos::MPBlockDiagonalPreconditioner(
			     mp_comm, num_mp_blocks, 
			     base_map, mp_map, prec_factory, 
			     params));
  }
  else if (prec_method == "Mean-based") {
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory = 
      buildPointPreconditionerFactory();
    mp_prec = Teuchos::rcp(new Stokhos::MPMeanBasedPreconditioner(
			     mp_comm, num_mp_blocks, 
			     base_map, mp_map, prec_factory, 
			     params));
  }
  else if (prec_method == "None")
    mp_prec = Teuchos::null;
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Unknown preconditioner method " << prec_method
		       << "." << std::endl);

  return mp_prec;
}

Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>
Stokhos::MPPreconditionerFactory::
buildPointPreconditionerFactory()
{
  std::string prec_name = 
    params->get("MP Preconditioner Type", "Ifpack");
  Teuchos::RCP<Teuchos::ParameterList> precParams = 
    Teuchos::rcp(&params->sublist("MP Preconditioner Parameters"),false);
  return Teuchos::rcp(new Stokhos::PreconditionerFactory(prec_name, 
							 precParams));
}
