// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_MPPreconditionerFactory.hpp"
#include "Stokhos_MPBlockDiagonalPreconditioner.hpp"
#include "Stokhos_MPMeanBasedPreconditioner.hpp"
#include "Stokhos_PreconditionerFactory.hpp"
#include "Teuchos_Assert.hpp"

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
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
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
