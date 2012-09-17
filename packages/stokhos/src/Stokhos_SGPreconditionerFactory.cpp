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

#include "Stokhos_ConfigDefs.h"
#include "Stokhos_SGPreconditionerFactory.hpp"
#include "Stokhos_MeanBasedPreconditioner.hpp"
#include "Stokhos_ApproxGaussSeidelPreconditioner.hpp"
#include "Stokhos_ApproxJacobiPreconditioner.hpp"
#include "Stokhos_ApproxSchurComplementPreconditioner.hpp"
#include "Stokhos_GaussSeidelPreconditioner.hpp"
#include "Stokhos_KroneckerProductPreconditioner.hpp"
#include "Stokhos_FullyAssembledPreconditioner.hpp"
#include "Stokhos_PreconditionerFactory.hpp"
#include "Teuchos_Assert.hpp"

Stokhos::SGPreconditionerFactory::
SGPreconditionerFactory(const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  params(params_)
{
  prec_method = params->get("Preconditioner Method", "Mean-based");
}

bool
Stokhos::SGPreconditionerFactory::
isPrecSupported() const
{
  return prec_method != "None";
}

Teuchos::RCP<Stokhos::SGPreconditioner> 
Stokhos::SGPreconditionerFactory::
build(const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk,
      const Teuchos::RCP<const Epetra_Map>& base_map,
      const Teuchos::RCP<const Epetra_Map>& sg_map)
{
  Teuchos::RCP<Stokhos::SGPreconditioner> sg_prec;
  
  if (prec_method == "Mean-based") {
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory = 
      buildMeanPreconditionerFactory();
    sg_prec = Teuchos::rcp(new Stokhos::MeanBasedPreconditioner(
			     sg_comm, sg_basis, epetraCijk, 
			     base_map, sg_map, prec_factory, 
			     params));
  }
  else if (prec_method == "Approximate Gauss-Seidel") {
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory = 
      buildMeanPreconditionerFactory();
    sg_prec = Teuchos::rcp(new Stokhos::ApproxGaussSeidelPreconditioner(
			     sg_comm, sg_basis, epetraCijk, 
			     base_map, sg_map, prec_factory, 
			     params));
  }
  else if (prec_method == "Approximate Jacobi") {
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory = 
      buildMeanPreconditionerFactory();
    sg_prec = Teuchos::rcp(new Stokhos::ApproxJacobiPreconditioner(
			     sg_comm, sg_basis, epetraCijk, 
			     base_map, sg_map, prec_factory, 
			     params));
  }
  else if (prec_method == "Approximate Schur Complement") {
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory = 
      buildMeanPreconditionerFactory();
    sg_prec = Teuchos::rcp(new Stokhos::ApproxSchurComplementPreconditioner(
			     sg_comm, sg_basis, epetraCijk, 
			     base_map, sg_map, prec_factory, 
			     params));
  }
#ifdef HAVE_STOKHOS_NOX
  else if (prec_method == "Gauss-Seidel") {
    Teuchos::RCP<NOX::Epetra::LinearSystem> det_solver = 
      params->get< Teuchos::RCP<NOX::Epetra::LinearSystem> >("Deterministic Solver");
    sg_prec = Teuchos::rcp(new Stokhos::GaussSeidelPreconditioner(
			     sg_comm, sg_basis, epetraCijk, 
			     base_map, sg_map, det_solver, 
			     params));
  }
#endif
  else if (prec_method == "Kronecker Product") {
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> mean_prec_factory = 
      buildMeanPreconditionerFactory();
    std::string G_prec_name = 
      params->get("G Preconditioner Type", "Ifpack");
    Teuchos::RCP<Teuchos::ParameterList> G_prec_params =
      Teuchos::rcp(&(params->sublist("G Preconditioner Parameters")),false);
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> G_prec_factory = 
      Teuchos::rcp(new Stokhos::PreconditionerFactory(G_prec_name, 
						      G_prec_params));
    sg_prec = Teuchos::rcp(new Stokhos::KroneckerProductPreconditioner(
			     sg_comm, sg_basis, epetraCijk, base_map, sg_map, 
			     mean_prec_factory, G_prec_factory, params));
  }
  else if (prec_method == "Fully Assembled") {
     std::string prec_name = 
       params->get("Fully Assembled Preconditioner Type", "Ifpack");
    Teuchos::RCP<Teuchos::ParameterList> prec_params =
      Teuchos::rcp(&(params->sublist("Fully Assembled Preconditioner Parameters")),false);
    Teuchos::RCP<Stokhos::AbstractPreconditionerFactory> prec_factory = 
      Teuchos::rcp(new Stokhos::PreconditionerFactory(prec_name, prec_params));
    sg_prec = Teuchos::rcp(new Stokhos::FullyAssembledPreconditioner(
			     prec_factory, params));
  }
  else if (prec_method == "None")
    sg_prec = Teuchos::null;
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Unknown preconditioner method " << prec_method
		       << "." << std::endl);

  return sg_prec;
}

Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>
Stokhos::SGPreconditionerFactory::
buildMeanPreconditionerFactory()
{
  std::string prec_name = 
    params->get("Mean Preconditioner Type", "Ifpack");
  Teuchos::RCP<Teuchos::ParameterList> precParams = 
    Teuchos::rcp(&params->sublist("Mean Preconditioner Parameters"),false);
  return Teuchos::rcp(new Stokhos::PreconditionerFactory(prec_name, 
							 precParams));
}
