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

#include "Stokhos_SGOperatorFactory.hpp"
#include "Stokhos_MatrixFreeOperator.hpp"
#include "Stokhos_KLMatrixFreeOperator.hpp"
#include "Stokhos_KLReducedMatrixFreeOperator.hpp"
#include "Stokhos_FullyAssembledOperator.hpp"
#include "Teuchos_TestForException.hpp"

Stokhos::SGOperatorFactory::
SGOperatorFactory(const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  params(params_)
{
}

Teuchos::RCP<Stokhos::SGOperator> 
Stokhos::SGOperatorFactory::
build(const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
      const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk,
      const Teuchos::RCP<const Epetra_Map>& domain_base_map,
      const Teuchos::RCP<const Epetra_Map>& range_base_map,
      const Teuchos::RCP<const Epetra_Map>& domain_sg_map,
      const Teuchos::RCP<const Epetra_Map>& range_sg_map)
{
  Teuchos::RCP<Stokhos::SGOperator> sg_op;
  std::string op_method = params->get("Operator Method", "Matrix Free");
  if (op_method == "Matrix Free") {
    sg_op = Teuchos::rcp(new Stokhos::MatrixFreeOperator(
			   sg_comm, sg_basis, epetraCijk, 
			   domain_base_map, range_base_map, 
			   domain_sg_map, range_sg_map, params));
  }
  else if (op_method == "KL Matrix Free") {
    sg_op = Teuchos::rcp(new Stokhos::KLMatrixFreeOperator(
			   sg_comm, sg_basis, epetraCijk, 
			   domain_base_map, range_base_map, 
			   domain_sg_map, range_sg_map, params));
  }
  else if (op_method == "KL Reduced Matrix Free") {
    sg_op = Teuchos::rcp(new Stokhos::KLReducedMatrixFreeOperator(
			   sg_comm, sg_basis, epetraCijk, 
			   domain_base_map, range_base_map, 
			   domain_sg_map, range_sg_map, params));
  }
  else if (op_method == "Fully Assembled") {
    Teuchos::RCP<const Epetra_CrsGraph> base_graph = 
      params->get< Teuchos::RCP<const Epetra_CrsGraph> >("Base Graph");
    sg_op = Teuchos::rcp(new Stokhos::FullyAssembledOperator(
			   sg_comm, sg_basis, epetraCijk, base_graph, 
			   domain_sg_map, range_sg_map, params));
  }
  else
    TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Unknown operator method " << op_method
		       << "." << std::endl);

  return sg_op;
}
