// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Stokhos_SGOperatorFactory.hpp"
#include "Stokhos_MatrixFreeOperator.hpp"
#include "Stokhos_KLMatrixFreeOperator.hpp"
#include "Stokhos_KLReducedMatrixFreeOperator.hpp"
#include "Stokhos_FullyAssembledOperator.hpp"
#include "Teuchos_Assert.hpp"

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
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
		       "Error!  Unknown operator method " << op_method
		       << "." << std::endl);

  return sg_op;
}
