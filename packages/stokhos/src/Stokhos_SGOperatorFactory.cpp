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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
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
