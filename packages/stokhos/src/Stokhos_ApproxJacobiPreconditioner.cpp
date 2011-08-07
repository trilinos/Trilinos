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

#include "Stokhos_ApproxJacobiPreconditioner.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Stokhos_SGOperatorFactory.hpp"

Stokhos::ApproxJacobiPreconditioner::
ApproxJacobiPreconditioner(
  const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk_,
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  const Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>& prec_factory_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  label("Stokhos Approximate Jacobi Preconditioner"),
  sg_comm(sg_comm_),
  sg_basis(sg_basis_),
  epetraCijk(epetraCijk_),
  base_map(base_map_),
  sg_map(sg_map_),
  prec_factory(prec_factory_),
  mean_prec(),
  useTranspose(false),
  num_iter(2),
  sg_op(),
  sg_poly(),
  rhs_block()
{
  num_iter = params_->get("Number of Jacobi Iterations", 2);

  Teuchos::RCP<Teuchos::ParameterList> sgOpParams =
    Teuchos::rcp(&(params_->sublist("Jacobi SG Operator")), false);
  sgOpParams->set("Include Mean", false);
  if (!sgOpParams->isParameter("Only Use Linear Terms"))
    sgOpParams->set("Only Use Linear Terms", true);

  // Build new parallel Cijk if we are only using the linear terms, Cijk
  // is distributed over proc's, and Cijk includes more than just the linear
  // terms (so we have the right column map; otherwise we will be importing
  // much more than necessary)
  if (sgOpParams->get<bool>("Only Use Linear Terms") && 
      epetraCijk->isStochasticParallel()) {
    int dim = sg_basis->dimension();
    if (epetraCijk->getKEnd() > dim+1)
      epetraCijk = 
	Teuchos::rcp(new EpetraSparse3Tensor(*epetraCijk, 1, dim+1));
					     
  }
  Stokhos::SGOperatorFactory sg_op_factory(sgOpParams);
  mat_free_op = sg_op_factory.build(sg_comm, sg_basis, epetraCijk, 
				    base_map, base_map, sg_map, sg_map);
}

Stokhos::ApproxJacobiPreconditioner::
~ApproxJacobiPreconditioner()
{
}

void
Stokhos::ApproxJacobiPreconditioner::
setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op_, 
		    const Epetra_Vector& x)
{
  sg_op = sg_op_;
  sg_poly = sg_op->getSGPolynomial();
  mean_prec = prec_factory->compute(sg_poly->getCoeffPtr(0));
  label = std::string("Stokhos Approximate Jacobi Preconditioner:\n") + 
    std::string("		***** ") + 
    std::string(mean_prec->Label());
  mat_free_op->setupOperator(sg_poly);
}

int 
Stokhos::ApproxJacobiPreconditioner::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;

  return 0;
}

int 
Stokhos::ApproxJacobiPreconditioner::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  return sg_op->Apply(Input, Result);
}

int 
Stokhos::ApproxJacobiPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Approximate Jacobi Time");
#endif

  // We have to be careful if Input and Result are the same vector.
  // If this is the case, the only possible solution is to make a copy
  const Epetra_MultiVector *input = &Input;
  bool made_copy = false;
  if (Input.Values() == Result.Values()) {
    input = new Epetra_MultiVector(Input);
    made_copy = true;
  } 

  int m = input->NumVectors();
  if (rhs_block == Teuchos::null || rhs_block->NumVectors() != m)
    rhs_block = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));

  // Extract blocks
  EpetraExt::BlockMultiVector input_block(View, *base_map, *input);
  EpetraExt::BlockMultiVector result_block(View, *base_map, Result);

  int myBlockRows = epetraCijk->numMyRows();
  result_block.PutScalar(0.0);
  for (int iter=0; iter<num_iter; iter++) {

    // Compute RHS
    if (iter == 0)
      rhs_block->Update(1.0, input_block, 0.0);
    else {
      mat_free_op->Apply(result_block, *rhs_block);
      rhs_block->Update(1.0, input_block, -1.0);
    }

    // Apply deterministic preconditioner
    for(int i=0; i<myBlockRows; i++) {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
      TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total AJ Deterministic Preconditioner Time");
#endif
      mean_prec->ApplyInverse(*(rhs_block->GetBlock(i)),
			      *(result_block.GetBlock(i)));
    }

  }

  if (made_copy)
    delete input;

  return 0; 
}

double 
Stokhos::ApproxJacobiPreconditioner::
NormInf() const
{
  return sg_op->NormInf();
}


const char* 
Stokhos::ApproxJacobiPreconditioner::
Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::ApproxJacobiPreconditioner::
UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::ApproxJacobiPreconditioner::
HasNormInf() const
{
  return sg_op->HasNormInf();
}

const Epetra_Comm & 
Stokhos::ApproxJacobiPreconditioner::
Comm() const
{
  return *sg_comm;
}
const Epetra_Map& 
Stokhos::ApproxJacobiPreconditioner::
OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::ApproxJacobiPreconditioner::
OperatorRangeMap() const
{
  return *sg_map;
}
