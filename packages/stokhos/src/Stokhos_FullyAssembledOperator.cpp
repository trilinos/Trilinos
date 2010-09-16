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

#include "Stokhos_FullyAssembledOperator.hpp"

Stokhos::FullyAssembledOperator::
FullyAssembledOperator(
  const Teuchos::RCP<const Epetra_CrsMatrix>& base_matrix,
  const Teuchos::RCP<const std::vector< std::vector<int> > >& rowStencil,
  const Teuchos::RCP<const std::vector<int> >& rowIndex,
  const Teuchos::RCP<const Epetra_Comm>& sg_comm,
  const Teuchos::RCP<Teuchos::ParameterList>& params) : 
  EpetraExt::BlockCrsMatrix(*base_matrix, *rowStencil, *rowIndex, *sg_comm),
  Cijk(),
  block_ops(),
  scale_op(true),
  include_mean(true),
  only_use_linear(false)
{
  if (params != Teuchos::null) {
    scale_op = params->get("Scale Operator by Inverse Basis Norms", true);
    include_mean = params->get("Include Mean", true);
    only_use_linear = params->get("Only Use Linear Terms", false);
  }
}

Stokhos::FullyAssembledOperator::
~FullyAssembledOperator()
{
}

void 
Stokhos::FullyAssembledOperator::
setupOperator(
   const Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> >& ops,
   const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_)
{
  block_ops = ops;
  Cijk = Cijk_;

  // Zero out matrix
  this->PutScalar(0.0);

  // Compute loop bounds
  int num_blocks = block_ops->size();
  int k_begin = 0;
  if (!include_mean)
    k_begin = 1;
  int k_end = num_blocks;
  int dim = block_ops->basis()->dimension();
  if (only_use_linear && num_blocks > dim+1)
    k_end = dim + 1;

  // Assemble matrix
  const Teuchos::Array<double>& norms = block_ops->basis()->norm_squared();
  for (int k=k_begin; k<k_end; k++) {
    Teuchos::RCP<Epetra_RowMatrix> block = 
      Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(block_ops->getCoeffPtr(k), 
						  true);
    int nj = Cijk->num_j(k);
    if (nj > 0) {
      const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
      for (int jj=0; jj<nj; jj++) {
	int j = j_indices[jj];
	const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,jj);
	const Teuchos::Array<double>& c_values = Cijk->values(k,jj);
	for (int ii=0; ii<i_indices.size(); ii++) {
	  int i = i_indices[ii];
	  double c = c_values[ii];
	  if (scale_op)
	    c /= norms[i];
	  this->SumIntoBlock(c, *block, i, j);
	}
      }
    }
  }
}

Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Operator> > 
Stokhos::FullyAssembledOperator::
getSGPolynomial()
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::VectorOrthogPoly<Epetra_Operator> > 
Stokhos::FullyAssembledOperator::
getSGPolynomial() const
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> > 
Stokhos::FullyAssembledOperator::
getTripleProduct() const
{
  return Cijk;
}
