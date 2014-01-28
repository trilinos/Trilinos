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

#include "Stokhos_FullyAssembledOperator.hpp"
#include "Teuchos_TimeMonitor.hpp"

Stokhos::FullyAssembledOperator::
FullyAssembledOperator(
  const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk_,
  const Teuchos::RCP<const Epetra_CrsGraph>& base_graph,
  const Teuchos::RCP<const Epetra_Map>& domain_sg_map_,
  const Teuchos::RCP<const Epetra_Map>& range_sg_map_,
  const Teuchos::RCP<Teuchos::ParameterList>& params) : 
  EpetraExt::BlockCrsMatrix(*base_graph, 
			    *(epetraCijk_->getStochasticGraph()), 
			    *sg_comm_),
  sg_comm(sg_comm_),
  sg_basis(sg_basis_),
  epetraCijk(epetraCijk_),
  domain_sg_map(domain_sg_map_),
  range_sg_map(range_sg_map_),
  Cijk(epetraCijk->getParallelCijk()),
  block_ops(),
  scale_op(true),
  include_mean(true),
  only_use_linear(false)
{
  scale_op = params->get("Scale Operator by Inverse Basis Norms", true);
  include_mean = params->get("Include Mean", true);
  only_use_linear = params->get("Only Use Linear Terms", false);
}

Stokhos::FullyAssembledOperator::
~FullyAssembledOperator()
{
}

void 
Stokhos::FullyAssembledOperator::
setupOperator(
   const Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly >& ops)
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: SG Fully Assembled Operator Assembly");
#endif

  block_ops = ops;

  // Zero out matrix
  this->PutScalar(0.0);

  // Compute loop bounds
  Cijk_type::k_iterator k_begin = Cijk->k_begin();
  Cijk_type::k_iterator k_end = Cijk->k_end();
  if (!include_mean && index(k_begin) == 0)
    ++k_begin;
  if (only_use_linear) {
    int dim = sg_basis->dimension();
    k_end = Cijk->find_k(dim+1);
  }

  // Assemble matrix
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
    int k = index(k_it);
    Teuchos::RCP<Epetra_RowMatrix> block = 
      Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(block_ops->getCoeffPtr(k), 
						  true);
    for (Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
	 j_it != Cijk->j_end(k_it); ++j_it) {
      int j = epetraCijk->GCID(index(j_it));
      for (Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	   i_it != Cijk->i_end(j_it); ++i_it) {
	int i = epetraCijk->GRID(index(i_it));
	double c = value(i_it);
	if (scale_op)
	  c /= norms[i];
	this->SumIntoGlobalBlock(c, *block, i, j);
      }
    }
  }

  this->FillComplete(*domain_sg_map, *range_sg_map);
}

Teuchos::RCP< Stokhos::EpetraOperatorOrthogPoly > 
Stokhos::FullyAssembledOperator::
getSGPolynomial()
{
  return block_ops;
}

Teuchos::RCP<const Stokhos::EpetraOperatorOrthogPoly > 
Stokhos::FullyAssembledOperator::
getSGPolynomial() const
{
  return block_ops;
}

int 
Stokhos::FullyAssembledOperator::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: SG Operator Apply");
#endif
  return Epetra_CrsMatrix::Apply(Input, Result);
}
