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

#include "Stokhos_KroneckerProductPreconditioner.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Epetra_LocalMap.h"
#include "EpetraExt_BlockMultiVector.h"

Stokhos::KroneckerProductPreconditioner::
KroneckerProductPreconditioner(
  const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk_,
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  const Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>& mean_prec_factory_,
  const Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>& G_prec_factory_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  label("Stokhos Kronecker Product Preconditioner"),
  sg_comm(sg_comm_),
  sg_basis(sg_basis_),
  epetraCijk(epetraCijk_),
  base_map(base_map_),
  sg_map(sg_map_),
  mean_prec_factory(mean_prec_factory_),
  G_prec_factory(G_prec_factory_),
  params(params_),
  mean_prec(),
  useTranspose(false),
  sg_op(),
  sg_poly(),
  Cijk(epetraCijk->getParallelCijk()),
  scale_op(true),
  only_use_linear(false)
{
  scale_op = params->get("Scale Operator by Inverse Basis Norms", true);
  only_use_linear = params_->get("Only Use Linear Terms", false);
  // Build new parallel Cijk if we are only using the linear terms, Cijk
  // is distributed over proc's, and Cijk includes more than just the linear
  // terms (so we have the right column map; otherwise we will be importing
  // much more than necessary)
  if (only_use_linear && epetraCijk->isStochasticParallel()) {
    int dim = sg_basis->dimension();
    if (epetraCijk->getKEnd() > dim+1)
      epetraCijk = 
	Teuchos::rcp(new EpetraSparse3Tensor(*epetraCijk, 1, dim+1));
					     
  }
}

Stokhos::KroneckerProductPreconditioner::
~KroneckerProductPreconditioner()
{
}

void
Stokhos::KroneckerProductPreconditioner::
setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op_, 
		    const Epetra_Vector& x)
{
  sg_op = sg_op_;
  sg_poly = sg_op->getSGPolynomial();
  mean_prec = mean_prec_factory->compute(sg_poly->getCoeffPtr(0));
  label = std::string("Stokhos Kronecker Product Preconditioner:\n") + 
    std::string("		***** ") + 
    std::string(mean_prec->Label());
  
  // Build graph of G matrix
  Teuchos::RCP<const Epetra_CrsGraph> graph = epetraCijk->getStochasticGraph();

  // Construct G matrix: G_{ij} = \sum tr(A'B)/tr(A'A)*<psi_alpha,psi_i,psi_j>.
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();
  G = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
  Teuchos::RCP<Epetra_CrsMatrix> A0 =
    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(sg_poly->getCoeffPtr(0), true);
  double traceAB0 = MatrixTrace(*A0, *A0);
  Cijk_type::k_iterator k_begin = Cijk->k_begin();
  Cijk_type::k_iterator k_end = Cijk->k_end();
  if (only_use_linear) {
    int dim = sg_basis->dimension();
    k_end = Cijk->find_k(dim+1);
  }
  for (Cijk_type::k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
    int k = index(k_it);
    Teuchos::RCP<Epetra_CrsMatrix> A_k =
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(sg_poly->getCoeffPtr(k), 
						  true);
    double traceAB = MatrixTrace(*A_k, *A0);
    for (Cijk_type::kj_iterator j_it = Cijk->j_begin(k_it); 
	 j_it != Cijk->j_end(k_it); ++j_it) {
      int j = epetraCijk->GCID(index(j_it));
      for (Cijk_type::kji_iterator i_it = Cijk->i_begin(j_it);
	   i_it != Cijk->i_end(j_it); ++i_it) {
        int i = epetraCijk->GRID(index(i_it));
	double c = value(i_it)*traceAB/traceAB0;
	if (scale_op)
	  c /= norms[i];
	G->SumIntoGlobalValues(i, 1, &c, &j);
      }
    }
  }
  G->FillComplete();

  // Build G preconditioner
  G_prec = G_prec_factory->compute(G);

  label = std::string("Stokhos Kronecker Product Preconditioner:\n") + 
    std::string("		***** ") + 
    std::string(mean_prec->Label()) + std::string("\n") + 
    std::string("		***** ") + 
    std::string(G_prec->Label());
}

int 
Stokhos::KroneckerProductPreconditioner::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  mean_prec->SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::KroneckerProductPreconditioner::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  EpetraExt::BlockMultiVector sg_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);

  const Epetra_Map& G_map = G->RowMap();
  int NumMyElements = G_map.NumMyElements();
  int vecLen = sg_input.GetBlock(0)->MyLength(); // Global length of vector.
  int m = sg_input.NumVectors();

  if (result_MVT == Teuchos::null || result_MVT->NumVectors() != vecLen*m) {
    result_MVT = Teuchos::rcp(new Epetra_MultiVector(G_map, vecLen*m));
  }

  // Preconditioner is P = (G x I)(I x A_0)

  // Apply I x A_0
  for (int i=0; i<NumMyElements; i++) {
    mean_prec->Apply(*(sg_input.GetBlock(i)), *(sg_result.GetBlock(i)));
  }

  Teuchos::RCP<Epetra_MultiVector> x;
  for (int irow=0 ; irow<NumMyElements; irow++) {
    x = sg_result.GetBlock(irow);
    for (int vcol=0; vcol<m; vcol++) {
      for (int icol=0; icol<vecLen; icol++) {
        (*result_MVT)[m*vcol+icol][irow] = (*x)[vcol][icol];
      }
    }
  }

  // Apply G x I
  G_prec->Apply(*result_MVT, *result_MVT);

  for (int irow=0; irow<NumMyElements; irow++) {
    x = sg_result.GetBlock(irow);
    for (int vcol=0; vcol<m; vcol++) {
      for (int icol=0; icol<vecLen; icol++) {
        (*x)[vcol][icol] = (*result_MVT)[m*vcol+icol][irow];
      }
    }
  }

  return 0;
}

int 
Stokhos::KroneckerProductPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
  TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Kronecker Product Prec Time");
#endif

  EpetraExt::BlockMultiVector sg_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);

  const Epetra_Map& G_map = G->RowMap();
  int NumMyElements = G_map.NumMyElements();
  int vecLen = sg_input.GetBlock(0)->MyLength(); // Global length of vector.
  int m = sg_input.NumVectors();

  if (result_MVT == Teuchos::null || result_MVT->NumVectors() != vecLen*m) {
    result_MVT = Teuchos::rcp(new Epetra_MultiVector(G_map, vecLen*m));
  }

  // Preconditioner is P^{-1} = (I x A_0^{-1})(G^{-1} x I)
  // => y = P^{-1}x => Y = A_0^{-1}(G^{-1}X^T)^T where X = multivec(x)

  Teuchos::RCP<Epetra_MultiVector> x;
  for (int irow=0 ; irow<NumMyElements; irow++) {
    x = sg_input.GetBlock(irow);
    for (int vcol=0; vcol<m; vcol++) {
      for (int icol=0; icol<vecLen; icol++) {
	(*result_MVT)[m*vcol+icol][irow] = (*x)[vcol][icol];
      }
    }
  }

  // Apply (G^{-1} x I)
  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
    TEUCHOS_FUNC_TIME_MONITOR("Stokhos: G Preconditioner Apply Inverse");
#endif
    G_prec->ApplyInverse(*result_MVT, *result_MVT);
  }
  
  for (int irow=0; irow<NumMyElements; irow++) {
    x = sg_result.GetBlock(irow);
    for (int vcol=0; vcol<m; vcol++) {
      for (int icol=0; icol<vecLen; icol++) {
	(*x)[vcol][icol] = (*result_MVT)[m*vcol+icol][irow];
      }
    }
  }

  // Apply (I x A_0^{-1})
  {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
    TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Mean Preconditioner Apply Inverse");
#endif
    for (int i=0; i<NumMyElements; i++) {
      mean_prec->ApplyInverse(*(sg_result.GetBlock(i)), 
			      *(sg_result.GetBlock(i)));
    }
  }
  
  return 0;
}

double 
Stokhos::KroneckerProductPreconditioner::
NormInf() const
{
  // I think this is only an upper bound
  return mean_prec->NormInf() * G_prec->NormInf();
}


const char* 
Stokhos::KroneckerProductPreconditioner::
Label() const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::KroneckerProductPreconditioner::
UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::KroneckerProductPreconditioner::
HasNormInf() const
{
  return mean_prec->NormInf() && G_prec->NormInf();
}

const Epetra_Comm & 
Stokhos::KroneckerProductPreconditioner::
Comm() const
{
  return *sg_comm;
}
const Epetra_Map& 
Stokhos::KroneckerProductPreconditioner::
OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::KroneckerProductPreconditioner::
OperatorRangeMap() const
{
  return *sg_map;
}

double 
Stokhos::KroneckerProductPreconditioner::
MatrixTrace(const Epetra_CrsMatrix& A, const Epetra_CrsMatrix& B) const {
  int n = A.NumMyRows(); // # of rows on the processor.
  double traceAB = 0.0;
  for (int i=0; i<n; i++) {
    int m = A.NumMyEntries(i); // # of non zero entries in row i.
    for (int j=0; j<m; j++) { 
      traceAB += A[i][j]*B[i][j];
    }
  }

  return traceAB;
}

