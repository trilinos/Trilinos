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

#include "Stokhos_MeanBasedPreconditioner.hpp"
#include "EpetraExt_BlockMultiVector.h"

Stokhos::MeanBasedPreconditioner::
MeanBasedPreconditioner(
  const Teuchos::RCP<const EpetraExt::MultiComm>& sg_comm_,
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
  const Teuchos::RCP<const Stokhos::EpetraSparse3Tensor>& epetraCijk_,
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  const Teuchos::RCP<Stokhos::AbstractPreconditionerFactory>& prec_factory_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  label("Stokhos Mean-Based Preconditioner"),
  sg_comm(sg_comm_),
  sg_basis(sg_basis_),
  epetraCijk(epetraCijk_),
  base_map(base_map_),
  sg_map(sg_map_),
  num_blocks(0),
  prec_factory(prec_factory_),
  mean_prec(),
  useTranspose(false)
{
}

Stokhos::MeanBasedPreconditioner::
~MeanBasedPreconditioner()
{
}

void
Stokhos::MeanBasedPreconditioner::
setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op, 
		    const Epetra_Vector& x)
{
   TEST_FOR_EXCEPTION(prec_factory == Teuchos::null, std::logic_error,
		      "Error!  setupPreconditioner() cannot be called when " <<
		      "prec_factory is null!" << std::endl);

   Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly > sg_poly = 
     sg_op->getSGPolynomial();
   mean_prec = prec_factory->compute(sg_poly->getCoeffPtr(0));
   label = std::string("Stokhos Mean-Based Preconditioner:\n") + 
     std::string("		***** ") + 
     std::string(mean_prec->Label());
   num_blocks = sg_basis()->size();
}

int 
Stokhos::MeanBasedPreconditioner::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  mean_prec->SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::MeanBasedPreconditioner::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  int myBlockRows = epetraCijk->numMyRows();
  EpetraExt::BlockMultiVector sg_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);
  for (int i=0; i<myBlockRows; i++) {
    mean_prec->Apply(*(sg_input.GetBlock(i)), *(sg_result.GetBlock(i)));
  }

  return 0;
}

int 
Stokhos::MeanBasedPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  int myBlockRows = epetraCijk->numMyRows();
  EpetraExt::BlockMultiVector sg_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);
  for (int i=0; i<myBlockRows; i++) {
    mean_prec->ApplyInverse(*(sg_input.GetBlock(i)), *(sg_result.GetBlock(i)));
  }

  return 0;
}

double 
Stokhos::MeanBasedPreconditioner::
NormInf() const
{
  return mean_prec->NormInf();
}


const char* 
Stokhos::MeanBasedPreconditioner::
Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::MeanBasedPreconditioner::
UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::MeanBasedPreconditioner::
HasNormInf() const
{
  return true;
}

const Epetra_Comm & 
Stokhos::MeanBasedPreconditioner::
Comm() const
{
  return *sg_comm;
}
const Epetra_Map& 
Stokhos::MeanBasedPreconditioner::
OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::MeanBasedPreconditioner::
OperatorRangeMap() const
{
  return *sg_map;
}
