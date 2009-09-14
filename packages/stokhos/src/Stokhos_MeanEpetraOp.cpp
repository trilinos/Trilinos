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

#include "Epetra_config.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_MeanEpetraOp.hpp"

Stokhos::MeanEpetraOp::MeanEpetraOp(
   const Teuchos::RCP<const Epetra_Map>& base_map_,
   const Teuchos::RCP<const Epetra_Map>& sg_map_,
   unsigned int num_blocks_,
   const Teuchos::RCP<Epetra_Operator>& mean_op_) :
  label("Stokhos::MeanEpetraOp"),
  base_map(base_map_),
  sg_map(sg_map_),
  mean_op(mean_op_),
  useTranspose(false),
  num_blocks(num_blocks_)
{
}

Stokhos::MeanEpetraOp::~MeanEpetraOp()
{
}

void
Stokhos::MeanEpetraOp::setMeanOperator(const Teuchos::RCP<Epetra_Operator>& op)
{
  mean_op = op;
}

Teuchos::RCP<Epetra_Operator>
Stokhos::MeanEpetraOp::getMeanOperator()
{
  return mean_op;
}

int 
Stokhos::MeanEpetraOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  mean_op->SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::MeanEpetraOp::Apply(const Epetra_MultiVector& Input, 
			     Epetra_MultiVector& Result) const
{
  EpetraExt::BlockMultiVector sg_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);
  for (unsigned int i=0; i<num_blocks; i++) {
    mean_op->Apply(*(sg_input.GetBlock(i)), *(sg_result.GetBlock(i)));
  }

  return 0;
}

int 
Stokhos::MeanEpetraOp::ApplyInverse(const Epetra_MultiVector& Input, 
				    Epetra_MultiVector& Result) const
{
  EpetraExt::BlockMultiVector sg_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);
  for (unsigned int i=0; i<num_blocks; i++) {
    mean_op->ApplyInverse(*(sg_input.GetBlock(i)), *(sg_result.GetBlock(i)));
  }

  return 0;
}

double 
Stokhos::MeanEpetraOp::NormInf() const
{
  return mean_op->NormInf();
}


const char* 
Stokhos::MeanEpetraOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::MeanEpetraOp::UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::MeanEpetraOp::HasNormInf() const
{
  return true;
}

const Epetra_Comm & 
Stokhos::MeanEpetraOp::Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
Stokhos::MeanEpetraOp::OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::MeanEpetraOp::OperatorRangeMap() const
{
  return *sg_map;
}
