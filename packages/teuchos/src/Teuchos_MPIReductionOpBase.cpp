// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_MPIReductionOpBase.hpp"

#ifdef HAVE_MPI

//
// This implementation of handling the reduction operator
// will work just fine in a single threaded program.
// For multi-threaded program this has to be reworked!
//

//
// Management of the reduction operator object
//

namespace {

Teuchos::RefCountPtr<const Teuchos::MPIReductionOpBase> the_reduct_op = Teuchos::null;

} // namespace

extern "C" {

void Teuchos_MPI_reduction_op(
  void              *invec
  ,void             *inoutvec
  ,int              *len
  ,MPI_Datatype    *datatype
  )
{
  Teuchos::get_reduct_op()->reduce(invec,inoutvec,len,datatype);
}

} // extern "C"

void Teuchos::set_reduct_op( const RefCountPtr<const MPIReductionOpBase>& reduct_op )
{
  TEST_FOR_EXCEPT( get_reduct_op() != null && reduct_op != null  );
  the_reduct_op = reduct_op;
}

Teuchos::RefCountPtr<const Teuchos::MPIReductionOpBase>
Teuchos::get_reduct_op()
{
  return the_reduct_op;
}

namespace Teuchos {

//
// Defintion of MPIReductionOpCreator
//

MPIReductionOpCreator::MPIReductionOpCreator( const Teuchos::RefCountPtr<const MPIReductionOpBase>& reduct_op )
{
  set_reduct_op(reduct_op);
  mpi_op_ = MPI_OP_NULL;
  TEST_FOR_EXCEPT(
    0!=MPI_Op_create( &Teuchos_MPI_reduction_op ,1 ,&mpi_op_ ) // Assume op is commutative?
    );
}

MPIReductionOpCreator::~MPIReductionOpCreator()
{
  MPI_Op_free(&mpi_op_);
  set_reduct_op( null );
}

} // namespace Teuchos

#endif // HAVE_MPI
